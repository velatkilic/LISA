#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lidar Scatterer Augmentation (LISA)
    - Given lidar point cloud and a rain rate generates corresponding noisy signal
    - Reflection data must be normalized to range [0 1] and
      range must be in units of meters
"""
import numpy as np
from scipy.special import gamma
from scipy.integrate import trapz
import PyMieScatt as ps


class LISA():
    def __init__(self,m=1.328,lam=905,rmax=200,rmin=1.5,bdiv=3e-3,dst=0.05,
                 dR=0.09,saved_model=False,atm_model='rain',mode='strongest'):
        '''
        Initialize LISA class
        Parameters
        ----------
        m           : refractive index contrast
        lam         : wavelength (nm)
        rmax        : max lidar range (m)
        rmin        : min lidar range (m)
        bdiv        : beam divergence angle (rad)
        dst         : droplet diameter starting point (mm)
        dR          : range accuracy (m)
        saved_model : use saved mie coefficients (bool)
        atm_model   : atmospheric model type
        mode        : lidar return mode: "strongest" or "last"

        Returns
        -------
        None.

        '''
        self.m    = m
        self.lam  = lam
        self.rmax = rmax   # max range (m)
        self.bdiv = bdiv  # beam divergence (rad)
        self.dst  = dst   # min rain drop diameter to be sampled (mm)
        self.rmin = rmin   # min lidar range (bistatic)
        self.dR   = dR
        self.mode = mode
        self.atm_model = atm_model
        
        
        if saved_model:
            # If Mie parameters are saved, use those
            dat = np.load('mie_q.npz')
            self.D     = dat['D']
            self.qext  = dat['qext']
            self.qback = dat['qback']
        else:
            try:
                dat = np.load('mie_q.npz')
                self.D     = dat['D']
                self.qext  = dat['qext']
                self.qback = dat['qback']
            except:
                # else calculate Mie parameters
                print('Calculating Mie coefficients... \nThis might take a few minutes')
                self.D,self.qext,self.qback = self.calc_Mie_params()
                print('Mie calculation done...')
        
        # Diameter distribution function based on user input
        if atm_model=='rain':
            self.N_model = lambda D, Rr    : self.N_MP_rain(D,Rr)
            self.N_tot   = lambda Rr,dst   : self.N_MP_tot_rain(Rr,dst)
            self.N_sam   = lambda Rr,N,dst : self.MP_Sample_rain(Rr,N,dst)
            
            # Augmenting function: hybrid Monte Carlo
            self.augment  = lambda pc,Rr : self.augment_mc(pc,Rr)
        
        elif atm_model=='snow':
            self.N_model = lambda D, Rr    : self.N_MG_snow(D,Rr)
            self.N_tot   = lambda Rr,dst   : self.N_MG_tot_snow(Rr,dst)
            self.N_sam   = lambda Rr,N,dst : self.MG_Sample_snow(Rr,N,dst)
            self.m       = 1.3031 # refractive index of ice
            
            # Augmenting function: hybrid Monte Carlo
            self.augment  = lambda pc,Rr : self.augment_mc(pc,Rr)
        
        elif atm_model=='chu_hogg_fog':
            self.N_model = lambda D : self.Nd_chu_hogg(D)
            
            # Augmenting function: average effects
            self.augment  = lambda pc : self.augment_avg(pc)
        
        elif atm_model=='strong_advection_fog':
            self.N_model = lambda D : self.Nd_strong_advection_fog(D)
            
            # Augmenting function: average effects
            self.augment  = lambda pc : self.augment_avg(pc)
        
        elif atm_model=='moderate_advection_fog':
            self.N_model = lambda D : self.Nd_moderate_advection_fog(D)
            
            # Augmenting function: average effects
            self.augment  = lambda pc : self.augment_avg(pc)
    
    def augment_mc(self,pc,Rr):
        '''
        Augment clean pointcloud for a given rain rate
        Parameters
        ----------
        pc : pointcloud (N,4) -> x,y,z,reflectivity
        Rr : rain rate (mm/hr)

        Returns
        -------
        pc_new : new noisy point cloud (N,5) -> x,y,z,reflectivity,label
                        label 0 -> lost point
                        label 1 -> randomly scattered point
                        label 2 -> not-scattered 
        '''
        shp    = pc.shape
        pc_new = np.zeros((shp[0],shp[1]+1))
        leng = len(pc)
        for i in range(leng):
            x    = pc[i,0]
            y    = pc[i,1]
            z    = pc[i,2]
            ref  = pc[i,3]
            if ref!=0:
                pc_new[i,:]  = self.lisa_mc(x,y,z,ref,Rr)            
        return pc_new
    
    def lisa_mc(self,x,y,z,ref,Rr):
        '''
        For a single lidar return, performs a hybrid Monte-Carlo experiment

        Parameters
        ----------
        x,y,z : coordinates of the point
        ref   : reflectivity [0 1]
        Rr    : rain rate (mm/hr)

        Returns
        -------
        x,y,z   : new coordinates of the noisy lidar point
        ref_new : new reflectivity
        '''
        rmax = self.rmax                      # max range (m)
        Pmin = 0.9*rmax**(-2)                 # min measurable power (arb units)
        
        bdiv = self.bdiv                      # beam divergence (rad)
        Db   = lambda x: 1e3*np.tan(bdiv)*x   # beam diameter (mm) for a given range (m)
        
        dst  = self.dst                       # min rain drop diameter to be sampled (mm)
        n    = self.m                         # refractive index of scatterer
        rmin = self.rmin                      # min lidar range (bistatic)
        
        
        Nd          = self.N_model(self.D,Rr) # density of rain droplets (m^-3)
        alpha, beta = self.alpha_beta(Nd)     # extinction coeff. (1/m)  
        
        ran   = np.sqrt(x**2 + y**2 + z**2)                               # range in m
        if ran>rmin:
            bvol  = (np.pi/3)*ran*(1e-3*Db(ran)/2)**2                         # beam volume in m^3 (cone)
            Nt    = self.N_tot(Rr,dst) * bvol                                 # total number of particles in beam path
            Nt    = np.int32(np.floor(Nt) + (np.random.rand() < Nt-int(Nt)))  # convert to integer w/ probabilistic rounding
        else:
            Nt = 0
            
        ran_r = ran*(np.random.rand(Nt))**(1/3) # sample distances from a quadratic pdf
        indx  = np.where(ran_r>rmin)[0]         # keep points where ranges larger than rmin
        Nt    = len(indx)                       # new particle number
        
        P0  = ref*np.exp(-2*alpha*ran)/(ran**2) # power
        snr = P0/Pmin # signal noise ratio
        if Nt>0:
            Dr    = self.N_sam(Rr,Nt,dst) # randomly sample Nt particle diameters
            ref_r = abs((n-1)/(n+1))**2   # Fresnel reflection at normal incidence
            ran_r = ran_r[indx]
            
            # Calculate powers for all particles       
            Pr = ref_r*np.exp(-2*alpha*ran_r)*np.minimum((Dr/Db(ran_r))**2,np.ones(Dr.shape))/(ran_r**2)
            if (self.mode=='strongest'):
                ind_r = np.argmax(Pr) # index of the max power
                
                if P0<Pmin and Pr[ind_r]<Pmin: # if all smaller than Pmin, do nothing
                    ran_new = 0
                    ref_new = 0
                    labl    = 0 # label for lost point
                elif P0<Pr[ind_r]: # scatterer has larger power
                    ran_new = ran_r[ind_r] # new range is scatterer range
                    ref_new = ref_r*np.exp(-2*alpha*ran_new)*np.minimum((Dr[ind_r]/Db(ran_r[ind_r]))**2,1) # new reflectance biased by scattering
                    labl    = 1 # label for randomly scattered point 
                else: # object return has larger power
                    sig     = self.dR/np.sqrt(2*snr)        # std of range uncertainty
                    ran_new = ran + np.random.normal(0,sig) # range with uncertainty added
                    ref_new = ref*np.exp(-2*alpha*ran)      # new reflectance modified by scattering
                    labl    = 2                             # label for a non-scattering point
            elif (self.mode=='last'):
                # if object power larger than Pmin, then nothing is scattered
                if P0>Pmin:
                    sig     = self.dR/np.sqrt(2*snr)        # std of range uncertainty
                    ran_new = ran + np.random.normal(0,sig) # range with uncertainty added
                    ref_new = ref*np.exp(-2*alpha*ran)      # new reflectance modified by scattering
                    labl    = 2                             # label for a non-scattering point
                # otherwise find the furthest point above Pmin
                else:
                    inds = np.where(Pr>Pmin)[0]
                    if len(inds) == 0:
                        ran_new = 0
                        ref_new = 0
                        labl    = 0 # label for lost point
                    else:
                        ind_r   = np.where(ran_r == np.max(ran_r[inds]))[0]
                        ran_new = ran_r[ind_r] # new range is scatterer range
                        ref_new = ref_r*np.exp(-2*alpha*ran_new)*np.minimum((Dr[ind_r]/Db(ran_r[ind_r]))**2,1) # new reflectance biased by scattering
                        labl    = 1 # label for randomly scattered point 
                    
            else:
                print("Invalid lidar return mode")
            
        else:
            if P0<Pmin:
                ran_new = 0
                ref_new = 0
                labl    = 0 # label for lost point
            else:
                sig     = self.dR/np.sqrt(2*snr)        # std of range uncertainty
                ran_new = ran + np.random.normal(0,sig) # range with uncertainty added
                ref_new = ref*np.exp(-2*alpha*ran)      # new reflectance modified by scattering
                labl    = 2                             # label for a non-scattering point
        
        # Angles are same
        if ran>0:
            phi = np.arctan2(y,x)  # angle in radians
            the = np.arccos(z/ran) # angle in radians
        else:
            phi,the=0,0
        
        # Update new x,y,z based on new range
        x = ran_new*np.sin(the)*np.cos(phi)
        y = ran_new*np.sin(the)*np.sin(phi)
        z = ran_new*np.cos(the)
        
        return x,y,z,ref_new,labl
    
    def augment_avg(self,pc):

        shp    = pc.shape      # data shape
        pc_new = np.zeros(shp) # init new point cloud
        leng   = shp[0]        # data length
        
        # Rename variables for better readability
        x    = pc[:,0]
        y    = pc[:,1]
        z    = pc[:,2]
        ref  = pc[:,3]          
        
        # Get parameters from class init
        rmax = self.rmax       # max range (m)
        Pmin = 0.9*rmax**(-2)  # min measurable power (arb units)
        rmin = self.rmin       # min lidar range (bistatic)
        
        # Calculate extinction coefficient from the particle distribution
        Nd          = self.N_model(self.D) # density of rain droplets (m^-3)
        alpha, beta = self.alpha_beta(Nd)  # extinction coeff. (1/m)  
        
        ran   = np.sqrt(x**2 + y**2 + z**2)  # range in m
        indx  = np.where(ran>rmin)[0]         # keep points where ranges larger than rmin
        
        P0        = np.zeros((leng,))                                  # init back reflected power
        P0[indx]  = ref[indx]*np.exp(-2*alpha*ran[indx])/(ran[indx]**2) # calculate reflected power
        snr       = P0/Pmin                                             # signal noise ratio
        
        indp = np.where(P0>Pmin)[0] # keep points where power is larger than Pmin
        
        sig        = np.zeros((leng,))                         # init sigma - std of range uncertainty
        sig[indp]  = self.dR/np.sqrt(2*snr[indp])               # calc. std of range uncertainty
        ran_new    = np.zeros((leng,))                         # init new range
        ran_new[indp]    = ran[indp] + np.random.normal(0,sig[indp])  # range with uncertainty added, keep range 0 if P<Pmin
        ref_new    = ref*np.exp(-2*alpha*ran)                   # new reflectance modified by scattering
        
        # Init angles
        phi = np.zeros((leng,))
        the = np.zeros((leng,))
        
        phi[indx] = np.arctan2(y[indx],x[indx])   # angle in radians
        the[indx] = np.arccos(z[indx]/ran[indx])  # angle in radians
        
        # Update new x,y,z based on new range
        pc_new[:,0] = ran_new*np.sin(the)*np.cos(phi)
        pc_new[:,1] = ran_new*np.sin(the)*np.sin(phi)
        pc_new[:,2] = ran_new*np.cos(the)
        pc_new[:,3] = ref_new
        
        return pc_new
    def msu_rain(self,pc,Rr):
        '''
        Lidar rain simulator from Goodin et al., 'Predicting the Influence of 
        Rain on LIDAR in ADAS', electronics 2019

        Parameters
        ----------
        pc : point cloud (N,4)
        Rr : rain rate in mm/hr

        Returns
        -------
        pc_new : output point cloud (N,4)

        '''
        shp    = pc.shape      # data shape
        pc_new = np.zeros(shp) # init new point cloud
        leng   = shp[0]        # data length
        
        # Rename variables for better readability
        x    = pc[:,0]
        y    = pc[:,1]
        z    = pc[:,2]
        ref  = pc[:,3]          
        
        # Get parameters from class init
        rmax = self.rmax       # max range (m)
        Pmin = 0.9*rmax**(-2)/np.pi  # min measurable power (arb units)
        
        # Calculate extinction coefficient from rain rate
        alpha = 0.01* Rr**0.6
        
        ran      = np.sqrt(x**2 + y**2 + z**2)  # range in m
        indv     = np.where(ran>0)[0] # clean data might already have invalid points
        P0       = np.zeros((leng,))
        P0[indv] = ref[indv]*np.exp(-2*alpha*ran[indv])/(ran[indv]**2) # calculate reflected power
        
        # init new ref and ran
        ran_new = np.zeros((leng,))
        ref_new = np.zeros((leng,))
        
        indp = np.where(P0>Pmin)[0] # points where power is greater than Pmin
        ref_new[indp] = ref[indp]*np.exp(-2*alpha*ran[indp]) # reflectivity reduced by atten
        sig = 0.02*ran[indp]* (1-np.exp(-Rr))**2
        ran_new[indp] = ran[indp] + np.random.normal(0,sig) # new range with uncertainty
        
        # Init angles
        phi = np.zeros((leng,))
        the = np.zeros((leng,))
        
        phi[indp] = np.arctan2(y[indp],x[indp])   # angle in radians
        the[indp] = np.arccos(z[indp]/ran[indp])  # angle in radians
        
        # Update new x,y,z based on new range
        pc_new[:,0] = ran_new*np.sin(the)*np.cos(phi)
        pc_new[:,1] = ran_new*np.sin(the)*np.sin(phi)
        pc_new[:,2] = ran_new*np.cos(the)
        pc_new[:,3] = ref_new
        
        return pc_new
    def haze_point_cloud(self,pts_3D,Rr=0):
        '''
        Modified from
        https://github.com/princeton-computational-imaging/SeeingThroughFog/blob/master/tools/DatasetFoggification/lidar_foggification.py

        Parameters
        ----------
        pts_3D : Point cloud
        Rr : Rain rate (mm/hr)

        Returns
        -------
        dist_pts_3d : Augmented point cloud
        '''
        n = []
        #Velodyne HDL64S2
        n = 0.05
        g = 0.35
        dmin = 2
            
        d = np.sqrt(pts_3D[:,0] * pts_3D[:,0] + pts_3D[:,1] * pts_3D[:,1] + pts_3D[:,2] * pts_3D[:,2])
        detectable_points = np.where(d>dmin)
        d = d[detectable_points]
        pts_3D = pts_3D[detectable_points]
        
        #######################################################################
        # This is the main modified part
        # For comparison we would like to calculate the extinction coefficient
        # from rain rate instead of sampling it from a distribution
        if (self.atm_model == 'rain') or (self.atm_model == 'snow'):
        	Nd  = self.N_model(self.D,Rr) # density of water droplets (m^-3)
        elif (self.atm_model == 'chu_hogg_fog') or (self.atm_model=='strong_advection_fog') or (self.atm_model=='moderate_advection_fog'):
        	Nd  = self.N_model(self.D) # density of water droplets (m^-3)
        else:
        	print('Warning: weather model not implemented')
        alpha, beta = self.alpha_beta(Nd)     # extinction coeff. (1/m)
        #######################################################################
    
        beta_usefull = alpha*np.ones(d.shape) # beta is the extinction coefficient (actually alpha)
        dmax = -np.divide(np.log(np.divide(n,(pts_3D[:,3] + g))),(2 * beta_usefull))
        dnew = -np.log(1 - 0.5) / (beta_usefull)
    
        probability_lost = 1 - np.exp(-beta_usefull*dmax)
        lost = np.random.uniform(0, 1, size=probability_lost.shape) < probability_lost
    
        cloud_scatter = np.logical_and(dnew < d, np.logical_not(lost))
        random_scatter = np.logical_and(np.logical_not(cloud_scatter), np.logical_not(lost))
        idx_stable = np.where(d<dmax)[0]
        old_points = np.zeros((len(idx_stable), 5))
        old_points[:,0:4] = pts_3D[idx_stable,:]
        old_points[:,3] = old_points[:,3]*np.exp(-beta_usefull[idx_stable]*d[idx_stable])
        old_points[:, 4] = np.zeros(np.shape(old_points[:,3]))
    
        cloud_scatter_idx = np.where(np.logical_and(dmax<d, cloud_scatter))[0]
        cloud_scatter = np.zeros((len(cloud_scatter_idx), 5))
        cloud_scatter[:,0:4] =  pts_3D[cloud_scatter_idx,:]
        cloud_scatter[:,0:3] = np.transpose(np.multiply(np.transpose(cloud_scatter[:,0:3]), np.transpose(np.divide(dnew[cloud_scatter_idx],d[cloud_scatter_idx]))))
        cloud_scatter[:,3] = cloud_scatter[:,3]*np.exp(-beta_usefull[cloud_scatter_idx]*dnew[cloud_scatter_idx])
        cloud_scatter[:, 4] = np.ones(np.shape(cloud_scatter[:, 3]))
    
    
        # Subsample random scatter abhaengig vom noise im Lidar
        random_scatter_idx = np.where(random_scatter)[0]
        scatter_max = np.min(np.vstack((dmax, d)).transpose(), axis=1)
        drand = np.random.uniform(high=scatter_max[random_scatter_idx])
        # scatter outside min detection range and do some subsampling. Not all points are randomly scattered.
        # Fraction of 0.05 is found empirically.
        drand_idx = np.where(drand>dmin)
        drand = drand[drand_idx]
        random_scatter_idx = random_scatter_idx[drand_idx]
        # Subsample random scattered points to 0.05%
        fraction_random = .05 # just set this according to the comment above^ rather than parsing arguments; also probably .05 not .05%
        subsampled_idx = np.random.choice(len(random_scatter_idx), int(fraction_random*len(random_scatter_idx)), replace=False)
        drand = drand[subsampled_idx]
        random_scatter_idx = random_scatter_idx[subsampled_idx]
    
    
        random_scatter = np.zeros((len(random_scatter_idx), 5))
        random_scatter[:,0:4] = pts_3D[random_scatter_idx,:]
        random_scatter[:,0:3] = np.transpose(np.multiply(np.transpose(random_scatter[:,0:3]), np.transpose(drand/d[random_scatter_idx])))
        random_scatter[:,3] = random_scatter[:,3]*np.exp(-beta_usefull[random_scatter_idx]*drand)
        random_scatter[:, 4] = 2*np.ones(np.shape(random_scatter[:, 3]))
    
        dist_pts_3d = np.concatenate((old_points, cloud_scatter,random_scatter), axis=0)
    
        return dist_pts_3d
    
    def calc_Mie_params(self):
        '''
        Calculate scattering efficiencies
        Returns
        -------
        D     : Particle diameter (mm)
        qext  : Extinction efficiency
        qback : Backscattering efficiency

        '''
        out   = ps.MieQ_withDiameterRange(self.m, self.lam, diameterRange=(1,1e7),
                                        nd=2000, logD=True)
        D     = out[0]*1e-6
        qext  = out[1]
        qback = out[6]
        
        # Save for later use since this function takes long to run
        np.savez('mie_q.npz',D=D,qext=qext,qback=qback)
        
        return D,qext,qback
    
    
    def alpha_beta(self,Nd):
        '''
        Calculates extunction and backscattering coefficients
        Parameters
        ----------
        Nd : particle size distribution, m^-3 mm^-1

        Returns
        -------
        alpha : extinction coefficient
        beta  : backscattering coefficient
        '''
        D  = self.D
        qe = self.qext
        qb = self.qback
        alpha = 1e-6*trapz(D**2*qe*Nd,D)*np.pi/4 # m^-1
        beta  = 1e-6*trapz(D**2*qb*Nd,D)*np.pi/4 # m^-1
        return alpha, beta
    
    # RAIN
    def N_MP_rain(self,D,Rr):
        '''
        Marshall - Palmer rain model

        Parameters
        ----------
        D  : rain droplet diameter (mm)
        Rr : rain rate (mm h^-1)

        Returns
        -------
        number of rain droplets for a given diameter (m^-3 mm^-1)
        '''
        return 8000*np.exp(-4.1*Rr**(-0.21)*D)
    
    def N_MP_tot_rain(self,Rr,dstart):
        '''
        Integrated Marshall - Palmer Rain model

        Parameters
        ----------
        Rr     : rain rate (mm h^-1)
        dstart : integral starting point for diameter (mm)

        Returns
        -------
        rain droplet density (m^-3) for a given min diameter
        '''
        lam = 4.1*Rr**(-0.21)
        return 8000*np.exp(-lam*dstart)/lam

    def MP_Sample_rain(self,Rr,N,dstart):
        '''
        Sample particle diameters from Marshall Palmer distribution

        Parameters
        ----------
        Rr     : rain rate (mm/hr)
        N      : number of samples
        dstart : Starting diameter (min diameter sampled)

        Returns
        -------
        diameters : diameter of the samples

        '''
        lmda      = 4.1*Rr**(-0.21)
        r         = np.random.rand(N)
        diameters = -np.log(1-r)/lmda + dstart
        return diameters
    
    # SNOW
    def N_MG_snow(self,D,Rr):
        '''
        Marshall - Palmer snow model

        Parameters
        ----------
        D  : snow diameter (mm)
        Rr : water equivalent rain rate (mm h^-1)

        Returns
        -------
        number of snow particles for a given diameter (m^-3 mm^-1)
        '''
        N0   = 7.6e3* Rr**(-0.87)
        lmda = 2.55* Rr**(-0.48)
        return N0*np.exp(-lmda*D)
    
    def N_MG_tot_snow(self,Rr,dstart):
        '''
        Integrated Marshall - Gunn snow model

        Parameters
        ----------
        Rr     : rain rate (mm h^-1)
        dstart : integral starting point for diameter (mm)

        Returns
        -------
        snow particle density (m^-3) for a given min diameter
        '''
        N0   = 7.6e3* Rr**(-0.87)
        lmda = 2.55* Rr**(-0.48)
        return N0*np.exp(-lmda*dstart)/lmda

    def MG_Sample_snow(self,Rr,N,dstart):
        '''
        Sample particle diameters from Marshall Palmer distribution

        Parameters
        ----------
        Rr     : rain rate (mm/hr)
        N      : number of samples
        dstart : Starting diameter (min diameter sampled)

        Returns
        -------
        diameters : diameter of the samples

        '''
        lmda      = 2.55* Rr**(-0.48)
        r         = np.random.rand(N)
        diameters = -np.log(1-r)/lmda + dstart
        return diameters
    # FOG
    def N_GD(self,D,rho,alpha,g,Rc):
        '''
        Gamma distribution model
        Note the parameters are NOT normalized to unitless values
        For example D^alpha term will have units Length^alpha
        It is therefore important to use exactly the same units for D as those
        cited in the paper by Rasshofer et al. and then perform unit conversion
        after an N(D) curve is generated
    
        D  : rain diameter
        Outputs number of rain droplets for a given diameter
        '''
        b = alpha/(g*Rc**g)
        
        Nd = g*rho*b**((alpha+1)/g)*(D/2)**alpha*np.exp(-b*(D/2)**g)/gamma((alpha+1)/g)
        
        return Nd
    # Coastal fog distribution
    # With given parameters, output has units cm^-3 um^-1 which is
    # then converted to m^-3 mm^-1 which is what alpha_beta() expects
    # so whole quantity is multiplied by (100 cm/m)^3 (1000 um/mm)
    def Nd_haze_coast(self,D):
        return 1e9*self.N_GD(D*1e3,rho=100,alpha=1,g=0.5,Rc=0.05e-3)
    
    # Continental fog distribution
    def Nd_haze_continental(self,D):
        return 1e9*self.N_GD(D*1e3,rho=100,alpha=2,g=0.5,Rc=0.07)
    
    # Strong advection fog
    def Nd_strong_advection_fog(self,D):
        return 1e9*self.N_GD(D*1e3,rho=20,alpha=3,g=1.,Rc=10)
    
    # Moderate advection fog
    def Nd_moderate_advection_fog(self,D):
        return 1e9*self.N_GD(D*1e3,rho=20,alpha=3,g=1.,Rc=8)
    
    # Strong spray
    def Nd_strong_spray(self,D):
        return 1e9*self.N_GD(D*1e3,rho=100,alpha=6,g=1.,Rc=4)
    
    # Moderate spray
    def Nd_moderate_spray(self,D):
        return 1e9*self.N_GD(D*1e3,rho=100,alpha=6,g=1.,Rc=2)
    
    # Chu/Hogg
    def Nd_chu_hogg(self,D):
        return 1e9*self.N_GD(D*1e3,rho=20,alpha=2,g=0.5,Rc=1)
    
