import numpy as np
from scipy.special import gamma
from scipy.integrate import trapz

def alpha_beta(nd, d, q_ext, q_back):
    '''
    Calculates extunction and backscattering coefficients
    Parameters
    ----------
    nd : particle size distribution, m^-3 mm^-1
    d  : particle diameter, mm
    q_ext: extinction efficiency
    q_back: back scattering efficiency

    Returns
    -------
    alpha : extinction coefficient
    beta  : backscattering coefficient
    '''
    alpha = 1e-6*trapz(d**2*q_ext*nd,d)*np.pi/4 # m^-1
    beta  = 1e-6*trapz(d**2*q_back*nd,d)*np.pi/4 # m^-1
    return alpha, beta

# RAIN
def n_mp_rain(d, Rr):
    '''
    Marshall - Palmer rain model

    Parameters
    ----------
    d  : rain droplet diameter (mm)
    Rr : rain rate (mm h^-1)

    Returns
    -------
    number of rain droplets for a given diameter (m^-3 mm^-1)
    '''
    return 8000*np.exp(-4.1*Rr**(-0.21)*d)

def n_mp_tot_rain(Rr, dstart):
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

def mp_sample_rain(Rr,N,dstart):
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
def n_mg_snow(d, Rr):
    '''
    Marshall - Gunn snow model

    Parameters
    ----------
    d  : snow diameter (mm)
    Rr : water equivalent rain rate (mm h^-1)

    Returns
    -------
    number of snow particles for a given diameter (m^-3 mm^-1)
    '''
    N0   = 7.6e3* Rr**(-0.87)
    lmda = 2.55* Rr**(-0.48)
    return N0*np.exp(-lmda*d)

def n_mg_tot_snow(Rr, dstart):
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

def mg_sample_snow(Rr, N, dstart):
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
def n_gd(d, rho, alpha, g, Rc):
    '''
    Gamma distribution model
    Note the parameters are NOT normalized to unitless values
    For example d^alpha term will have units Length^alpha
    It is therefore important to use exactly the same units for d as those
    cited in the paper by Rasshofer et al. and then perform unit conversion
    after an N(d) curve is generated

    d  : rain diameter
    Outputs number of rain droplets for a given diameter
    '''
    b = alpha/(g*Rc**g)
    
    Nd = g*rho*b**((alpha+1)/g)*(d/2)**alpha*np.exp(-b*(d/2)**g)/gamma((alpha+1)/g)
    
    return Nd

# Coastal fog distribution
# With given parameters, output has units cm^-3 um^-1 which is
# then converted to m^-3 mm^-1 which is what alpha_beta() expects
# so whole quantity is multiplied by (100 cm/m)^3 (1000 um/mm)
def nd_haze_coast(d):
    return 1e9*n_gd(d*1e3,rho=100,alpha=1,g=0.5,Rc=0.05e-3)

# Continental fog distribution
def nd_haze_continental(d):
    return 1e9*n_gd(d*1e3,rho=100,alpha=2,g=0.5,Rc=0.07)

# Strong advection fog
def nd_strong_advection_fog(d):
    return 1e9*n_gd(d*1e3,rho=20,alpha=3,g=1.,Rc=10)

# Moderate advection fog
def nd_moderate_advection_fog(d):
    return 1e9*n_gd(d*1e3,rho=20,alpha=3,g=1.,Rc=8)

# Strong spray
def nd_strong_spray(d):
    return 1e9*n_gd(d*1e3,rho=100,alpha=6,g=1.,Rc=4)

# Moderate spray
def nd_moderate_spray(d):
    return 1e9*n_gd(d*1e3,rho=100,alpha=6,g=1.,Rc=2)

# Chu/Hogg
def nd_chu_hogg(d):
    return 1e9*n_gd(d*1e3,rho=20,alpha=2,g=0.5,Rc=1)