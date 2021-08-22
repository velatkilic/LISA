# Lidar Light Scattering Augmentation (LISA)
LISA is a physics based augmentation method that models effects of adverse weather conditions on lidar data. 
It treats effects of small scatterers on average using the Beer-Lambert law where extinction coefficients are
calculated thorugh Mie theory and particle size distributions. Rain and snow models use a hybrid Monte-Carlo 
method to augment clean lidar point clouds for a given rain rate.  LISA can be easily extended to other
scatterers with different droplet size models or materials (i.e dust).

One of these effects is reduced range due to attenuation of the laser signal due to scatterers. Note sparsity
of the rainy scenes at large distances from the sensor:
![Reduced range](/images/rain.png)

Particles near the sensor can backscatter enough light and cause false detections (blue box). Also note "fuzziness"
of the scan lines under adverse weather arising from reduced SNR (signal to noise ratio) which leads to range uncertainty.
![Randomly scattered points near sensor and lost points](/images/fog_snow.png)

# Using LISA
Option 1: Build from the source and install using pip:

    pip install .

Option 2: (todo) This option is not working yet

    pip install lisa
    
Basic usage
Create a Lisa object with default options

    augm = lisa.Lisa()

or create a Lisa object with non default options

    lidar = lisa.Lidar() # lidar object
    water = lisa.Water() # material object
    rain  = lisa.MarshallPalmerRain() # particle distribution model
    
    augm  = lisa.Lisa(lidar, water, rain)
    
Next feed point cloud data (2D array) to the augmentor
    
    # Point cloud format: x,y,z,reflectance
    pc    = np.array([[1.0, 0., 0., 0.5]]) 
    
    # pcnew should be [0,0,0,0,0] since range = sqrt(1+0+0) = 1 is smaller than min lidar range
    pcnew = augm.augment(pc) 

## Reference
Cite as 

[V. Kilic, D. Hegde, V. Sindagi, A.B. Cooper, M.A. Foster and V.M. Patel,
"Lidar Light Scattering Augmentation (LISA): Physics-based Simulation of Adverse Weather Conditions for 3D Object Detection",
arXiv:2107.07004, (2021).](https://arxiv.org/abs/2107.07004)
