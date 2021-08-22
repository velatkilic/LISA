# Lidar Scattering Augmentation (LISA)
LISA treats effects of small scatterers on average using the Beer-Lambert law and physical models of a bistatic lidar. Rain and snow models use a hybrid Monte-Carlo method to augment clean lidar point clouds for a given rain rate. LISA can be easily extended to other scatterers with different droplet size models or materials (i.e dust).

Requirements: Numpy, Scipy, PyMieScatt

Install PyMieScatt: _pip install PyMieScatt_

## Example usage
Copy _atmos_models_ into your project

Initialize lisa object and pass relevant parameters (see class description for details):
> lisa = LISA(atm_model='rain')

Call augment method to add noise to the clean point cloud data
> data_r33 = lisa.augment(data_clear,33) # rain rate 33 mm/hr

## Reference
This version was used in the following paper:
[V. Kilic, D. Hegde, V. Sindagi, A.B. Cooper, M.A. Foster and V.M. Patel,
"Lidar Light Scattering Augmentation (LISA): Physics-based Simulation of Adverse Weather Conditions for 3D Object Detection",
arXiv:2107.07004, (2021).](https://arxiv.org/abs/2107.07004)