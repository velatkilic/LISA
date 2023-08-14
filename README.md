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

For original paper code can be found [here](/python_old/). A faster version is currently being developed.

## Reference
Cite as
```
@article{kilic2021lidar,
  title={Lidar light scattering augmentation (lisa): Physics-based simulation of adverse weather conditions for 3d object detection},
  author={Kilic, Velat and Hegde, Deepti and Sindagi, Vishwanath and Cooper, A Brinton and Foster, Mark A and Patel, Vishal M},
  journal={arXiv preprint arXiv:2107.07004},
  year={2021}
}

@inproceedings{hegde2023source,
  title={Source-free Unsupervised Domain Adaptation for 3D Object Detection in Adverse Weather},
  author={Hegde, Deepti and Kilic, Velat and Sindagi, Vishwanath and Cooper, A Brinton and Foster, Mark and Patel, Vishal M},
  booktitle={2023 IEEE International Conference on Robotics and Automation (ICRA)},
  pages={6973--6980},
  year={2023},
  organization={IEEE}
}
```
