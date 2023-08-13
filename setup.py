from setuptools import setup, find_packages

with open("README.md", "r") as file:
    long_description = file.read()

setup(
    name='pylisa',
    version='0.0.2',
    author='Velat Kilic',
    description='Lidar light scattering augmentation',
    long_description=long_description,
    packages=find_packages(),
    keywords=['lidar', 'augmentation', 'scattering', 'lidar augmentation'],
    classifiers = [
        "Development Status :: 2 - Pre-Alpha",
    ]
)