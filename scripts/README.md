# Read, Process and Write Augmented Data
`run.py` is a script to read a cvbin file with corresponding `.yaml` file and output a cvbin file with augmented point cloud using simulator LISA.

## Using the Script
An example config file is located at `LISA/scripts/lisa_example.yaml`, you need to configure it by yourself and all the information is should be inside it. Then by running:
>  python run.py --config scripts/lisa_example.yaml
you can get the output cvbin file.