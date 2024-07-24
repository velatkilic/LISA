import rerun as rr
import json
import os
import struct
import yaml
import numpy as np
import argparse
import cv_pb2
import atmos_models
from scipy.spatial.transform import Rotation as R


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', help='path to script config file')
    # parser.add_argument('--context-cvbin', help='path to context cvbin file with sensor calibration parameters')
    args = parser.parse_args()

    return args

def parse_config(config_file):
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    
    return config

def parse_cvbin_paths(orignal_cvbin_path, out_dir):
    assert orignal_cvbin_path.endswith('.cvbin'), "This file does not ends with .cvbin"
    original_cvbin_name = orignal_cvbin_path.split('/')[-1]
    rainy_cvbin_path = os.path.join(out_dir, original_cvbin_name[:-6] + "_augmented" + original_cvbin_name[-6:])

    return orignal_cvbin_path, rainy_cvbin_path

def parse_transforms(files_path):
    transform_matrix = np.identity(4)
    for file_path, choice in files_path:
        assert choice in ['forward', 'inverse'], "transform choice can only be \' forward \' or \'inverse\'"

        matrix = parse_transform(file_path, choice)
        transform_matrix = matrix.dot(transform_matrix)
    
    return transform_matrix

def parse_transform(file_path, choice):
    matrix = np.eye(4)

    try:
        with open(file_path, "r") as file:
            data = yaml.safe_load(file)

        transform = data.get('transform', {})

        rotation = transform.get('rotation', {})
        rotation = R.from_quat([rotation['x'], rotation['y'], rotation['z'], rotation['w']]).as_matrix()

        translation = transform.get('translation', {})
        translation = [translation['x'], translation['y'], translation['z']]

    except Exception:
        print(f"Error occurs. Try to ignore loading {file}")

        rotation = np.identity(3)
        translation = np.zeros(3)

    matrix[:3, :3] = rotation
    matrix[:3, 3] = translation

    if choice == "forward":
        matrix = np.linalg.inv(matrix)

    return matrix

def deserialize_cvbin_content(cvbin_path):
    """Deserialize the content of a cvbin file."""
    cvbin_content = cv_pb2.CVBinContent()
    cvbin_bytes = open(cvbin_path, 'rb').read()
    cvbin_content.ParseFromString(cvbin_bytes)

    return cvbin_content

def serialize_cvbin_content(cvbin_content, output_path): 
    """Serialize a cvbin content to a file.""" 
    with open(output_path, 'wb') as f: 
        f.write(cvbin_content.SerializeToString())

def decode_lidar_data(lidar_data, values_per_point = 4):
    """
    Decode the raw bytes of the LIDAR data as a list of points.
    A single point in the LIDAR data consists of the 4 or 5 fields: x, y, z, and intensity with optional label. Each field is a 4-byte float.
    """
    bytes_per_value = 4

    num_values = len(lidar_data) // bytes_per_value
    points = struct.unpack(f'{num_values}f', lidar_data)
    # Reshape to a 2D array whre each row is a point.
    points = np.array(points).reshape(-1, values_per_point)

    return points

def transform_lidar(pcs, transform_matrix):
    for idx in range(len(pcs)):
        coordinate = np.append(pcs[idx][:3], 1)
        pcs[idx][:3] = transform_matrix.dot(coordinate)[:3]

    return pcs

def simulate_lisa(pcs, config_lisa):
    lisa = atmos_models.LISA(m=config_lisa['refractive_index_contrast'],
                                  lam=config_lisa['wavelength'],
                                  rmax=config_lisa['max_lidar_range'],
                                  rmin=config_lisa['min_lidar_range'],
                                  bdiv=config_lisa['beam_divergence_angle'],
                                  dst=config_lisa['droplet_diameter_starting_point'],
                                  dR=config_lisa['range_accuracy'],
                                  saved_model=config_lisa['use_saved_mie_coefficients'],
                                  atm_model=config_lisa['atmospheric_model_type'],
                                  mode=config_lisa['lidar_return_mode'])
    
    pcs_augmented = lisa.augment_mc(pcs, config_lisa['rain_rate'])
    return pcs_augmented

def encode_lidar_data(lidar_data):
    return lidar_data.astype(np.float32).tobytes()

def main():
    args = parse_arguments()

    config = parse_config(args.config)

    orignal_cvbin_path, rainy_cvbin_path = parse_cvbin_paths(config['cvbin_path'], config['out_dir'])
    cvbin_data = deserialize_cvbin_content(orignal_cvbin_path)

    transform_matrix = parse_transforms(config['operations'])

    lidars_rain = []
    for lidar in cvbin_data.multiple_sensor_scene.lidars:
        '''for each channel'''

        decoded_lidar = decode_lidar_data(lidar.raw)
        decoded_lidar = transform_lidar(decoded_lidar, transform_matrix)

        lidar_augmented = simulate_lisa(decoded_lidar, config['lisa'])

        if config['visualize']:
            # points_lost = np.array([point for point in decoded_lidar_ if point[-1] == 0])
            points_scattered = np.array([point for point in lidar_augmented if point[-1] == 1])
            points_unchanged = np.array([point for point in lidar_augmented if point[-1] == 2])

            rr.init("point_cloud", spawn=True)
            # rr.log("/points_lost", rr.Points3D(points_lost[:,:3], colors=[255, 0, 0]))
            rr.log("/points_original", rr.Points3D(decoded_lidar[:,:3], colors=[255, 0, 0]))
            rr.log("/points_scattered", rr.Points3D(points_scattered[:,:3], colors=[0, 255, 0]))
            rr.log("/points_unchanged", rr.Points3D(points_unchanged[:,:3], colors=[0, 0, 255]))
        
        lidar_augmented = transform_lidar(lidar_augmented, np.linalg.inv(transform_matrix))
        lidar.raw = encode_lidar_data(lidar_augmented)

        lidar.field.append("label")
        lidars_rain.append(lidar)

    del cvbin_data.multiple_sensor_scene.lidars[:]
    cvbin_data.multiple_sensor_scene.lidars.extend(lidars_rain)
    
    serialize_cvbin_content(cvbin_data, rainy_cvbin_path)

if __name__ == "__main__":
    main()