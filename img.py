#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def normalize_column(col):
    return (col - np.min(col)) / (np.max(col) - np.min(col))

def main(input_file, output_image, size=256, use_labels=False):
    # Load the data
    data = np.loadtxt(input_file, delimiter=None)

    if data.shape[1] < 4:
        raise ValueError("Input file must have at least four columns: x y data label")

    x_raw, y_raw, d_raw, label_raw = data[:, 0], data[:, 1], data[:, 2], data[:, 3].astype(int)

    # Normalize x and y
    x_norm = normalize_column(x_raw)
    y_norm = normalize_column(y_raw)

    # Create grid
    grid_x, grid_y = np.meshgrid(
        np.linspace(0, 1, size),
        np.linspace(0, 1, size)
    )

    if use_labels:
        # Define a color map as a dictionary of RGB tuples (values in range 0–1)
        label_colors = {
            255: (1.0, 0.0, 0.0),   # too many iterations, red
            6: (1.0, 0.0, 1.0),   # numerical problems, magenta
            3: (0.0, 1.0, 0.0),   # 
            2: (0.0, 0.0, 1.0),   # miss, blue
            1: (1.0, 1.0, 0.0),   # hit, yellow
            4: (0.0, 0.0, 0.0),   # horizon cross, black
            5: (0.0, 0.0, 0.0),   # horizon cross, black
            7: (0.3, 0.3, 0.3)    # escape to infinity, gray
        }

        # Prepare RGB channels
        grid_r = griddata((x_norm, 1-y_norm), [label_colors[l][0] for l in label_raw], (grid_x, grid_y), method='nearest', fill_value=0)
        grid_g = griddata((x_norm, 1-y_norm), [label_colors[l][1] for l in label_raw], (grid_x, grid_y), method='nearest', fill_value=0)
        grid_b = griddata((x_norm, 1-y_norm), [label_colors[l][2] for l in label_raw], (grid_x, grid_y), method='nearest', fill_value=0)

        # Stack to RGB image
        rgb_image = np.stack([grid_r, grid_g, grid_b], axis=2)
        plt.imsave(output_image, rgb_image)
    
    else:
        # Grayscale mode: 1 - normalized(data)
        d_inverted = 1 - normalize_column(d_raw)
        grid_d = griddata((x_norm, 1-y_norm), d_inverted, (grid_x, grid_y), method='cubic', fill_value=0)
        plt.imsave(output_image, grid_d, cmap='gray', vmin=0, vmax=1)

if __name__ == "__main__":
    # Set `use_labels` to True to enable color labeling
    main("output.txt", "output.png", size=2048, use_labels=True)

