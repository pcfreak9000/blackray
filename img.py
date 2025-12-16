#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.colors as mcolors
from sys import argv
from multiprocessing import Pool

def normalize_column(col):
    return (col - np.min(col)) / (np.max(col) - np.min(col))

def gengrid(args):
    return griddata(args[0], args[1], args[2], args[3], args[4])

def main(preptuple):
    data = preptuple[0]
    x_raw = preptuple[1]
    y_raw = preptuple[2]
    d_raw = preptuple[3]
    label_raw = preptuple[4]
    x_norm = preptuple[5]
    y_norm = preptuple[6]
    grid_x = preptuple[7]
    grid_y = preptuple[8]
    output_image = preptuple[9]
    use_labels = preptuple[10]

    if use_labels:
        # Define a color map as a dictionary of RGB tuples (values in range 0–1)
        label_colors = {
            255: (1.0, 0.0, 0.0),   # too many iterations, red
            6: (1.0, 0.0, 1.0),   # numerical problems, magenta
            130: (1.0, 0.647, 0.0),   # backside hit, orange
            131: (0.9, 0.583, 0.0),  # bakside hit, different
            2: (0.0, 0.0, 1.0),   # miss, blue
            128: (1.0, 1.0, 0.0),   # hit, yellow
            129: (0.9, 0.9, 0.0),  # hit, different element
            4: (0.0, 0.0, 0.0),   # horizon cross, black
            5: (0.0, 0.2, 0.0),   # horizon cross, green-black
            7: (0.3, 0.3, 0.3)    # escape to infinity, gray
        }
        interpolmeth='linear'
        # Prepare RGB channels
        #grid_r = gengrid(((x_norm, 1.0-y_norm), [label_colors[l][0] for l in label_raw], (grid_x, grid_y), interpolmeth, 0.0))
        #grid_g = gengrid(((x_norm, 1.0-y_norm), [label_colors[l][1] for l in label_raw], (grid_x, grid_y), interpolmeth, 0.0))
        #grid_b = gengrid(((x_norm, 1.0-y_norm), [label_colors[l][2] for l in label_raw], (grid_x, grid_y), interpolmeth, 0.0))
        lbs = []
        for i in range(3):
            lbs.append(((x_norm, 1.0-y_norm), [label_colors[l][i] for l in label_raw], (grid_x, grid_y), interpolmeth, 0.0))
        with Pool(3) as p:
            grid_r, grid_g, grid_b = p.map(gengrid, lbs)
        
        # Stack to RGB image
        rgb_image = np.stack([grid_r, grid_g, grid_b], axis=2)
        rgb_image_cl = np.clip(rgb_image, min=0.0, max=1.0)
        plt.imsave(output_image, rgb_image_cl)

    else:
        grid_d = gengrid(((x_norm, y_norm), np.log10(d_raw), (grid_x, grid_y), 'nearest', 0.0))

        #absmax = np.nanmax(np.abs(grid_d))
        #norm = mcolors.TwoSlopeNorm(vmin=-absmax, vcenter=0.0, vmax=absmax)
        #should yield the same result as the above twoslopenorm:
        norm = mcolors.CenteredNorm()

        #lmax = np.nanmax(grid_d)
        #lmin = np.nanmin(grid_d)
        #logmax = np.abs(np.log10(lmax))
        #logmin = np.abs(np.log10(lmin))
        #if logmax > logmin:
        #    lmin = 1.0/lmax;
        #else:
        #    lmax = 1.0/lmin;
        #norm = mcolors.LogNorm(vmin=lmin, vmax=lmax)

        plt.figure(figsize=(6, 5))
        img = plt.imshow(grid_d, cmap='seismic_r', norm=norm, origin='lower', extent=[0, 1, 0, 1])
        cbar = plt.colorbar(img)
        cbar.set_label(r'$\log_{10}(g)$')
        #plt.title("Interpolated Data (log scale)")
        plt.xlabel("Normalized X")
        plt.ylabel("Normalized Y")
        plt.tight_layout()
        plt.savefig(output_image, dpi=600)
        plt.close()

if __name__ == "__main__":
    outtxt = argv[1]
    outg = argv[2]
    outi = argv[3]
    
    size=2048
    
    # Set `use_labels` to True to enable color labeling
    # Load the data
    data = np.loadtxt(outtxt, delimiter=None)

    if data.shape[1] < 4:
        raise ValueError("Input file must have at least four columns: x y data label")

    x_raw, y_raw, d_raw, label_raw = data[:, 0], data[:, 1], data[:, 2], data[:, 3].astype(int)
    
    if len(x_raw)<4:
        print("img.py: not enough data for interpolation, aborting")
        exit()

    # Normalize x and y
    x_norm = normalize_column(x_raw)
    y_norm = normalize_column(y_raw)

    # Create grid
    grid_x, grid_y = np.meshgrid(
        np.linspace(0, 1, size),
        np.linspace(0, 1, size)
    )
    
    main((data, x_raw, y_raw, d_raw, label_raw, x_norm, y_norm, grid_x, grid_y, outg, False))
    main((data, x_raw, y_raw, d_raw, label_raw, x_norm, y_norm, grid_x, grid_y, outi, True))

