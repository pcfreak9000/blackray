#!/bin/python

import os
import sys
import numpy as np

# ----------------------------------------------------------------------------------------------------
# INPUT PARAMETERS
# ----------------------------------------------------------------------------------------------------

path_to_disk_data_file = "disk.csv"
path_to_xillver_file = "/home/roman/xillver-a-Ec5.fits"
gamma = 1.8
afe = 1
logxi = 1.
ecut = 300
incl = 50

spin = 0.0
a13 = 0
a22 = 0
epsi3 = 0
a52 = 0
alpha = -3

Nr = 200 #1000
Nph = 80 #400
Rmax = 80
rstep  = np.e**((np.log(Rmax)-np.log(np.cos(np.pi*incl/180)))/Nr);
pstep  = 2*np.pi/Nph;

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

rt_command = "./rt %f %f %f %f %f %f %f %f %f %s"%(spin, incl, a13, a22, a52, epsi3, alpha, rstep, pstep, path_to_disk_data_file)
# os.system("cd raytracer/")
print("Start of ray-tracing part ...")
os.system(rt_command)
print("End of ray-tracing part")

xillver_command = "python get_xillver_spectrum.py %f %f %f %f %f %s"%(gamma, afe, logxi, ecut, incl, path_to_xillver_file)
os.system(xillver_command)

# photons_data_filename = 'data/photons_data_a9.40000e-01.i7.00e+01.e_0.00e+00.a13_0.00e+00.a22_0.00e+00.a52_0.00e+00.dat'

photons_data_filename = "data/photons_data_a%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%.05Lf.dat"%(spin,incl,epsi3,a13,a22,a52)
final_filename = "output/spectrum_a_%.5f_i_%.3f_a13_%.5f_a22_%.5f_gam_%.2f_afe_%.2f_xi_%.2f_ecut_%.2f.dat"%(spin, incl, a13, a22, gamma, afe, logxi, ecut)
xill_spec_filename = "data/xill_spec_gam_%.2f_afe_%.2f_xi_%.2f_ecut_%.2f_incl_%.2f.dat" % (gamma, afe, logxi, ecut, incl)
convolver_command = "./conv %s %s %s %f %f %f %f" % (photons_data_filename, xill_spec_filename, final_filename, incl, logxi, alpha, rstep)
print("Start of convolution ...")
os.system(convolver_command)
print("End of convolution! The resulting spectrum is in output/")
