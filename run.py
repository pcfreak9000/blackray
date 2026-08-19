#! /usr/bin/env python

import os
import sys
import numpy as np
import math
import argparse

def main(**kwargs):
    # ----------------------------------------------------------------------------------------------------
    # INPUT PARAMETERS
    # ----------------------------------------------------------------------------------------------------
    tempdir = kwargs['temp']+"/"
    os.system("mkdir "+tempdir)
    os.system("mkdir "+tempdir+"data")
    os.system("mkdir "+tempdir+"output")
    os.system("mkdir "+tempdir+"ironline_data")
    path_to_disk_data_file = kwargs['disk']
    path_to_xillver_file = kwargs['xillver'] #"/home/roman/xillver-a-Ec5.fits"
    gamma = 1.8
    afe = 1
    logxi = 3.1
    ecut = 300
    incl = kwargs['incl']

    spin = kwargs['spin']
    a13 = 0
    a22 = 0
    epsi3 = 0
    a52 = 0
    alpha = kwargs['alpha']

    Nr = 1000 #1000
    Nph = 400 #400
    Rmax = 200
    rstep  = np.e**((np.log(Rmax)-np.log(np.abs(np.cos(np.pi*incl/180))))/Nr);
    pstep  = 2*np.pi/Nph;
    print("rstep: "+str(rstep))
    print("pstep: "+str(pstep))
    # ----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------
    outtxt = tempdir+"output/a_%.5f_i_%.3f_a13_%.5f_a22_%.5f_gam_%.2f_afe_%.2f_xi_%.2f_ecut_%.2f_debugdata.txt"%(spin, incl, a13, a22, gamma, afe, logxi, ecut)
    outgtxt = tempdir+"output/a_%.5f_i_%.3f_a13_%.5f_a22_%.5f_gam_%.2f_afe_%.2f_xi_%.2f_ecut_%.2f_g.png"%(spin, incl, a13, a22, gamma, afe, logxi, ecut)
    outitxt = tempdir+"output/a_%.5f_i_%.3f_a13_%.5f_a22_%.5f_gam_%.2f_afe_%.2f_xi_%.2f_ecut_%.2f_i.png"%(spin, incl, a13, a22, gamma, afe, logxi, ecut)
    rt_command = "ulimit -c unlimited && ./rt %f %f %f %f %f %f %f %f %f %s %s %s"%(spin, incl, a13, a22, a52, epsi3, alpha, rstep, pstep, path_to_disk_data_file,tempdir,outtxt)
    # os.system("cd raytracer/")
    print("Start of ray-tracing part ...")
    dddd = os.system(rt_command)
    if dddd != 0:
        print(dddd)
    print("End of ray-tracing part")
    
    print("Starting debug image processes...")
    img_command = "python img.py %s %s %s &"%(outtxt,outgtxt,outitxt)
    os.system(img_command)
    
    print("Fetching xillver spectrum")
    xillver_command = "python get_xillver_spectrum.py %f %f %f %f %f %s %s"%(gamma, afe, logxi, ecut, incl, path_to_xillver_file,tempdir)
    os.system(xillver_command)

    # photons_data_filename = 'data/photons_data_a9.40000e-01.i7.00e+01.e_0.00e+00.a13_0.00e+00.a22_0.00e+00.a52_0.00e+00.dat'

    photons_data_filename = tempdir+"data/photons_data_a%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%.05Lf.dat"%(spin,incl,epsi3,a13,a22,a52)
    final_filename = tempdir+"output/a_%.5f_i_%.3f_a13_%.5f_a22_%.5f_gam_%.2f_afe_%.2f_xi_%.2f_ecut_%.2f_spectrum.dat"%(spin, incl, a13, a22, gamma, afe, logxi, ecut)
    xill_spec_filename = tempdir+"data/xill_spec_gam_%.2f_afe_%.2f_xi_%.2f_ecut_%.2f_incl_%.2f.dat" % (gamma, afe, logxi, ecut, incl)
    convolver_command = "./conv %s %s %s %f %f %f %f" % (photons_data_filename, xill_spec_filename, final_filename, incl, logxi, alpha, rstep)
    print("Start of convolution ...")
    aaaa = os.system(convolver_command)
    if aaaa != 0:
        print(aaaa)
    print("End of convolution! The resulting spectrum is in output/")
    

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('disk', type=str, help="Path to disk file")
    parser.add_argument('spin', type=float, help="BH spin")
    parser.add_argument('incl', type=float, help="Inclination in degrees")
    parser.add_argument('alpha', type=float, help="alpha parameter")
    parser.add_argument('xillver', type=str, help="Path to xillver table file")
    parser.add_argument('temp', type=str, help="Path to temp dir")
    #parser.add_argument('input_filename',
    #                    type=str,
    #                    help='base name of files to be converted, including directory')
    #parser.add_argument('output_filename',
    #                    type=str,
    #                    help='name of new files to be saved, including directory')
    #parser.add_argument('start',
    #                    type=int,
    #                    help='first file number to be converted')
    #parser.add_argument('end',
    #                    type=int,
    #                    help='last file number to be converted')
    #parser.add_argument('stride',
    #                    type=int,
    #                    default=0,
    #                    help='stride in file numbers to be converted')
    #parser.add_argument('--max', type=float, default=None, help='limit used values to max')
    #parser.add_argument('-q', '--quantities',
    #        type=str,
    #        nargs='+',
    #        help='quantities to be averaged')
    #parser.add_argument('--lvl', type=int, default=0, help='Min level to avg over')
    args = parser.parse_args()
    main(**vars(args))
