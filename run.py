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
#    os.system("mkdir "+tempdir+"data")
#    os.system("mkdir "+tempdir+"output")
#    os.system("mkdir "+tempdir+"ironline_data")
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

    Nr = 200 #1000
    Nph = 80 #400
    Rmax = 200
    rstep  = np.e**((np.log(Rmax)-np.log(np.abs(np.cos(np.pi*incl/180))))/Nr);
    pstep  = 2*np.pi/Nph;
    print("rstep: "+str(rstep))
    print("pstep: "+str(pstep))
    # ----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------
    outtxt = tempdir+"rays_debug.dat"
    outgtxt = tempdir+"gfactor.png"
    outitxt = tempdir+"hitindex.png"
    illumfilename = "data_lp_pr.dat"
    rt_command = "ulimit -c unlimited && ./rt Metric,spin=%f,Raytracing,incl=%f,Metric,a13=%f,a22=%f,a52=%f,epsi3=%f,Corona,alpha=%f,Raytracing,rstep=%f,pstep=%f %s %s %s %s"%(spin, incl, a13, a22, a52, epsi3, alpha, rstep, pstep, path_to_disk_data_file,tempdir,outtxt,illumfilename)
    print("Start of ray-tracing part ...")
    dddd = os.system(rt_command)
    if dddd != 0:
        print(dddd)
    print("End of ray-tracing part")

    print("Starting debug image processes...")
    img_command = "python img.py %s &"%(tempdir)
    os.system(img_command)

    print("Fetching xillver spectrum")
    xillver_command = "python get_xillver_spectrum.py %f %f %f %f %f %s %s"%(gamma, afe, logxi, ecut, incl, path_to_xillver_file,tempdir)
    os.system(xillver_command)

    photons_data_filename= tempdir+"photons.dat"
    final_filename = tempdir+"spectrum.dat"

    xill_spec_filename = tempdir+"xillver_spectrum.dat"
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
