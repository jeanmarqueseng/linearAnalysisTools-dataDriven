# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#
# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#
# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#
# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#

# This code processes the snapshot matrices and computes:
# ----- POD singular values and vectors
# ----- DMD eigenvalues and vectors
# ----- not the spatial modes, these are done afterwards on outputModes.py

# ----- Department of Mechanical and Aerospace Engineering
# ----- University of California, Los Angeles
# ----- Author: Jean Helder Marques Ribeiro
# ----- email:  jeanmarques@g.ucla.edu
# ----- Date: August/2021
import os,sys,getopt
sys.path.append('./labTools')

import numpy as np
import scipy.linalg as sp
from mpi4py import MPI

from dataTools import *

from scipy.io import loadmat,savemat

import argparse
import time

if __name__ == '__main__':

    # communications for MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    start_time = time.time()

    # ---------------------- INPUT ARGS ---------------------#
    # ---------------------- INPUT ARGS ---------------------#
    # ---------------------- INPUT ARGS ---------------------#
    # ---------------------- INPUT ARGS ---------------------#
    parser = argparse.ArgumentParser()
    parser.add_argument('-domain',type=int, nargs=1,help='<whole> only')
    parser.add_argument('-outputdir',type=str,help='Output folder for the binary files')
    parser.add_argument('-spod',type=int,help='Sieber-POD filter size: 0 (no filter), # of snapshots (DFT)')
    args = parser.parse_args()

    if(args.domain):
        domain = args.domain
    else:
        domain = False

    if(args.outputdir):
        outdir = args.outputdir
    else:
        outdir = './flowfield/'

    if(args.spod):
        spod = args.spod
    else:
        spod = 0

    # -------------- POD and DMD modes -----------------#   
    # -------------- POD and DMD modes -----------------#
    # -------------- POD and DMD modes -----------------# 
    # -------------- POD and DMD modes -----------------#
    if (rank == 0):
        # number of variables
        qVar = np.arange(0,3,1)
        nVar = qVar.shape[0]                                        
        # check number of snapshots with the first correlation matrix for RHO
        if (domain):
            matLoad    = loadmat(outdir + '/%s.%02d.autoCorrelation.mat' % (var(qVar[0]),0))
        else:
            matLoad    = loadmat(outdir + '/%s.autoCorrelation.mat' % (var(qVar[0])))
        T = matLoad['C'].shape[0]
        if (spod > T): spod = T-1
        
        matLoad    = loadmat(outdir + '/%s.spectralCorrelation.mat' % (var(qVar[0])))
        F = matLoad['CSD'].shape[0]
        N = matLoad['CSD'].shape[1]
        
        # Read correlation matrix for each variable
        # if domain splitting in autoCorrelation.py, remember to sum all C's
        C     = np.zeros((np.amax(qVar)+1,T,T),dtype=np.float32)
        for iv in qVar:
            if (domain):
                for i in range(domain[0]):
                    print('... reading %s' % (outdir + '/%s.%02d.autoCorrelation.mat' % (var(qVar[iv]),i)))
                    matLoad    = loadmat(outdir + '/%s.%02d.autoCorrelation.mat' % (var(qVar[iv]),i))
                    C[iv,:,:] += matLoad['C']
            else:
                matLoad    = loadmat(outdir + '/%s.autoCorrelation.mat' % (var(qVar[iv])))
                C[iv,:,:] = matLoad['C']

        # --------------- STANDARD POD and DMD ----------------------#
        # --------------- STANDARD POD and DMD ----------------------#
        # --------------- STANDARD POD and DMD ----------------------#
        # --------------- STANDARD POD and DMD ----------------------#
        # Compute POD/DMD norm
        Cnorm = C[0,:,:] + C[1,:,:] + C[2,:,:]              # kinetic energy norm
        #Cnorm = C[3,:,:] + C[4,:,:] + C[5,:,:]              # enstrophy norm
        #Cnorm = C[6,:,:]                                    # pressure norm

        # apply SPOD filter (SPOD by Sieber et al. (JFM, 2016))
        Cnorm    = spod_filter(Cnorm,spod)
        cc       = Cnorm[0:-1,0:-1]                         # Correlation matrix related to X
        cp       = Cnorm[0:-1,1:]                           # Correlation matrix related to X'
        
        # compute SVD of normalized correlation matrix
        uc,sc,vc = sp.svd(cc)
        sc       = np.sqrt(sc)

        savemat(outdir + '/singularValues.mat',    {'S':sc})
        savemat(outdir + '/singularVectors.mat',   {'U':uc})
        savemat(outdir + '/correlationMatrix.mat', {'C':cc})
        
        # I did not truncate the POD modes, but you can do that later if you want just by doing this:
        # trunc = 5
        # sc = sc[0:trunc]
        # uc = uc[0:trunc,:]
        # vc = vc[0:trunc,:]

        # compute DMD modes (check J. Tu Thesis if you want, but really, it's just this)
        # Transition matrix A
        A   =  np.diag(1/sc)@vc@cp@vc.transpose().conj()@np.diag(1/sc)
        
        mu,v = sp.eig(A) 

        savemat(outdir + '/eigenValues.mat',  {'mu':mu})
        savemat(outdir + '/eigenVectors.mat', {'V':v})
        
        # --------------- SPECTRAL POD ----------------------#
        # --------------- SPECTRAL POD ----------------------#
        # --------------- SPECTRAL POD ----------------------#
        # --------------- SPECTRAL POD ----------------------#
        CSD   = np.zeros((np.amax(qVar)+1,F,N,N),dtype=np.complex64)
        svec  = np.zeros((F,N,N),dtype=np.complex64)
        sval  = np.zeros((F,N),dtype=np.complex64)
        for iv in range(nVar):
            matLoad    = loadmat(outdir + '/%s.spectralCorrelation.mat' % (var(iv)))
            CSD[iv,:,:,:] = matLoad['CSD']

        for f in range(F):
            # Compute SPOD norm
            Cnorm = CSD[0,f,:,:] + CSD[1,f,:,:] + CSD[2,f,:,:]              # kinetic energy norm
            #Cnorm = CSD[3,f,:,:] + CSD[4,f,:,:] + CSD[5,f,:,:]              # enstrophy norm
            #Cnorm = CSD[6,f,:,:]                                            # pressure norm

            # compute SVD of normalized correlation matrix
            us,ss,vs = sp.svd(Cnorm)
            ss       = np.sqrt(ss)

            svec[f,:,:] = vs
            sval[f,:]   = ss

        savemat(outdir + '/spectralValues.mat',    {'Ss':sval})
        savemat(outdir + '/spectralVectors.mat',   {'Vs':svec})
        
        # Wrapping up everything and writing output files
        end_time = time.time()

        print('\n JOB FINISHED \n')
        print(' Total time: %.2f \n' % (end_time - start_time))

    comm.barrier()
