# Preprocess to get the required matrices
import os,sys,getopt
sys.path.append('./labTools')

import numpy as np
import scipy.linalg as sp
from mpi4py import MPI

from genGrid import *

import argparse
import time

def spod_filter(C,Nfilt):
    T     = C.shape[0]
    Cbar  = np.zeros((T,T))
    # Gaussian function
    f = np.exp(-np.linspace(-2.285,2.285,2*Nfilt+1)**2) 
    f = f/np.sum(f)
    if (Nfilt > 0):
        for i in range(T):
            for j in range(T):
                dummy = 0.
                for k in range(-Nfilt,Nfilt+1):
                    dummy += C[i+k,j+k]*f[k+Nfilt] 
                Cbar[i,j] = dummy
    else:
        Cbar = C
    return Cbar

def var(x):
    return {
        0: 'rho',
        1: 'ux',
        2: 'uy',
        3: 'uz',
        4: 'vx',
        5: 'vy',
        6: 'vz',
        7: 'p',
    }[x]

if __name__ == '__main__':

    # communications for MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('-range',type=int, nargs=3,help='<istart> <iend> <istep>')
    parser.add_argument('-outputdir',type=str,help='Output folder for the binary files')
    parser.add_argument('-spod',type=int,help='Sieber-POD filter size: 0 (no filter), # of snapshots (DFT)')
    parser.add_argument('-dt',type=float,help='Time step between data files')
    args = parser.parse_args()

    if(args.range):
        frange = args.range
        files = np.arange(frange[0],frange[1]+frange[2],frange[2])
    else:
        raise RuntimeError('You must tell me the range of files you want: <istart> <iend> <istep>')

    if(args.outputdir):
        outdir = args.outputdir
    else:
        outdir = './'

    if(args.dt):
        dt = args.dt
    else:
        dt = 0.1

    if(args.spod):
        spod = args.spod
    else:
        spod = 0
    
    if (rank == 0):
        nVar = 8
        T    = files.shape[0]
        if (spod > T): spod = T-1

        # Read correlation matrix and read variable
        C     = np.zeros((nVar,T,T),dtype=np.float32)
        for iv in range(nVar):
            filename = outdir + '/%s.autoCorrelation.dat' % (var(iv))
            C[iv,:,:] = np.loadtxt(filename,delimiter=' ')

        # Compute POD/DMD norm
        Cnorm = C[1,:,:] + C[2,:,:] + C[3,:,:]              # kinetic energy norm
        #Cnorm = C[4,:,:] + C[5,:,:] + C[6,:,:]              # enstrophy norm
        #Cnorm = C[7,:,:]              # pressure norm

        # apply SPOD filter (SPOD by Sieber et al. (JFM, 2016))
        Cnorm    = spod_filter(Cnorm,spod)
        cc       = Cnorm[0:-1,0:-1]                         # Correlation matrix related to X
        cp       = Cnorm[0:-1,1:]                           # Correlation matrix related to X'
        
        # compute SVD of normalized correlation matrix
        uc,sc,vc = sp.svd(cc)
        sc       = np.sqrt(sc)

        np.savetxt(outdir + '/singularValues.dat', sc, delimiter=' ')
        np.savetxt(outdir + '/singularVectors.dat', uc, delimiter=' ')
        np.savetxt(outdir + '/correlationMatrix.dat', cc, delimiter=' ')
        
        # compute DMD modes (check J. Tu Thesis if you want, but really, it's just this)
        # I did not truncate the POD modes, but you can do that later if you want just by doing this:
        # trunc = 5
        # sc = sc[0:trunc]
        # uc = uc[0:trunc,:]
        # vc = vc[0:trunc,:]
        # Compute transition matrix A
        A   =  np.diag(1/sc)@vc@cp@vc.transpose().conj()@np.diag(1/sc)
        
        mu,v = sp.eig(A) 

        np.savetxt(outdir + '/eigenValues.dat', mu)
        np.savetxt(outdir + '/eigenVectors.dat', v)
        
        # Wrapping up everything and writing output files
        end_time = time.time()

        print('\n JOB FINISHED \n')
        print(' Total time: %.2f \n' % (end_time - start_time))

    comm.barrier()