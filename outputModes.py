# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#
# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#
# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#
# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#

# This code generates the spatial modes for POD and DMD

# ----- Department of Mechanical and Aerospace Engineering
# ----- University of California, Los Angeles
# ----- Author: Jean Helder Marques Ribeiro
# ----- email:  jeanmarques@g.ucla.edu
# ----- Date: August/2021import os,sys,getopt
import os,sys,getopt
sys.path.append('./labTools')

import numpy as np
from numpy.fft import fft, ifft

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
    parser.add_argument('-range',type=int, nargs=3,help='<istart> <istep> <iend>')
    parser.add_argument('-modes',type=int, nargs=3,help='<istart> <istep> <iend>')
    parser.add_argument('-saveAll',type=bool,help='Save all SPOD modes')
    parser.add_argument('-outputdir',type=str,help='Output folder for the binary files')
    args = parser.parse_args()

    if(args.range):
        frange = args.range
        files = np.arange(frange[0],frange[1]+frange[2],frange[1])
    else:
        raise RuntimeError('You must tell me the range of FILES you want: <istart> <iend> <istep>')

    if(args.modes):
        mrange = args.modes
        mrange = np.arange(mrange[0],mrange[1]+mrange[2],mrange[1])
    else:
        raise RuntimeError('You must tell me the range of MODES you want: <istart> <iend> <istep>')

    if(args.outputdir):
        outdir = args.outputdir
    else:
        outdir = './flowfield/'

    if(args.saveAll):
        saveflag = args.saveAll
    else:
        saveflag = False

    # ---------------------- Data Organization ---------------------#
    # ---------------------- Data Organization ---------------------#
    # ---------------------- Data Organization ---------------------#
    # ---------------------- Data Organization ---------------------#
    # number of snapshots and range of variables in DMD+POD modes
    T = files.shape[0]
    M = mrange.shape[0]
    qVar = np.arange(1,4,1)                 # prints only velocities

    matLoad    = loadmat(outdir + '/%s.spectralCorrelation.mat' % (var(0)))
    F = matLoad['CSD'].shape[0]
    N = matLoad['CSD'].shape[1]
    nblk  = (N+1)//2
    B     = T//nblk

    # Indexing FFT blocks for Spectral POD (Schmidt and Colonius)
    nblkAll = 2*nblk - 1
    blkPart = np.zeros((nblkAll,),dtype=np.int32)
    blkIdx  = np.zeros((nblkAll,2),dtype=np.int32)
    blkSize = T // nblk
    for n in range(0,nblkAll,2):                # loop for main blocks
        blkIdx[n,0] = np.sum(blkPart)
        blkPart[n]  = blkSize
        blkIdx[n,1] = blkIdx[n,0] + blkSize
    for n in range(1,nblkAll,2):                # loop for overlap blocks
        blkIdx[n,0] = blkIdx[n-1,1] // 2
        blkPart[n]  = blkSize
        blkIdx[n,1] = blkIdx[n,0] + blkSize


    # Load singular values and vectors, eigenvalues and vectors
    matLoad = loadmat(outdir + '/singularValues.mat');       ss = matLoad['S']
    matLoad = loadmat(outdir + '/singularVectors.mat');      u  = matLoad['U']
    matLoad = loadmat(outdir + '/correlationMatrix.mat');    C  = matLoad['C']
    matLoad = loadmat(outdir + '/eigenValues.mat');          mu = matLoad['mu']
    matLoad = loadmat(outdir + '/eigenVectors.mat');         v  = matLoad['V']
    s = np.zeros((T-1,T-1))
    for i in range(T-1):
        s[i,i] = 1./ss[0][i]

    # ----------- Get info from GRID.CGNS -------------------------#
    # ----------- Get info from GRID.CGNS -------------------------#
    # ----------- Get info from GRID.CGNS -------------------------#
    # ----------- Get info from GRID.CGNS -------------------------#
    index_file,nbases = CGNS.open_file_read('flowfield/grid.cgns')
    nno, ncv, naux = CGNS.zone_size_read(index_file, 1, 1)
    CGNS.close_file(index_file)
    
    # --------- MPI organization per 'domain' only ----------------#
    # --------- MPI organization per 'domain' only ----------------#
    # --------- MPI organization per 'domain' only ----------------#
    # --------- MPI organization per 'domain' only ----------------#
    minRank     = ncv//size
    addCvtoRank = ncv%size
    ncvRank     = np.zeros((size,),dtype=np.int32)
    qIdx        = np.zeros((size,2),dtype=np.int32)
    for r in range(size):
        qIdx[r,0] = np.sum(ncvRank)
        if(r<addCvtoRank):
            ncvRank[r] = minRank + 1
        else:
            ncvRank[r] = minRank
        qIdx[r,1] = qIdx[r,0] + ncvRank[r]
    if (np.sum(ncvRank)!=ncv):
        raise RuntimeError('For some reason MPI distribution is not working well, please check in the code')

    if (rank == 0): 
        for i in mrange:
            CGNS.write_unstructured_solndata_3d_cgns(nno,ncv,'CellCenter','flowfield/podmodes_k%03d.cgns' % i)
            CGNS.write_unstructured_solndata_3d_cgns(nno,ncv,'CellCenter','flowfield/dmdmodes_k%03d.cgns' % i)
        for f in range(F):
            if saveflag:
                for n in range(N):
                    CGNS.write_unstructured_solndata_3d_cgns(nno,ncv,'CellCenter','flowfield/spodmodes_f%05d_k%03d.cgns' % (f,n))
            else:
                CGNS.write_unstructured_solndata_3d_cgns(nno,ncv,'CellCenter','flowfield/spodmodes_f%05d_k%03d.cgns' % (f,0))
                
                
    comm.barrier()

    # -------------- ESTIMATED RAM MEMORY PER PROCESS -------------------------#
    # -------------- ESTIMATED RAM MEMORY PER PROCESS -------------------------#
    # -------------- ESTIMATED RAM MEMORY PER PROCESS -------------------------#
    # -------------- ESTIMATED RAM MEMORY PER PROCESS -------------------------#
    if (rank == 0): 
        memall = ( (ncv)*8. +
                   (ncvRank[rank])*16. +
                   (ncvRank[rank]*T)*8. +
                   (ncvRank[rank]*B)*16. +
                   (ncvRank[rank]*F*N)*16. ) / 1024 / 1024 / 1024 
        print('\n... ESTIMATED RAM MEMORY TO BE ALLOCATED (per process): %.3f Gb' % (memall))

    # Allocate q data matrices for POD/DMD modes
    qread = np.zeros((ncv),dtype=np.float64)
    qmod  = np.zeros((ncvRank[rank]),  dtype=np.complex128)            # initialize POD/DMD matrix
    q     = np.zeros((ncvRank[rank],T),dtype=np.float64)            # initialize data matrix
    
    comm.barrier()

    for iv in qVar:
        # -------------- COMPUTATION OF POD/DMD MODES -------------------------#
        # -------------- COMPUTATION OF POD/DMD MODES -------------------------#
        # -------------- COMPUTATION OF POD/DMD MODES -------------------------#
        # -------------- COMPUTATION OF POD/DMD MODES -------------------------#
        if (rank == 0): print('\n... STANDARD POD and DMD:')
        varname = 'binFiles/%s' % (var(iv))
        tsIdx = 0                                                   # initialize timestep index
        for idx in files:
            if (rank == 0): print('... ... reading %s variable, file: %d' % (var(iv), idx))
            qread      = readbinfiles('binFiles/%s' % (var(iv)),idx)
            q[:,tsIdx] = qread[np.arange(qIdx[rank,0],qIdx[rank,1])]
            tsIdx      += 1    

        # Remove mean flow
        qm = np.mean(q,axis=1)                                      # Mean flow
        for i in range(T):
            q[:,i] = q[:,i] - qm

        # Write DMD modes
        for j in mrange:
            qmod  = ( q[:,1:]@u@s@v[:,j] ) / mu[0,j]
            qsend = np.real(qmod).ravel()
            comm.Gatherv(sendbuf=qsend, recvbuf=(qread, ncvRank), root=0)      # Gather to plot on root node
            if (rank == 0): 
                print('... ... ... computing %s variable, DMD mode: %d' % (var(iv), j))
                CGNS.write_unstructured_fielddata_3d_cgns(qread/np.amax(np.abs(qread)),iv,'flowfield/dmdmodes_k%03d.cgns' % j)
        
        # Write POD modes
        for j in mrange:
            qmod  = ( q[:,0:-1]@u[:,j] )  * s[j,j]
            qsend = np.real(qmod).ravel()
            comm.Gatherv(sendbuf=qsend, recvbuf=(qread, ncvRank), root=0)      # Gather to plot on root node
            if (rank == 0): 
                print('... ... ... computing %s variable, POD mode: %d' % (var(iv), j))
                CGNS.write_unstructured_fielddata_3d_cgns(qread/np.amax(np.abs(qread)),iv,'flowfield/podmodes_k%03d.cgns' % j)

        # -------------- COMPUTATION OF SPOD MODES -------------------------#
        # -------------- COMPUTATION OF SPOD MODES -------------------------#
        # -------------- COMPUTATION OF SPOD MODES -------------------------#
        # -------------- COMPUTATION OF SPOD MODES -------------------------#
 
        # Allocate q data matrices for SPOD modes
        if (rank == 0): print('\n... SPECTRAL POD:')
        qf    = np.zeros((ncvRank[rank],B),dtype=np.complex128)         # read  FFT output
        qfft  = np.zeros((ncvRank[rank],F,N),dtype=np.complex128)       # store FFT-data matrix for SPOD
    
        # Load singular values and vectors
        matLoad = loadmat(outdir + '/spectralValues.mat');       sval = matLoad['Ss']
        matLoad = loadmat(outdir + '/spectralVectors.mat');      svec = matLoad['Vs']    
        smval = np.zeros((F,N,N))
        for f in range(F):
            for i in range(N):
                smval[f,i,i] = 1./np.abs(sval[f,i])

        for n in range(N):
            qf   = fft(q[:,np.arange(blkIdx[n,0],blkIdx[n,1])], axis=1)
            qfft[:,:,n] = qf[:,:F]

        for f in range(1,F):
            if saveflag:
                for j in range(N):
                    qmod = ( qfft[:,f,:]@svec[f,:,j] )  * smval[f,j,j]
                    qsend = np.real(qmod).ravel()
                    comm.Gatherv(sendbuf=qsend, recvbuf=(qread, ncvRank), root=0)      # Gather to plot on root node
                    if (rank == 0): 
                        print('... ... ... computing %s variable, SPOD mode %d, frequency %d out of %d' % (var(iv), j, f, F-1))
                        CGNS.write_unstructured_fielddata_3d_cgns(qread/np.amax(np.abs(qread)),iv,'flowfield/spodmodes_f%05d_k%03d.cgns' % (f,j))            

            else:
                qmod = ( qfft[:,f,:]@svec[f,:,0] )  * smval[f,0,0]
                qsend = np.real(qmod).ravel()
                comm.Gatherv(sendbuf=qsend, recvbuf=(qread, ncvRank), root=0)      # Gather to plot on root node
                if (rank == 0): 
                    print('... ... ... computing %s variable, SPOD mode: %d, frequency %d out of %d' % (var(iv), 0, f, F-1))
                    CGNS.write_unstructured_fielddata_3d_cgns(qread/np.amax(np.abs(qread)),iv,'flowfield/spodmodes_f%05d_k%03d.cgns' % (f,0))            


    # Wrapping up everything and writing output files
    end_time = time.time()

    if (rank == 0):
        for j in mrange:
            CGNS.write_link_cgns('flowfield/podmodes_k%03d.cgns' % j,'grid.cgns')
            CGNS.write_link_cgns('flowfield/dmdmodes_k%03d.cgns' % j,'grid.cgns')
        for f in range(F):
            if saveflag:
                for n in range(N):
                    CGNS.write_link_cgns('flowfield/spodmodes_f%05d_k%03d.cgns' % (f,n),'grid.cgns')
            else:
                CGNS.write_link_cgns('flowfield/spodmodes_f%05d_k%03d.cgns' % (f,0),'grid.cgns')

        print('\n JOB FINISHED \n')
        print(' Total time: %.2f \n' % (end_time - start_time))

