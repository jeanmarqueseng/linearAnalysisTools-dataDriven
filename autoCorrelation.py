# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#
# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#
# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#
# ---------- WELCOME TO JEAN'S DATA DRIVEN WORLD !!!! -------------#

# This code computes the correlation SNAPSHOT MATRIX for each input flow variable
# You can also use this to process any unstructured data from CharLES
# Hey, do not forget to install the pyCGNS Library, okay? Let me know if you have any problems with it

# ----- Department of Mechanical and Aerospace Engineering
# ----- University of California, Los Angeles
# ----- Author: Jean Helder Marques Ribeiro
# ----- email:  jeanmarques@g.ucla.edu
# ----- Date: August/2021

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
    parser.add_argument('-variable',type=int, nargs=2,help='<part> of <whole>')
    parser.add_argument('-opt',type=str,help='<meanflow> or <correlation>')
    parser.add_argument('-nblocks',type=int,help='number of FFT blocks')
    parser.add_argument('-outputdir',type=str,help='Output folder for the binary files')
    args = parser.parse_args()

    if(args.range):
        frange = args.range
        files = np.arange(frange[0],frange[1]+frange[2],frange[1])
    else:
        raise RuntimeError('You must tell me the range of files you want: <istart> <iend> <istep>')

    if(args.variable):
        variable = args.variable
    else:
        variable = False

    if(args.opt):
        opt = args.opt
    else:
        opt = 'correlation'

    if(args.nblocks):
        nblk = args.nblocks
    else:
        nblk = 2


    if(args.outputdir):
        outdir = args.outputdir
    else:
        outdir = './flowfield/'

    # ---------------------- Data Organization ---------------------#
    # ---------------------- Data Organization ---------------------#
    # ---------------------- Data Organization ---------------------#
    # ---------------------- Data Organization ---------------------#

    # number of variables we analyze: U-X, U-Y, U-Z, VORT-X, VORT-Y, VORT-Z, P, RHO in this order
    qVar = np.arange(0,3,1)
    nVar = qVar.shape[0]                                        
    # Size of autocorrelation snapshot matrix, it's the number of snapshots we have
    T = files.shape[0]
    # Size of SPOD blocks
    N = 2*nblk - 1  # allocate number of blocks
    B = T // nblk   # allocate sizes of block
    F = B // 2      # allocate sizes for frequencies
    
    # ------------- Grid files generation with CGNS ----------------#
    # ------------- Grid files generation with CGNS ----------------#
    # ------------- Grid files generation with CGNS ----------------#
    # ------------- Grid files generation with CGNS ----------------#
    if (rank == 0):
        print('...Generating grid files CGNS...')
        nodes = 'binFiles/nodalCoordinates.dat'
        elems = 'binFiles/elementConnectivity.dat'
        airfoil = True      # False if you do not want to creat the wing.cgns file
        closed  = True
        if (airfoil):
            makeWing(30.0,0.0,4.0,2,closed)     # in order <aoa>,<sweep>,<wingspan>,<# of nodes in z>,<closed or open wing tip>
        nno = genGrid(nodes, elems)
    comm.barrier()

    # ------------ Cell-Centered info (from postpro)----------------#
    # ------------ Cell-Centered info (from postpro)----------------#
    # ------------ Cell-Centered info (from postpro)----------------#
    # ------------ Cell-Centered info (from postpro)----------------#
    vol   = np.fromfile('binFiles/cellVolume.dat', np.float64)
    #x_cv  = np.fromfile('binFiles/cellCentroid.dat', np.float64)
    ncv   = vol.shape[0]
    #x_cv  = np.reshape(x_cv,(ncv,3))

    # --------- MPI organization per 'domain' and 'variable' -------#
    # --------- MPI organization per 'domain' and 'variable' -------#
    # --------- MPI organization per 'domain' and 'variable' -------#
    # --------- MPI organization per 'domain' and 'variable' -------#
    if (variable):
        part         = variable[0]
        whole        = variable[1]+1  # just to work with python zero-indexing
        minPart      = nVar // whole
        addVartoPart = nVar % whole
        varPart      = np.zeros((whole,),dtype=np.int32)
        varIdx       = np.zeros((whole,2),dtype=np.int32)
        for r in range(whole):
            varIdx[r,0] = np.sum(varPart)
            if(r<addVartoPart):
                varPart[r] = minPart + 1
            else:
                varPart[r] = minPart
            varIdx[r,1] = varIdx[r,0] + varPart[r]
        if (np.sum(varPart)!=nVar):
            raise RuntimeError('For some reason MPI distribution is not working well, please check in the code')
    else:
        part            = 0
        varPart         = np.zeros((1,),dtype=np.int32)    
        varIdx          = np.zeros((1,2),dtype=np.int32)
        varPart[part]   = nVar
        varIdx[part,0]  = 0 
        varIdx[part,1]  = varPart[part]

    minRank      = ncv // size
    addCvtoRank  = ncv % size
    ncvRank      = np.zeros((size,),dtype=np.int32)
    qIdx         = np.zeros((size,2),dtype=np.int32)
    for r in range(size):
        qIdx[r,0] = np.sum(ncvRank)
        if(r<addCvtoRank):
            ncvRank[r] = minRank + 1
        else:
            ncvRank[r] = minRank
        qIdx[r,1] = qIdx[r,0] + ncvRank[r]
    if (np.sum(ncvRank)!=ncv):
        raise RuntimeError('For some reason MPI distribution is not working well, please check in the code')

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

    # ------------------------- CORE -------------------------#
    # ------------------------- CORE -------------------------#
    # ------------------------- CORE -------------------------#
    # ------------------------- CORE -------------------------#
    if (opt == 'meanflow'):

        if (rank == 0): 
            memall = ( (ncv)*8. +
                       (ncvRank[rank]*T)*8. ) / 1024 / 1024 / 1024 
            print('\n... ESTIMATED RAM MEMORY TO BE ALLOCATED (per process): %.3f Gb' % (memall))

        
        qread = np.zeros((ncv),dtype=np.float64)
        if (rank == 0): CGNS.write_unstructured_solndata_3d_cgns(nno,ncv,'CellCenter','flowfield/timeAveraged.cgns')
        comm.barrier()

        for iv in qVar:

            # ------------------------- MEAN FLOW -------------------------#
            # ------------------------- MEAN FLOW -------------------------#
            # ------------------------- MEAN FLOW -------------------------#
            # ------------------------- MEAN FLOW -------------------------#
            if (rank == 0): print('\n... MEAN FLOW:')
            varname = 'binFiles/%s' % (var(iv))
            qm      = np.zeros((ncvRank[rank]),dtype=np.float64)
            for idx in files:
                if (rank == 0): print('... ... reading %s variable, file: %d' % (var(iv), idx))
                qread = readbinfiles('binFiles/%s' % (var(iv)),idx)
                qm   += qread[np.arange(qIdx[rank,0],qIdx[rank,1])]
                
            # Mean flow
            qm = qm/T 
            comm.Gatherv(sendbuf=qm, recvbuf=(qread, ncvRank), root=0)      # Gather to plot on root node
            if (rank == 0): CGNS.write_unstructured_fielddata_3d_cgns(qread,iv,'flowfield/timeAveraged.cgns')
            comm.barrier()

        if (rank == 0): CGNS.write_link_cgns('flowfield/timeAveraged.cgns','grid.cgns')


    elif (opt == 'correlation'):

        if (rank == 0): 
            memall = ( (T*T)*4. +
                       (F*N*N)*8. +
                       (ncv)*8. +
                       (ncvRank[rank]*T)*8. +
                       (ncvRank[rank]*B)*16. +
                       (ncvRank[rank]*F*N)*16. ) / 1024 / 1024 / 1024 
            print('\n... ESTIMATED RAM MEMORY TO BE ALLOCATED (per process): %.3f Gb' % (memall))

        # memory allocation
        C     = np.zeros((T,T),dtype=np.float32)
        CSD   = np.zeros((F,N,N),dtype=np.complex64)
        qread = np.zeros((ncv),dtype=np.float64)                                   # read  temporal data
        q     = np.zeros((ncvRank[rank],T),dtype=np.float64)                       # store temporal data matrix
        qf    = np.zeros((ncvRank[rank],B),dtype=np.complex128)                    # read  FFT output
        qfft  = np.zeros((ncvRank[rank],F,N),dtype=np.complex128)                  # store FFT-data matrix for SPOD

        for iv in range(varIdx[part,0],varIdx[part,1]):

            # -------------------- SNAPSHOT MATRIX --------------------------#
            # -------------------- SNAPSHOT MATRIX --------------------------#
            # -------------------- SNAPSHOT MATRIX --------------------------#
            # -------------------- SNAPSHOT MATRIX --------------------------#
            if (rank == 0): print('\n... SNAPSHOT MATRIX:')
            varname = 'binFiles/%s' % (var(qVar[iv]))
            tsIdx = 0                                                       # initialize timestep index
            for idx in files:
                qread      = readbinfiles('binFiles/%s' % (var(qVar[iv])),idx)
                q[:,tsIdx] = qread[np.arange(qIdx[rank,0],qIdx[rank,1])]
                tsIdx     += 1                                              # timestep index increment

            # Remove mean flow
            qm = np.mean(q,axis=1)                                          # Mean flow
            for i in range(T):
                q[:,i] = q[:,i] - qm
                q[:,i] = q[:,i]*np.sqrt(vol[np.arange(qIdx[rank,0],qIdx[rank,1])])

            # Compute autoCorrelation matrix C for each variable
            for i in range(T):
                if (rank == 0): print('... ... computing %s variable, row: %d' % (var(qVar[iv]), i))
                for j in range(i,T):
                    C[i,j] = np.dot(q[:,i],q[:,j])
                    if (j>i):
                        C[j,i] = C[i,j]

            C = comm.allreduce(C, op=MPI.SUM)
    
            if (rank == 0):
                outputFile = outdir + '/%s.autoCorrelation.mat' % (var(qVar[iv]))
                savemat(outputFile, {"C":C})
            comm.barrier()

            # -------------------- SPECTRAL POD --------------------------#
            # -------------------- SPECTRAL POD --------------------------#
            # -------------------- SPECTRAL POD --------------------------#
            # -------------------- SPECTRAL POD --------------------------#
            if (rank == 0): print('\n... SPECTRAL POD:')
            for n in range(N):
                qf   = fft(q[:,np.arange(blkIdx[n,0],blkIdx[n,1])], axis=1)
                qfft[:,:,n] = qf[:,:F]

            for f in range(F):
                for n in range(N):
                    qfft[:,f,n] = qfft[:,f,n]*np.sqrt(vol[np.arange(qIdx[rank,0],qIdx[rank,1])])

            # Compute Cross-Spectral Density matrix CS for each variable
            for f in range(F):
                if (rank == 0): print('... ... computing %s variable, frequency %d out of %d' % (var(qVar[iv]), f, F-1))
                for i in range(N):
                    for j in range(i,N):
                        CSD[f,i,j] = np.dot(qfft[:,f,i],qfft[:,f,j].conj())
                        if (j>i):
                            CSD[f,j,i] = CSD[f,i,j].conj().transpose()

            CSD = comm.allreduce(CSD, op=MPI.SUM)
    
            if (rank == 0):
                outputFile = outdir + '/%s.spectralCorrelation.mat' % (var(qVar[iv]))
                savemat(outputFile, {"CSD":CSD})
            comm.barrier()

    # Wrapping up everything
    end_time = time.time()
    if (rank == 0):
        print('\n JOB FINISHED \n')
        print(' Total time: %.2f \n' % (end_time - start_time))

