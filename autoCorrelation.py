# Preprocess to get the required matrices
import os,sys,getopt
sys.path.append('./labTools')

import numpy as np
from mpi4py import MPI

from genGrid import *

import argparse
import time

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
    parser.add_argument('-mpiopt',type=str,help='split MPI by <domain> or <variable>')
    parser.add_argument('-range',type=int, nargs=3,help='<istart> <iend> <istep>')
    parser.add_argument('-outputdir',type=str,help='Output folder for the binary files')
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

    if(args.mpiopt):
        mpiopt = args.mpiopt
    else:
        mpiopt = 'domain'

    # MPI distribution by 'variables' only works if nVar = size
    nVar = 8                                        
    if (size != nVar): mpiopt = 'domain'

    # Size of autocorrelation snapshot matrix
    T = files.shape[0]

    # Grid files generation with CGNS
    if (rank == 0):
        print('...Generating grid files CGNS...')
        nodes = 'binFiles/nodalCoordinates.dat'
        elems = 'binFiles/elementConnectivity.dat'
        airfoil = True      # False if you do not want to creat the wing.cgns file
        closed  = True
        if (airfoil):
            makeWing(26.0,0.0,0.01,2,closed)     # in seq. <aoa>,<sweep>,<wingspan>,<# of nodes in z>,<closed or open wing tip>
        nno = genGrid(nodes, elems)
    comm.barrier()

    # Read cell centered information
    vol   = np.fromfile('binFiles/cellVolume.dat', np.float64)
    x_cv  = np.fromfile('binFiles/cellCentroid.dat', np.float64)
    ncv   = vol.shape[0]
    x_cv  = np.reshape(x_cv,(ncv,3))

    # MPI distribution according to parallelization per 'domain' or 'variable'
    if (mpiopt == 'domain'):
        # Compute number of CVs per rank (MPI)
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
    elif (mpiopt == 'variable'):
        ncvRank     = ncv*np.ones((size,),dtype=np.int32)

    # Allocate correlation matrix and read variable
    C     = np.zeros((nVar,T,T),dtype=np.float32)
    qread = np.zeros((ncv),dtype=np.float64)
    
    if (rank == 0): CGNS.write_unstructured_solndata_3d_cgns(nno,ncv,'CellCenter','flowfield/timeAveraged.cgns')
    comm.barrier()

    if (mpiopt == 'domain'):
        # Stack variables in matrix
        for iv in range(nVar):

            varname = 'binFiles/%s' % (var(iv))

            tsIdx = 0                                                       # initialize timestep index
            q = np.zeros((ncvRank[rank],T),dtype=np.float64)                # initialize data matrix
            
            for idx in files:
                if (idx >= 10000000):
                    filename = '%s.%08d.dat' % (varname,idx)
                elif (idx >= 1000000):
                    filename = '%s.%07d.dat' % (varname,idx)
                else:
                    filename = '%s.%06d.dat' % (varname,idx)

                qread      = np.fromfile(filename, np.float64)
                q[:,tsIdx] = qread[np.arange(qIdx[rank,0],qIdx[rank,1])]
                tsIdx += 1                                                  # timestep index increment

            # Remove mean flow
            qm = np.mean(q,axis=1)                                          # Mean flow
            comm.Gatherv(sendbuf=qm, recvbuf=(qread, ncvRank), root=0)      # Gather to plot on root node
            if (rank == 0): CGNS.write_unstructured_fielddata_3d_cgns(qread,iv,'flowfield/timeAveraged.cgns')
            comm.barrier()

            for i in range(T):
                q[:,i] = q[:,i] - qm

            # Compute autoCorrelation matrix C for each variable
            for i in range(T):
                for j in range(i,T):
                    C[iv,i,j] = np.dot(q[:,i]*vol[np.arange(qIdx[rank,0],qIdx[rank,1])],q[:,j])
                    if (j>i):
                        C[iv,j,i] = C[iv,i,j]

            C[iv,:,:] = comm.allreduce(C[iv,:,:], op=MPI.SUM)
        
    elif (mpiopt == 'variable'):
        varname = 'binFiles/%s' % (var(rank))

        tsIdx = 0                                                       # initialize timestep index
        q = np.zeros((ncvRank[rank],T),dtype=np.float64)                # initialize data matrix

        for idx in files:
            if (idx >= 10000000):
                filename = '%s.%08d.dat' % (varname,idx)
            elif (idx >= 1000000):
                filename = '%s.%07d.dat' % (varname,idx)
            else:
                filename = '%s.%06d.dat' % (varname,idx)

            qread      = np.fromfile(filename, np.float64)
            q[:,tsIdx] = qread
            tsIdx += 1                                                  # timestep index increment

        # Remove mean flow
        qm = np.mean(q,axis=1)                                          # Mean flow
        CGNS.write_unstructured_field_3d_cgns(qm,rank,'flowfield/timeAveraged.cgns')

        for i in range(T):
            q[:,i] = q[:,i] - qm

        # Compute autoCorrelation matrix C for each variable
        for i in range(T):
            for j in range(i,T):
                C[rank,i,j] = np.dot(q[:,i]*vol,q[:,j])
                if (j>i):
                    C[rank,j,i] = C[rank,i,j]

        C = comm.allreduce(C, op=MPI.SUM)
    
    # Wrapping up everything and writing output files
    end_time = time.time()

    if (rank == 0):
        CGNS.write_link_cgns('flowfield/timeAveraged.cgns','grid.cgns')
        for iv in range(8):
            outputFile = outdir + '/%s.autoCorrelation.dat' % (var(iv))
            np.savetxt(outputFile, C[iv,:,:], delimiter=' ')

        print('\n JOB FINISHED \n')
        print(' Total time: %.2f \n' % (end_time - start_time))

