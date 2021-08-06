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
    parser.add_argument('-range',type=int, nargs=3,help='<istart> <iend> <istep>')
    parser.add_argument('-modes',type=int, nargs=3,help='<istart> <iend> <istep>')
    parser.add_argument('-outputdir',type=str,help='Output folder for the binary files')
    args = parser.parse_args()

    if(args.range):
        frange = args.range
        files = np.arange(frange[0],frange[1]+frange[2],frange[2])
    else:
        raise RuntimeError('You must tell me the range of FILES you want: <istart> <iend> <istep>')

    if(args.modes):
        mrange = args.modes
        mrange = np.arange(mrange[0],mrange[1]+mrange[2],mrange[2])
    else:
        raise RuntimeError('You must tell me the range of MODES you want: <istart> <iend> <istep>')

    if(args.outputdir):
        outdir = args.outputdir
    else:
        outdir = './'

    # number of snapshots and range of variables in DMD+POD modes
    T = files.shape[0]
    M = mrange.shape[0]
    nVar = np.arange(1,4,1)                 # prints only velocities

    # Load singular values and vectors, eigenvalues and vectors
    s      = np.loadtxt(outdir + '/singularValues.dat',delimiter=' ')
    u      = np.loadtxt(outdir + '/singularVectors.dat',delimiter=' ')
    C      = np.loadtxt(outdir + '/correlationMatrix.dat',delimiter=' ')
    mu     = np.loadtxt(outdir + '/eigenValues.dat',dtype=np.complex128)
    v      = np.loadtxt(outdir + '/eigenVectors.dat',dtype=np.complex128)

    # Get important grid information
    index_file,nbases = CGNS.open_file_read('flowfield/grid.cgns')
    nno, ncv, naux = CGNS.zone_size_read(index_file, 1, 1)
    CGNS.close_file(index_file)
    
    # MPI: ompute number of CVs per rank (MPI)
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

    # Allocate correlation matrix and read variable
    qread = np.zeros((ncv),dtype=np.float64)
    
    if (rank == 0): 
        for i in mrange:
            CGNS.write_unstructured_solndata_3d_cgns(nno,ncv,'CellCenter','flowfield/podmodes_k%03d.cgns' % i)
            CGNS.write_unstructured_solndata_3d_cgns(nno,ncv,'CellCenter','flowfield/dmdmodes_k%03d.cgns' % i)
    comm.barrier()

    # Stack variables in matrix
    for iv in nVar:

        varname = 'binFiles/%s' % (var(iv))

        tsIdx = 0                                                       # initialize timestep index
        q     = np.zeros((ncvRank[rank],T),dtype=np.float64)            # initialize data matrix
        qpod  = np.zeros((ncvRank[rank],M),dtype=np.float64)            # initialize POD matrix
        qdmd  = np.zeros((ncvRank[rank],M),dtype=np.float64)            # initialize DMD matrix
        
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
        for i in range(T):
            q[:,i] = q[:,i] - qm

        # Write DMD modes
        qdmd = q[:,1:]@u@v@np.diag(1/s)
        for j in mrange:
            qsend = np.real(qdmd[:,j]).ravel()
            comm.Gatherv(sendbuf=qsend, recvbuf=(qread, ncvRank), root=0)      # Gather to plot on root node
            if (rank == 0): 
                CGNS.write_unstructured_fielddata_3d_cgns(qread,iv,'flowfield/podmodes_k%03d.cgns' % j)
        comm.barrier()

        # Write POD modes
        qpod = q[:,0:-1]@u@np.diag(1/s)
        for j in mrange:
            qsend = np.real(qpod[:,j]).ravel()
            comm.Gatherv(sendbuf=qsend, recvbuf=(qread, ncvRank), root=0)      # Gather to plot on root node
            if (rank == 0): 
                CGNS.write_unstructured_fielddata_3d_cgns(qread,iv,'flowfield/podmodes_k%03d.cgns' % j)
        comm.barrier()


    # Wrapping up everything and writing output files
    end_time = time.time()

    if (rank == 0):
        for j in mrange:
            CGNS.write_link_cgns('flowfield/podmodes_k%03d.cgns' % j,'grid.cgns')
            CGNS.write_link_cgns('flowfield/dmdmodes_k%03d.cgns' % j,'grid.cgns')

        print('\n JOB FINISHED \n')
        print(' Total time: %.2f \n' % (end_time - start_time))

