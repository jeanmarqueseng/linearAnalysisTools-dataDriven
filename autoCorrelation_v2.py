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
from mpi4py import MPI

from genGrid import *

from scipy.io import loadmat,savemat

import argparse
import time

# ------------- Read files in binFiles folder -----------#
# ------------- Read files in binFiles folder -----------#
# ------------- Read files in binFiles folder -----------#
# ------------- Read files in binFiles folder -----------#
def readbinfiles(varname,idx):
    # is this strange? yes. But this is the way CharLES does it, so keep it like that
    if (idx >= 10000000):
        filename = '%s.%08d.dat' % (varname,idx)
    elif (idx >= 1000000):
        filename = '%s.%07d.dat' % (varname,idx)
    else:
        filename = '%s.%06d.dat' % (varname,idx)

    qread      = np.fromfile(filename, np.float64)
    return qread

# ------------- Index-to-Variable Function --------------#
# ------------- Index-to-Variable Function --------------#
# ------------- Index-to-Variable Function --------------#                
# ------------- Index-to-Variable Function --------------#
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

    # ---------------------- INPUT ARGS ---------------------#
    # ---------------------- INPUT ARGS ---------------------#
    # ---------------------- INPUT ARGS ---------------------#
    # ---------------------- INPUT ARGS ---------------------#
    parser = argparse.ArgumentParser()
    parser.add_argument('-range',type=int, nargs=3,help='<istart> <istep> <iend>')
    parser.add_argument('-domain',type=int, nargs=2,help='<part> of <whole>')
    parser.add_argument('-opt',type=str,help='<meanflow> or <correlation>')
    parser.add_argument('-outputdir',type=str,help='Output folder for the binary files')
    args = parser.parse_args()

    if(args.range):
        frange = args.range
        files = np.arange(frange[0],frange[1]+frange[2],frange[1])
    else:
        raise RuntimeError('You must tell me the range of files you want: <istart> <iend> <istep>')

    if(args.domain):
        domain = args.domain
    else:
        domain = False

    if(args.opt):
        opt = args.opt
    else:
        opt = 'correlation'

    if(args.outputdir):
        outdir = args.outputdir
    else:
        outdir = './flowfield/'

    # ---------------------- Data Organization ---------------------#
    # ---------------------- Data Organization ---------------------#
    # ---------------------- Data Organization ---------------------#
    # ---------------------- Data Organization ---------------------#

    # number of variables we analyze: RHO, U-X, U-Y, U-Z, VORT-X, VORT-Y, VORT-Z, P in this order
    nVar = 8                                        
    # Size of autocorrelation snapshot matrix, it's the number of snapshots we have
    T = files.shape[0]

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
            makeWing(26.0,0.0,0.01,2,closed)     # in order <aoa>,<sweep>,<wingspan>,<# of nodes in z>,<closed or open wing tip>
        nno = genGrid(nodes, elems)
    comm.barrier()

    # ------------ Cell-Centered info (from postpro)----------------#
    # ------------ Cell-Centered info (from postpro)----------------#
    # ------------ Cell-Centered info (from postpro)----------------#
    # ------------ Cell-Centered info (from postpro)----------------#
    vol   = np.fromfile('binFiles/cellVolume.dat', np.float64)
    x_cv  = np.fromfile('binFiles/cellCentroid.dat', np.float64)
    ncv   = vol.shape[0]
    x_cv  = np.reshape(x_cv,(ncv,3))

    # --------- MPI organization per 'domain' and 'variable' -------#
    # --------- MPI organization per 'domain' and 'variable' -------#
    # --------- MPI organization per 'domain' and 'variable' -------#
    # --------- MPI organization per 'domain' and 'variable' -------#
    if (domain):
        part        = domain[0]
        whole       = domain[1]+1  # just to work with python zero-indexing
        minPart     = ncv // whole
        addCvtoPart = ncv % whole
        ncvPart     = np.zeros((whole,),dtype=np.int32)
        qIdx        = np.zeros((whole,2),dtype=np.int32)
        for r in range(whole):
            qIdx[r,0] = np.sum(ncvPart)
            if(r<addCvtoPart):
                ncvPart[r] = minPart + 1
            else:
                ncvPart[r] = minPart
            qIdx[r,1] = qIdx[r,0] + ncvPart[r]
        if (np.sum(ncvPart)!=ncv):
            raise RuntimeError('For some reason MPI distribution is not working well, please check in the code')
    else:
        part          = 0
        ncvPart       = np.zeros((1,),dtype=np.int32)    
        qIdx          = np.zeros((1,2),dtype=np.int32)
        ncvPart[part] = ncv
        qIdx[part,0]  = 0 
        qIdx[part,1]  = ncvPart[part]

    minRank      = nVar // size
    addVartoRank = nVar % size
    varRank      = np.zeros((size,),dtype=np.int32)
    varIdx       = np.zeros((size,2),dtype=np.int32)
    for r in range(size):
        varIdx[r,0] = np.sum(varRank)
        if(r<addVartoRank):
            varRank[r] = minRank + 1
        else:
            varRank[r] = minRank
        varIdx[r,1] = varIdx[r,0] + varRank[r]
    if (np.sum(varRank)!=nVar):
        raise RuntimeError('For some reason MPI distribution is not working well, please check in the code')

    # -------------------- SNAPSHOT MATRIX --------------------------#
    # -------------------- SNAPSHOT MATRIX --------------------------#
    # -------------------- SNAPSHOT MATRIX --------------------------#
    # -------------------- SNAPSHOT MATRIX --------------------------#
    if (opt == 'meanflow'):
        
        qread = np.zeros((ncv),dtype=np.float64)
        
        if (rank == 0): CGNS.write_unstructured_solndata_3d_cgns(nno,ncv,'CellCenter','flowfield/timeAveraged.cgns')
        comm.barrier()

        for iv in range(nVar):

            varname = 'binFiles/%s' % (var(iv))
            qm      = np.zeros((ncvPart[part]),dtype=np.float64)
            for idx in files:
                if (rank == 0): print('...reading %s variable, file: %d' % (var(iv), idx))
                qread = readbinfiles('binFiles/%s' % (var(iv)),idx)
                qm   += qread[np.arange(qIdx[part,0],qIdx[part,1])]
                
            # Mean flow
            qm = qm/T 
            comm.Gatherv(sendbuf=qm, recvbuf=(qread, ncvPart), root=0)      # Gather to plot on root node
            if (rank == 0): CGNS.write_unstructured_fielddata_3d_cgns(qread,iv,'flowfield/timeAveraged.cgns')
            comm.barrier()

    elif (opt == 'correlation'):

        C     = np.zeros((nVar,T,T),dtype=np.float32)
        qread = np.zeros((ncv),dtype=np.float64)
        q = np.zeros((ncvPart[part],T),dtype=np.float64)                       # initialize data matrix

        for iv in range(varIdx[rank][0],varIdx[rank][1]):

            varname = 'binFiles/%s' % (var(iv))
            tsIdx = 0                                                       # initialize timestep index
            for idx in files:
                qread      = readbinfiles('binFiles/%s' % (var(iv)),idx)
                q[:,tsIdx] = qread[np.arange(qIdx[part,0],qIdx[part,1])]
                tsIdx     += 1                                              # timestep index increment

            # Remove mean flow
            qm = np.mean(q,axis=1)                                          # Mean flow
            for i in range(T):
                q[:,i] = q[:,i] - qm

            # Compute autoCorrelation matrix C for each variable
            for i in range(T):
                if (rank == 0): print('...computing %s variable, row: %d' % (var(iv), i))
                for j in range(i,T):
                    C[iv,i,j] = np.dot(q[:,i]*vol[np.arange(qIdx[part,0],qIdx[part,1])],q[:,j])
                    if (j>i):
                        C[iv,j,i] = C[iv,i,j]


        C = comm.allreduce(C, op=MPI.SUM)
    
    # Wrapping up everything and writing output files
    end_time = time.time()

    if (rank == 0):
        if (opt == 'meanflow'): 
            CGNS.write_link_cgns('flowfield/timeAveraged.cgns','grid.cgns')
        else:
            for iv in range(8):
                if (domain):
                    outputFile = outdir + '/%s.%02d.autoCorrelation.mat' % (var(iv),part)
                else:
                    outputFile = outdir + '/%s.autoCorrelation.mat' % (var(iv))
                savemat(outputFile, {"C":C[iv,:,:]})

        print('\n JOB FINISHED \n')
        print(' Total time: %.2f \n' % (end_time - start_time))

