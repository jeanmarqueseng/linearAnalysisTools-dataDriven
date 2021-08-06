# Generate grid files
import os,sys,getopt
import numpy as np
import CGNS

def genGrid(nodes,elems):

    # Read binary files
    x_no  = np.fromfile(nodes, np.float64)

    # Read ASCII files
    no_conn = np.loadtxt(elems) 

    # Number of nodes and CVs
    nno        = x_no.shape[0]//3
    ncv        = no_conn.shape[0]
    cell_shape = no_conn.shape[1]

    # Reshape nodal Coordinates
    x_no  = np.reshape(x_no,(nno,3))

    # Build grid in CGNS
    CGNS.create_file_cgns('flowfield/grid.cgns',3)
    CGNS.write_unstructured_grid_3d_cgns(nno,ncv,x_no[:,0],x_no[:,1],x_no[:,2],no_conn.T,'flowfield/grid.cgns')

    return nno


def makeWing(aoa,sw,dz,no3dcells,closed):

    # Generate wing geometry file:
    NACA   = 15             # Chord Thickness (NACA00??)
    points = 1001              # Number of points on the airfoil
    dx     = 1.0/(points-1)
    j = 0
    x0 = np.arange(0,1,dx) 
    y0 = np.zeros(x0.shape)
    for j in range(x0.shape[0]):
        xi    = x0[j]/1.0 
        y0[j] = 5.0 *NACA/100 *( + 0.29690*xi**0.5 
                              - 0.12600*xi**1 
                              - 0.35160*xi**2 
                              + 0.28430*xi**3 
                              - 0.10360*xi**4 )

    x0 = x0-1
    Alpha = -aoa
    Xact  = [] 
    Yact  = [] 

    xrt = np.zeros(x0.shape)
    yrt = np.zeros(x0.shape)
    xrb = np.zeros(x0.shape)
    yrb = np.zeros(x0.shape)
    for i in range(j):
        xrt[i] =  + x0[i]
        yrt[i] =  + y0[i]
        
        if (np.abs(xrt[i] + 0.917) <= 0.01):
            Xact.append(xrt[i]) 
            Yact.append(yrt[i]) 
        
        xrb[i] =  + x0[i]
        yrb[i] =  - y0[i]

    XA  = np.concatenate((xrt,xrb), axis=0) 
    YA  = np.concatenate((yrt,yrb), axis=0)
    ZA  = np.zeros(XA.shape[0])

    XA  = np.concatenate((XA,XA), axis=0) 
    YA  = np.concatenate((YA,YA), axis=0)
    ZA  = np.concatenate((ZA,dz*no3dcells*np.ones(ZA.shape[0])),axis=0)

    XA = XA + 1

    for ino in range(XA.shape[0]):
        dummyx = XA[ino]*np.cos(np.deg2rad(Alpha)) - YA[ino]*np.sin(np.deg2rad(Alpha))
        dummyy = XA[ino]*np.sin(np.deg2rad(Alpha)) + YA[ino]*np.cos(np.deg2rad(Alpha))
        XA[ino] = dummyx 
        YA[ino] = dummyy

    # QUAD_4
    if (closed):
        ns = 4
    else:
        ns = 2

    # QUAD_4
    ELEM = np.zeros((ns*(x0.shape[0]-1),4))

    for i in range(x0.shape[0]-1):
        # Suction side
        ii = 0*(x0.shape[0]-1) + i
        ELEM[ii,0] = np.int(i)
        ELEM[ii,1] = np.int(i + 1)
        ELEM[ii,2] = np.int(2*(points-1) + i + 1)
        ELEM[ii,3] = np.int(2*(points-1) + i)
        # Pressure side
        ii = 1*(x0.shape[0]-1) + i
        ELEM[ii,0] = np.int((points-1) + i)
        ELEM[ii,1] = np.int((points-1) + i + 1)
        ELEM[ii,2] = np.int(3*(points-1) + i + 1)
        ELEM[ii,3] = np.int(3*(points-1) + i)
        if (closed):
            # Side surface
            ii = 2*(x0.shape[0]-1) + i
            ELEM[ii,0] = np.int(i)
            ELEM[ii,1] = np.int(i+1)
            ELEM[ii,2] = np.int((points-1) + i + 1)
            ELEM[ii,3] = np.int((points-1) + i)
            # Side surface
            ii = 3*(x0.shape[0]-1) + i
            ELEM[ii,0] = np.int(2*(points-1) + i)
            ELEM[ii,1] = np.int(2*(points-1) + i + 1)
            ELEM[ii,2] = np.int(3*(points-1) + i + 1)
            ELEM[ii,3] = np.int(3*(points-1) + i)

    ELEM += np.int(1)

    # Apply sweep angle:
    for ino in range(XA.shape[0]):
        if (ZA[ino] > 1.0e-10):
            XA[ino] += ZA[ino]*np.tan(np.deg2rad(sw))

    CGNS.create_file_cgns('flowfield/wing.cgns',2)
    CGNS.write_unstructured_grid_2d_cgns(XA.shape[0],ELEM.shape[0],XA,YA,ZA,ELEM.T,'flowfield/wing.cgns')

    return 0

