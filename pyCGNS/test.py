import CGNS
import numpy as np

x_no = [0,0,1,1,2,2,0,0,1,1,2,2,0,1,1,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1]
x_no = np.reshape(x_no,(3,12)).T

no_conn = np.array([[1,2,3,4,7,8,9,10],[4,3,5,6,10,9,11,12]])

nno = x_no.shape[0]
ncv = no_conn.shape[0]
cell_shape = no_conn.shape[1]

q = [1,2]


CGNS.create_file_cgns('grid.cgns','3D')

CGNS.write_unstructured_grid_3d_cgns(1,nno,ncv,cell_shape,x_no[:,0],x_no[:,1],x_no[:,2],no_conn.T,'grid.cgns')

CGNS.write_unstructured_soln_3d_cgns(1,nno,ncv,cell_shape,q,'soln.cgns')

CGNS.write_link_cgns(1,'soln.cgns','grid.cgns')