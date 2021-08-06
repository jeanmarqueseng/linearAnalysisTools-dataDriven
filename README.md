# linearAnalysisTools-dataDriven 
A unstructured data-driven framework for Professor Kunihiko Taira Group at UCLA

'qsub installCGNS.sub' on Hoffman2 to build CGNS library for I/O.
copy pyCGNS/CGNS.so within the folder you are running the analysis and within labTools folder.

execute with:
mpiexe -n $NSLOTS autoCorrelation.py -options

for each .py file.

Code is divided in 3 parts:

autoCorrelation.py: computes snapshot matrix for each variable

modalAnalysis.py: computes POD singular values and vectors as well as DMD eigenvalues and vectors

outputMode.py: generates CGNS modes for visualization.

autoCorrelation.py and outputModes.py are parallelized in spatial domain via MPI.

autoCorrelation.py inputs:

  -range <istart> <iend> <istep>:  range of files (built with np.arange())
        <istart> first file
        <iend> last file
        <istart> step btw files
          
  -outputdir :
        folder to put correlation matrices for each variable
          
          
modalAnalysis.py inputs:

  -range <istart> <iend> <istep>: range of files (built with np.arange())
        <istart> first file
        <iend> last file
        <istart> step btw files
          
  -outputdir :
        folder to put correlation matrices for each variable
          

outputModes.py inputs:

  -range <istart> <iend> <istep>: range of files (built with np.arange())
        <istart> first file
        <iend> last file
        <istart> step btw files
  -modes <istart> <iend> <istep>: range of modes to output both POD and DMD (built with np.arange())
        <istart> first mode
        <iend> last mode
        <istart> step btw modes
          
  -outputdir :
        folder to put correlation matrices for each variable

 
