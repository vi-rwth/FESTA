# structure-from-fes

### Requirements

The metadynamics run has to be 3-dimensional (exactly 2 collective variables (CVs) + energy). 

The scripts use numpy, os, sys, time, itertools, shapely, copy, multiprocessing, MDAnalysis, matplotlib and tqdm modules. These have to be installed before running the scripts.

The script has to be started in the parent directory of the MD-output directory.<br>
The MD-output directory has to include a FES-file, a trajectory-file and a COLVAR-file. The COLVAR-file has to be named COLVAR and the FES-file must be named "fes.dat". The FES-file has to include the CV1, CV2 and ENERGY in the first three columns. The COLVAR-file has to include the STEPS, CV1 and CV2 in the first three columns. Only one COLVAR-file is expected.<br>
The COLVAR and trajectory-file must have a similar MD-step length.

Most trajectory-formats supported by MDAnalysis are supported here as well. For a complete list see:<br>
https://docs.mdanalysis.org/stable/documentation_pages/coordinates/init.html#id2 <br>
Trajectory-formats must have multiframe capabilities. Not supported formats are: PDBQT, PQR, MMTF, NAMDBIN, IN, FHAIMS, COOR, CRD.<br>
If the trajectory-file does not offer topological information, a separate topology-file has to be provided.

### Input parameters

Input parameters: {MD-output directory} {threshold/"auto"} {trajectory-file} OPTIONAL:{topology-file}.<br>

The threshold-values correspond to the free energy specified in the FES-file. If "auto" is chosen, the lowest 2/7 of the total free energy span are considered as the energy of the minima.

### Created files

The minimum frames will be printed in (for each distinctive minimum) separate trajectory-files inside the MD-output directory in a "minima" directory. The format will be chosen based on the trajectory-file format. The minima directory contains a FES-visualization PNG-file with a white outline showing the minima-areas.<br>
Average CV1 and CV2 for each minimum will be printed in a separate "min_overview.txt"-file inside the "minima" directory.<br>
ONLY PDB: Structure will be fixed in the minimum trajectory-file relative to the first atom, avoiding jumps between frames.<br>
ONLY PDB: Only x,y,z, atom-type and atom-number columns will be copied from the trajectory PDB-file to the minimum PDB-file. Other columns will be filled with zeros.<br>
Periodic boundary conditions will be neglected.<br>

### Usage advice

It is advised to use FES-histograms with sufficient bins, since the tolerance of the minimum frame identification and mapping of minima depend on the bin quantity.<br>
The multiprocessing version might show significant use of memory. In case this behavior is problematic for the used machine it is advised to use the serial script or upgrade the memory of the machine.<br>
The scripts were tested on Windows and Linux. The output may lack certain visuals if used on Windows.
