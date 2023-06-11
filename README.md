# structure-from-fes

### Requirements

The metadynamics run has to be 3-dimensional (exactly 2 collective variables (CVs) + energy). 

The scripts use numpy, os, sys, time, itertools, shapely, copy, multiprocessing and tqdm modules. These have to be installed before running the scripts.

The script has to be started in the parent directory of the MD-output directory.<br>
The parent directory has to include a FES-file named fes_{MD-output directory}_{other}.dat. The FES-file has to include the CV1, CV2 and ENERGY in the first three columns.<br>
The MD-output directory has to include a PDB-file and a COLVAR-file. If multiple PDB-files are present, the largest one will be used. The COLVAR-file has to be named COLVAR. Only one COLVAR-file is expected. The COLVAR-file has to include the STEPS, CV1 and CV2 in the first three columns.<br>
The COLVAR and PDB-file must have a similar MD-step length.

### Input parameters

Input parameters: {MD-output directory} {threshold/"auto"}.<br>
All input parameters must be specified.

The threshold-values correspond to the free energy specified in the FES-file. If "auto" is chosen, the lowest 2/7 of the total free energy span are considered as the energy of the minima.

### Created PDB-files

The minimum frames will be printed in (for each distinctive minimum) separate PDB-files inside the MD-output directory in a "minima" directory.<br>
Average CV1 and CV2 for each minimum will be printed in the TITLE line inside the PDB-files.<br>
Structure will be fixed in the minimum PDB-file relative to the first atom, avoiding jumps between frames.<br>
Periodic boundary conditions will be neglected.<br>
Only x,y,z, atom-type and atom-number columns will be copied from the trajectory PDB-file to the minimum PDB-file. Other columns will be filled with zeros.

### Usage advice

It is advised to use FES-histograms with sufficient bins, since the tolerance of the minimum frame identification and mapping of minima
depend on the bin quantity.<br>
The multiprocessing version shows significant use of memory. In case this behavior is problematic for the used machine it is advised to use the serial script or upgrade the memory of the machine.<br>
The scripts were tested on Windows and Linux. The output may lack certain visuals if used on Windows.
