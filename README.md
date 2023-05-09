# structure-from-fes

The scripts use numpy, os, sys, time, itertools, shapely, copy, muliprocessing and tqdm modules. These have to be installed before running the scripts.

The script has to be started in the parent directory of the MD-output directory.
The Parent directory has to include a FES-file named fes_{MD-output directory}_{other}.dat.
The MD-output directory has to include a PDB-file and a COLVAR-file. If multiple PDB-files are present, the largest one will be used.
Only x,y,z, atom-type and atom-number columns will be copied from the trajectory PDB-file to the minimum PDB-file. Other columns will be filled with zeros.

Input parameters: {MD-output directory} {threshold/auto} {frozen atom number}.
All input parameters have to be specified.

The Minima frames will be printed in (for each distinctive minimum) separate PDB-files inside the MD-output directory in "minima" directory.
Average a and b for each minimum will be printed in the TITLE line inside the PDB-files.
Frozen atom specified in input will be used to fix this atom in the minimum PDB-file, avoiding jumps between frames.
Periodic boundary conditions will be neglected.

It is advised to use FES-histograms with sufficient bins, since the tolerance of the minimum frame identification and mapping of minima
depend on the bin quantity.
The multiprocessing version might show significant use of memory. In case this behavior is problematic for the used machine it is advised to use the
serial version of the program.
