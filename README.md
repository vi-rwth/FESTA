# structure-from-fes

Program uses numpy, os, sys, time, itertools, shapely, copy, muliprocessing and tqdm modules. These have to be installed before running the program.
Program has to be started in parent directory of MD-output directory.
Parent directory has to include FES-file named fes_<MD-output directory>_<other>.dat
MD-output directory has to include PDB-file and COLVAR-file.
Only x,y,z, atom-type and atom-number columns will be copied from trajectory PDB-file to minimum PDB-file. Other columns will be filled with zeros.

Input parameters: <MD-output directory> <threshold/auto> <Frozen_atom_number>
All input parameters have to be specified.

Minima frames will be printed in (for each distinctive minimum) separate PDB-files inside MD-output directory in "minima" directory.
Average a and b for each minimum will be printed in the TITLE line inside the PDB-files.
Frozen atom specified in input will be used to fix this atom in the minimum PDB-file, avoiding jumps between frames.
Periodic boundary conditions will be neglected.

It is advised to use FES-histograms with a sufficient quantity of bins, since the tolerance of the minimum frame identification and mapping of minima
depend on the bin quantity.
The multiprocessing version might show significant use of memory. In case this behavior is problematic for the used machine it is advised to use the
serial version of the program.
