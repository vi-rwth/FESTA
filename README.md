# Metadynamics Structure Extraction Tool (MSET)

### Requirements

The metadynamics run must be 3-dimensional (exactly 2 collective variables (CVs) + energy). 

The scripts use numpy, os, argparse, shutil, psutil, time, itertools, shapely, copy, multiprocessing (only for "parallel" script), MDAnalysis, matplotlib, operator and tqdm modules. These have to be installed before running the scripts.

The MD-output directory must include a FES-file, a trajectory-file and a COLVAR-file. The COLVAR- and FES-files are expected to be created using PLUMED or have PLUMED-like appearance.
The COLVAR- and trajectory-file must have the same MD-step length.

Most trajectory-formats supported by MDAnalysis are supported here as well. For a complete list see:<br>
https://docs.mdanalysis.org/stable/documentation_pages/coordinates/init.html#id2 <br>
Trajectory-formats must have multiframe capabilities. If the trajectory-file does not offer topological information, a separate topology-file must be provided.<br>
Verified NOT SUPPORTED formats are: PDBQT, PQR, MMTF, NAMDBIN, IN, FHAIMS, COOR, CRD, GRO.<br>
Verified SUPPORTED formats are: ARC, DCD, ENT, GMS, H5MD, LAMMPS, MOL2, NC, NCDF, PDB, TRJ, TRR, TRZ, XTC, XYZ.<br>

### Input parameters

Use "-h" command to access this explanation from the script.<br>
<pre>
"-traj"    --- MD trajectory-file name in the MD-output-directory. Format is also used for output-files.<br>
               REQUIRED FOR RUNNING SCRIPT<br>
"-md"      --- MD-output-directory path.<br>
               DEFAULT: Current directory path.<br>
"-thresh"  --- Specifies threshold for assigning. Input value must correspond with values in the FES-file.<br>
               DEFAULT: Lowest 1/12 of the total energy span.<br>
"-topo"    --- MD topology-file name in the MD-output-directory, if trajectory-file does not specify topology.<br>
               DEFAULT: None.<br>
"-fes"     --- FES-file name in the MD-output-directory.<br>
               DEFAULT: "fes.dat".<br>
"-colv"    --- COLVAR-file in the MD-output-directory.<br> 
               DEFAULT: "COLVAR".<br>
"-png"     --- Specifies whether a PNG-visualization of the FES should be created. Expects True/False.<br>
               DEFAULT: True.<br>
"-nopbc"   --- Suppresses the automatic periodicity search (triggered when the minima touch the edges). Expects True/False.<br>
               DEFAULT: False.<br>
"-mindist" --- Smallest allowed distance, at which areas are considered separate minima (unit: bins of FES-histogram). Must be larger than 1.<br>
               DEFAULT: 10.<br>
"-stride"  --- Reads only every n-th frame of trajectory.<br>
               DEFAULT: 1.<br>
</pre>
### Created files

The minimum frames will be printed in (for each distinctive minimum) separate trajectory-files inside the MD-output directory in a "minima" directory. The format will be chosen based on the trajectory-file format. If not declared otherwise, the minima directory will contain a FES-visualization PNG-file with a white outline showing the minima-areas.<br>
Average CV1 and CV2 for each minimum will be printed in a separate "min_overview.txt"-file inside the "minima" directory.<br>
Periodicity will be detected and considered automatically, if not declared otherwise.<br>
ONLY (CP2K) PDB: Structure will be fixed in the minimum trajectory-file relative to the first atom, avoiding jumps between frames.<br>
ONLY (CP2K) PDB: Only x,y,z, atom-type and atom-number columns will be copied from the trajectory PDB-file to the minimum PDB-file. Other columns will be filled with zeros.<br>

### Usage advice

It is advised to use FES-histograms with sufficient bins, since the tolerance of the minimum frame identification and mapping of minima depend on the bin quantity.<br>
The scripts were tested on Windows and Linux. The output may lack certain visuals if used on Windows.
