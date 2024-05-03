# Metadynamics Structure Extraction Tool (MSET)

The goal of this project is to provide a helpful tool for fast and accurate extraction of trajectory frames from 2D metadynamics simulations, by identifying minima based on a free-energy threshold set by the user.

## Before running the script

### Requirements

- Metadynamics simulation with 2 collective variables
- single output directory with:
    +  Free-Energy-Surface (FES) histogram file
    +  COLVAR file
    +  Trajectory file
    +  Topology file (_if no topology information in trajectory file_)

### Formats

- PLUMED-like appearence expected for COLVAR- and FES-files <br>
- Multiframe capable trajectory, supported by MDAnalysis. Complete list: <br>
  https://docs.mdanalysis.org/stable/documentation_pages/coordinates/init.html#id2 <br>
  **Not Supported:** PDBQT, PQR, MMTF, NAMDBIN, IN, FHAIMS, COOR, CRD, GRO

### Installation

The following libraries have to be manually installed (e.g. using `pip`):

- numpy
- psutil
- tqdm
- matplotlib
- shapely
- MDAnalysis

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
