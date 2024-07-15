<a href="url"><img src="https://github.com/vi-rwth/FESTA/assets/133004724/b2c86b88-cc14-4ff3-af7b-c8e693cea11c" align="middle" height="400" width="1.562*400" ></a>

# Free Energy Surface Trajectory Analysis (FESTA)

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

The following libraries must be manually installed (e.g. using `pip`):

- numpy
- tqdm
- matplotlib
- shapely
- MDAnalysis

## Running the script

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

- directory "minima" with:
   + trajectory files (input format or XYZ) for distinctive minima
   + overview file with average CVs for each minimum
   + FES PNG-file with white outline highlighting minima areas (_if not declared otherwise_)

### Further information

- Periodicity will be detected and considered automatically. (_if not declared otherwise_)
- PDB-files created with `CP2K` cannot be opened with MDAnalysis and will therefore be handled by a separate reader.

## Usage advice

It is advised to use FES-histograms with sufficient bins, since the tolerance of the minimum frame identification and mapping of minima depend on the bin quantity.<br>
The script were tested on Windows and Linux. Performance can be significantly worse on Windows.
