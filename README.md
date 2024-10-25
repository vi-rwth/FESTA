# Free Energy Surface Trajectory Analysis (FESTA)
>by Valentin Istomin and GiovanniMaria Piccini

<p align="center">
  <img width="1.562*300" height="300" src="https://github.com/vi-rwth/FESTA/assets/133004724/b2c86b88-cc14-4ff3-af7b-c8e693cea11c">
</p>

The goal of this project is to provide a helpful tool for fast and accurate extraction of trajectory frames from 2D CV-based enhanced-sampling molecular dynamics simulations, by identifying minima via a free-energy threshold set by the user.

## Before running the script
### Requirements
- Metadynamics simulation with 2 collective variables (CVs)
- single output directory with:
    +  Free-Energy-Surface (FES) histogram file \
       Generable using PLUMED `sum_hills` module: https://www.plumed.org/doc-v2.9/user-doc/html/sum_hills.html
    +  COLVAR file
    +  Trajectory file
    +  Topology file (_if no topology information in trajectory file_)
- COLVAR file must be written to with the same frequency as the trajectory file \
  Regenerable using PLUMED `driver` module: https://www.plumed.org/doc-v2.9/user-doc/html/driver.html

### Formats
- PLUMED-like appearence expected for COLVAR- and FES-files 
- Multiframe capable trajectory, supported by MDAnalysis. \
  Complete list: https://docs.mdanalysis.org/stable/documentation_pages/coordinates/init.html#id2 \
  **Not Supported:** PDBQT, PQR, MMTF, NAMDBIN, IN, FHAIMS, COOR, CRD, GRO
- Some formats might be readable but not writeable, in which case the output is written in XYZ format.

### Installation
The following libraries must be manually installed (e.g. using `pip`):

- shapely
- MDAnalysis

## Running the script
### Input parameters
Use "-h" command to access this explanation from the script.
```
"-traj"    ->  Name of the trajectory file in the MD-output-directory. Format is also used for output-files.  
               REQUIRED FOR RUNNING SCRIPT

"-fes"     ->  Name of the FES histogram bin file in the MD-output-directory.  
               DEFAULT: "fes.dat".

"-colv"    ->  Name of the COLVAR file in the MD-output-directory.  
               DEFAULT: "COLVAR".

"-topo"    ->  MD topology-file name in the MD-output-directory, if trajectory-file does not specify topology.
               DEFAULT: None.

"-thresh"  ->  Specifies threshold for assigning. Input value must correspond with values in the FES-file.
               DEFAULT: Lowest 1/12 of the total energy span.

"-mindist" ->  Smallest allowed distance, at which areas are considered separate minima (in CV units). Must be larger than the diagonal of a single bin.
               DEFAULT: 2% of the FES diagonal.

"-stride"  ->  Reads only every n-th frame of trajectory.
               DEFAULT: 1.

"-md"      ->  MD-output-directory path.
               DEFAULT: Current directory path.

"-png"     ->  Specifies whether a PNG-visualization of the FES should be created. Expects True/False.
               DEFAULT: True.

"-nopbc"   ->  Suppresses the automatic periodicity search (activated when the minima touch the edges). Expects True/False.
               DEFAULT: False.
```

### Created files
- directory "minima" featuring:
   + trajectory files (input format or XYZ) for distinctive minima
   + overview file with average CVs for each minimum
   + FES PNG-file with white outline highlighting minima areas (_if not declared otherwise_)

### Further information
- Periodicity will be detected and considered automatically. (_if not declared otherwise_)
- MDAnalysis frequently fails to open PDB-files created with `CP2K`. A separate PDBReader is implemented to assist in such cases.

## How to cite
If you use FESTA, please cite:

> Istomin, V.; Piccini, GM. FESTA: A Polygon-Based Approach for Extracting Relevant Structures from Free Energy Surfaces Obtained in Molecular Simulations. *J. Chem. Inf. Model.* **2024**, *Article ASAP*. DOI: 10.1021/acs.jcim.4c01022

Link: https://pubs.acs.org/doi/full/10.1021/acs.jcim.4c01022

## Usage advice
It is advised to use FES-histograms with a sufficient bin count (1 million bins resulted in errors below 1%), since the tolerance of the frame extraction and the mapping of the minima depend on the bin quantity.
The script has been tested on Windows and Linux. The script might show reduced performance on Windows workstations.

