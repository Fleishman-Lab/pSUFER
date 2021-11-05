# pSUFER
protein Strain, Unsatisfactoriness, and Frustration findER


READ THIS FILE BEFORE RUNNING

Goal: Compute a side-chain suboptimality analysis for PDBs. Suboptimal positions
are defined as those for which Rosetta mutation scanning (filterscan)
identified at least 5 identities with lower energies (<=0) than the native

# Code requirements
All xmls and flags can be found in their respective folders. 
To run the whole code, you are required to install:

1. python > 3.5
2. Rosetta (free for academic users)
3. snakemake > 5.7.4

# Running the code

To start, place your pdb files in in/<my dir>/<PDB>.pdb. **Filenames shouldn't include quiestions marks (?)**. 
Also, adjust path to Rostta executable in pSUFER_config.yaml, under "<rosetta_exec>"
  
Next, run the pSUFFER.snakemake code. An example of a commandline would be:

```
snakemake --cores 2000 --snakefile pSUFER.snakefile --configfile pSUFER_config.yaml --latency-wait 600
```
It is highly recomended adding a "--cluster" flag and relevant parameters to run the code on a computer cluster. 

The code will perform the following:
  1. Create 4 refinement trajectories for each pdb in ``in/<my dir>/<PDB>.pdb``. All trajectories and scorefiles would be found in ``temp`` directory, and will be deleted in the end of the run. The The lowest energy structure will be placed under ``relax/<my dir>/<PDB>.pdb``.  
  2. Run a mutational scan (filterscan) for each position in the protein. This will result in several new directories:  
   a. ``temp_resfiles``: temporary repository for resfiles. Most of it is deleted at the end of the run, but the score logs are saved
   b. ``resfiles``: The final resfiles for each PDB
   c. ``pymol_pSUFERed``: a set of pymol scripts to load the relaxed PDB files and label the suboptimal positions.
   
 
You can control the pSUFER analysis through pSUFER_config.yaml:
1. trajectories: how many relax trajectories to perform.
2. energy_cutoff_for_pymol: What is the energy cutoff by which to decide look for suboptimal positions
3. identities_at_frustrated_pos: how many identities at a position count for being suboptimal.

In case you wish to avoid the relaxation step, you can place the structures directly under ``relax/<my dir>/<PDB>.pdb``.
  
# Results
  All suboptimal positions are summerized in pSUFERED.csv
  To visualized the suboptimal positions, two options are available:
    1. open the relaxed structures in pymol (under ``relax/<my dir>/<PDB>.pdb``) and drop the file ``pymol_pSUFERed/<my dir>/<PDB>/load_struct.pml`` on the session.
    2. open the ``pymol_pSUFERed/master.pse`` to see all pdbs in the same pymol session. 
  
