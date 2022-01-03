# pSUFER
protein Strain, Unsatisfactoriness, and Frustration findER

Goal: Compute a side-chain suboptimality analysis for PDBs. Suboptimal positions
are defined as those for which Rosetta mutation scanning (filterscan)
identified at least 5 identities with lower energies (<=0) than the native

# Code requirements
All xmls and flags can be found in their respective folders. 
To run the whole code, you are required to install:

1. python > 3.5
2. Rosetta (free for academic users)
3. snakemake > 5.7.4 < 6.8.3
* snakemake [installation instruction](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)  

# Running the code

You can control the pSUFER analysis through pSUFER_config.yaml:
1. rosetta_exec: a path to Rosetta executable
2. trajectories: how many relax trajectories to perform. default is 4
3. energy_cutoff_for_pymol: ΔΔG cutoff below which a mutation is defined as "improving". default is 0, i.e better than WT 
4. identities_at_frustrated_pos: how many identities at a position count for being suboptimal. default 5
5. jump: adjust if jump is needed for rosetta fold tree. default 0  

To start, place your pdb files in ``in/<my dir>/<PDB>.pdb``. **Filenames shouldn't include question marks (?)**. 

Next, run the pSUFER.snakemake code. An example of a commandline would be:

```
snakemake --cores 2000 --snakefile pSUFER.snakefile --configfile pSUFER_config.yaml --latency-wait 600
```
It is highly recommended adding a "--cluster" flag and relevant parameters to run the code on a computer cluster. 

The code will perform the following:
  1. Create 4 refinement trajectories for each pdb in ``in/<my dir>/<PDB>.pdb``. All trajectories and scorefiles would be found in ``temp`` directory, and will be deleted at the end of the run. The lowest energy structure will be placed under ``relax/<my dir>/<PDB>.pdb``.  
  2. Run a mutational scan (filterscan) for each position in the protein. This will result in several new directories:  
  
    a. temp_resfiles: temporary repository for resfiles. Most of it is deleted at the end of the run, but the score logs are saved
    b. resfiles: The final resfiles for each PDB
    c. pymol_pSUFERed: a set of pymol scripts to load the relaxed PDB files and label the suboptimal positions.
   
 
In case you wish to avoid the relaxation step, you can place the structures directly under ``relax/<my dir>/<PDB>.pdb``.

# Results
  All suboptimal positions are summarized in pSUFERED.csv
  
  To visualize the suboptimal positions, two options are available:
  1. open the relaxed structures in pymol (under ``relax/<my dir>/<PDB>.pdb``) and drop the file ``pymol_pSUFERed/<my dir>/<PDB>/load_struct.pml`` on the session.
  2. open the ``pymol_pSUFERed/master.pse`` to see all pdbs in the same pymol session. 


# Running without snakemake
If you wish to avoid using snakemake, pSUFER can be run manually using the xmls and flags provided in the respective folders. 
The basic steps would be:
1. relax the structure several time and take the lowest energy (total_score) variant:
``{rosetta_exec} -parser:protocol xmls/refine.xml -s {input.pdb}  @flags/flags_refinement -parser:script_vars jump={jump_value}`` 

2. run filterscan on all positions in the structure.  
``{rosetta_exec} @flags/filterscan.flags -s {input.pdb} -script_vars current_res={pos} scores_path={path} files_path={path} resfiles_path={path} jump={jump_value}``
The ``unify_res_file`` function in ``scripts/resfily.py`` can be used to order all the individual resfiles created in this step according to the different ΔΔG cutoffs specified in flags/filterscan.flags.
3. check which positions include more than a desired number of mutations (excluding the WT identity) in the resfile created in step 2.
