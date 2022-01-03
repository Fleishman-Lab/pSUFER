import os
import glob
import json
import shutil
from scripts.resfile import unify_res_file

#To run on cluster:
#snakemake --cluster "bsub -u /dev/null -q new-short -R rusage[mem=1024] -C 1 -G fleishman-wx-grp-lsf -o logs/{rule}.{wildcards}.out -e logs/{rule}.{wildcards}.err" --cores 2000 --snakefile pSUFER.snakefile --configfile pSUFER_config.yaml --latency-wait 600

DIRS,PDBS,=glob_wildcards('in/{DIR}/{PDB}.pdb')
targets=[]
for i in range(len(DIRS)):
    targets.append(f'pymol_pSUFERed/{DIRS[i]}/{PDBS[i]}/load_struct.pml')
#print(targets)

localrules: pymol_scripts
rule all:
  input:
    targets

def count_residues(file):
    res_num=0
    pdb=open(file,'r')
    for line in pdb.readlines():
      fields=line.split()
      if(len(fields)>2 and fields[0]=="ATOM" and fields[2]=="CA"):
        res_num=res_num+1
    return res_num

def resfile_names(wildcards):
    resnums=count_residues(f'in/{wildcards.D}/{wildcards.P}.pdb')
    list=[]
    for res in range(1, resnums+1):
        list.append(f'temp_resfiles/{wildcards.D}/{wildcards.P}/res_{res}')
    return(list)

def expand_resfile_names_with_thresholds(names,dir):
    expanded_names=[]
    for n in names:
        for t in config['energy_cutoffs_rosetta']:
            file_name=f'{dir}/{n}.{t}'
            if(os.path.exists(file_name)):
                expanded_names.append(file_name)
    return(expanded_names)

#read a resfile and return only the lines that indicate suboptimality
def frustrated_lines(dir, pdb, config):
    with open(resfile_name(dir, pdb, config),'r') as resfile_handle:
        frustrated=list(filter(lambda l: len(l.split())>2 and len(l.split()[-1])>=config['identities_at_frustrated_pos'],
        resfile_handle.readlines()))
    return(frustrated)

def resfile_name(dir,pdb, config):
    resfile=f'resfiles/{dir}/{pdb}/designable_aa_resfile.{config["energy_cutoff_for_pymol"]}'
    return(resfile)

rule rosetta_relax:
  input:
    pdb='in/{D}/{A}.pdb'
  output:
    pdb=temporary('temp/{D}/{num}?{A}_0001.pdb')
  params:
    dir='temp/{D}'
  run:
    command = f'{config["rosetta_exec"]} -parser:protocol xmls/refine.xml -s {input.pdb}  @flags/flags_refinement -out:prefix {params.dir}/{wildcards.num}? -out:file:scorefile {params.dir}/{wildcards.A}.sc -parser:script_vars jump={config["jump"]}'
    print(command)
    os.system(command)

rule select_best:
  input:
    pdb=expand('temp/{{D}}/{num}?{{A}}_0001.pdb', num=range(config["trajectories"]))
  output:
    pdb='relax/{D}/{A}.pdb'
  params:
    score='temp/{D}/{A}.sc'
  run:
    command = f'cp `python scripts/parse_scores.py --scorefile {params.score} --operations "sort total_score Asc cols description" | ' + "awk 'NR==2{print $0 \".pdb\"}'`"+" {out_pdb}".format(out_pdb=output.pdb)
    print(command)
    os.system(command)


rule filter_scan:
  input:
    pdb='relax/{D}/{P}.pdb'
  output:
    resfile=touch('temp_resfiles/{D}/{P}/res_{pos}') #of the form: temp_resfiles/D/P/res_{pos}
  run:
    path=os.path.dirname(output.resfile)
    pos=output.resfile.split('_')[-1]
    command = f'{config["rosetta_exec"]} @flags/filterscan.flags -s {input.pdb} -script_vars current_res={pos} scores_path={path} files_path={path} resfiles_path={path} jump={config["jump"]}'
    print(command)
    os.system(command)

rule unify_resfiles:
    input:
        resfile=resfile_names #of the form: temp_resfiles/D/P/res_{pos}
    output:
        touch('resfiles/{D}/{P}/designable_aa_resfile.done')
    run:
        curr_dir=os.getcwd()
        resfiles_thresholds=expand_resfile_names_with_thresholds(input.resfile,curr_dir)
        unify_res_file(resfiles_thresholds, pdb=f'{curr_dir}/relax/{wildcards.D}/{wildcards.P}.pdb', out_path=f'resfiles/{wildcards.D}/{wildcards.P}/')

rule pymol_scripts:
    input:
        resfiles='resfiles/{D}/{P}/designable_aa_resfile.done'
    output:
        pymol_script='pymol_pSUFERed/{D}/{P}/load_struct.pml'
    run:
        r=input.resfiles
        target_resfile=r.split('.')[-2] + f'.{config["energy_cutoff_for_pymol"]}'
        pymol_script_file=open(output.pymol_script,'a')
        pymol_script_file.write(f'load relax/{wildcards.D}/{wildcards.P}.pdb\nshow lines\nremove hydrogens\n')
        pymol_script_file.write('color cyan, chain A and {p}\ncolor gray, chain B and {p}\n'.format(p=wildcards.P)) #in case of two chains
        frustrated=frustrated_lines(wildcards.D, wildcards.P, config)
        selection=list(map(lambda l: f'{wildcards.P} AND chain {l.split()[1]} AND resi {l.split()[0]}', frustrated))
        selection_string=' OR '.join(selection)
        selection_name=f'{wildcards.P}_frustrated'
        pymol_script_file.write(f'select {selection_name}, {selection_string}\n')
        pymol_script_file.write(f'show sticks, {selection_name}\n')
        pymol_script_file.write(f'color yellow, {selection_name}\n')
        pymol_script_file.write('color atomic, (not elem C)\n')
        pymol_script_file.write('hide sticks, name c+n+o\n')
        pymol_script_file.write('hide lines, name c+n+o\n')
        pymol_script_file.close()

        frustrated_pos=','.join(list(map(lambda f: f'{f.split()[0]}{f.split()[1]}', frustrated)))
        with open('pSUFERed.csv','a') as frustrated_table:
            frustrated_table.write(f'{wildcards.P},{frustrated_pos}\n')

        pymol_master_script=open('pymol_pSUFERed/master.pml','a')
        pymol_master_script.write(f'load {output.pymol_script}\n')
        pymol_master_script.close()
onsuccess:
    print('writing pymol sessions into pymol_pSUFERed/')
    pymol_master_script=open('pymol_pSUFERed/master.pml','a')
    pymol_master_script.write('order *, yes\nsave pymol_pSUFERed/master.pse\n')
    pymol_master_script.close()
    os.system('pymol -cq -r pymol_pSUFERed/master.pml')
    print('done writing pymol sessions. Open pymol_pSUFERed/master.pse')
    print('deleting unnecessary files in temp_resfiles and logs.\n')
    os.system('nohup find temp_resfiles/ -type f -not -name "*.log" -exec rm -f {} \;')
    print('tarzipping this directory\n')
    os.system('tar -cvf temp_resfiles.tar temp_resfiles')
    os.system('gzip -f temp_resfiles.tar')
    os.system('tar -cvf logs.tar logs')
    os.system('gzip -f logs.tar')
    print('deleting the directories\n')
    shutil.rmtree('temp_resfiles', ignore_errors=True)
    shutil.rmtree('logs', ignore_errors=True)
