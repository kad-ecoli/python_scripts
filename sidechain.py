#!/usr/bin/env python
docstring='''
sidechain.py [option] backbone.pdb > full.pdb
    use backbone structure "backbone.pdb" to reconstruct full atom model
    "full.pdb". PULCHRA, MODELLER can start from C-alpha trace or 
    main-chain model. Scrwl and RASP must start from main-chain model. 

options:
    -algo={pulchra,modeller,scrwl,remo,modrefiner}
        Software Program to use for full atom reconstruction
        pulchra    - PULCHRA 0.99, 0.999, 3.04, 30.6
        modeller   - MODELLER 9.14, 9.15, 9.16, 9.17
        scrwl      - Scwrl4, RASP1.90
    -execpath=/usr/bin/mod9.16
        Path to executable. Default values:
        If -algo=modeller, searching "mod*" in $PATH .
        If -algo=pulchra, "pulchra" at the directory of this script or in 
            $PATH.
        If -algo=scwrl, "Scwrl*" or "RASP" at the directory of this script
            or in $PATH.
    -super_algo={TMalign,TMscore,MMalign,
            pymol-super,pymol-cealign,pymol-align}
        Algorithm for superposition of models to first template structure
        Default is not performming any superposition
    -super_execpath=/usr/bin/pymol
        Path to executable for superposition. Default is searching in $PATH
'''
import sys,os
import subprocess
import shutil
try: # check if HomologyModelling module is importable
    from HomologyModelling import HomologyModelling,locate_modeller
    importable_HomologyModelling_module=True
except ImportError:
    importable_HomologyModelling_module=False

def locate_execpath(algo,PATH=''):
    '''locate modeller executable in "PATH", which defaults to environmental 
    variable $PATH '''
    execpath=""
    if not PATH:
        PATH=os.getenv("PATH")

    mod_pattern=re.compile("mod\d{1,2}[.v]\d{1,2}")
    for d in PATH.split(os.path.pathsep):
        file_lst=[e for e in sorted(os.listdir(d)) if \
            os.path.isfile(os.path.join(d,e)) and mod_pattern.match(e)]
        if file_lst:
            execpath=os.path.join(d,file_lst[-1])
            break
    return execpath

if __name__=="__main__":
    # Default parameters for options
    algo='pulchra'    # pulchra, modeller, scwrl, rasp 
    execpath=''   
    super_algo=''
    super_execpath=''

    # parse arguments
    argv=[] # input FASTA format alignment files
    for arg in sys.argv[1:]:
        if not arg.startswith('-'):
            argv.append(arg)
        elif arg.startswith('-algo='):
            algo=arg[len("-algo="):].lower()
        elif arg.startswith('-execpath='):
            execpath=arg[len("-execpath="):]
        elif arg.startswith('-super_algo='):
            super_algo=arg[len("-super_algo="):].lower()
        elif arg.startswith('-super_execpath='):
            super_execpath=arg[len("-super_execpath="):]
        else:
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()

    if not len(argv):
        sys.stderr.write(docstring)
        exit()

    # check if superpose module is importable
    importable_superpose_module=False
    if super_algo:
        try:
            from superpose import superpose
            importable_superpose_module=True
        except Exception,e:
            sys.stderr.write(e+"\nERROR! Cannot import superpose module\n")
            exit()

    # check if modeller can be used
    if algo=="modeller":
        if not importable_HomologyModelling_module:
            sys.stderr.write(e+"\nERROR! Cannot import HomologyModelling module\n")
            exit()
        elif not execpath:
            execpath=locate_modeller()
