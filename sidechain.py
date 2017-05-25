#!/usr/bin/env python
docstring='''
sidechain.py [option] backbone.pdb > full.pdb
    use backbone structure "backbone.pdb" to reconstruct full atom model
    "full.pdb". PULCHRA, MODELLER can start from C-alpha trace or 
    main-chain model. Scrwl and RASP must start from main-chain model. 

options:
    -algo={pulchra,scrwl,modeller,foldx,remo,modrefiner}
        Software for full atom reconstruction
        pulchra    - PULCHRA 0.99, 0.999, 3.04, 30.6
        modeller   - MODELLER 9.14, 9.15, 9.16, 9.17
        scrwl      - Scwrl4, RASP1.90
    -execpath=/usr/bin/mod9.16
        Path to executable. Default values:
        If -algo=modeller, searching "mod*" in $PATH .
        If -algo=pulchra, "pulchra" at the current folder or script folder 
        If -algo=scwrl, "Scwrl*" or "RASP" at current folder or script folder
    -super_algo={MMalign,TMscore}
        Algorithm for superposition of models to first template structure
        Default is not performming any superposition
    -super_execpath=/usr/bin/pymol
        Path to executable for superposition. Default is search at current
        folder or script folder
    -seq=seq.fasta
        fasta sequence of output PDB. For -algo={scrwl,modeller}
        each entry is sequence for one chain
'''
import sys,os
import subprocess
import shutil
import random
import re
from string import Template

from pdb2fasta import pdbtxt2seq

# Template for importable modeller script. parameters: $tmp_dir
modeller_template=Template('''#!/usr/bin/env python
import modeller
import modeller.automodel

def comparative_modelling():
    modeller.log.minimal()
    env = modeller.environ()

    # list of path to PDB format template structures
    env.io.atom_files_directory = ["$tmp_dir"]

    a = modeller.automodel.automodel(
        env, 
        alnfile =  "$tmp_dir/alignment.ali", # PIR format alignment file
        knowns =   "template",  # list of template structures
        sequence = "target",    # target sequence
    )

    a.starting_model= 1
    a.ending_model = 1 # number of models to generate
    a.make() # perform the actual homology modeling

if __name__=="__main__":
    comparative_modelling()
''')

# Template for PIR format alignment file. parameter:
# $chain1, $chain2, $target_sequence, $template_sequence
aln_template=Template('''
>P1;template
structure:template:.:$chain1:.:$chain2::::
$template_sequence*

>P1;target
sequence:target:.:$chain1:.:$chain2::::
$target_sequence*
''')

def locate_modeller(PATH=''):
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

def locate_execpath(algo):
    '''locate pulchra/Scwrl*/RASP executable in current folder or folder of 
    script'''
    execpath=""
    if os.path.isfile(algo):
        execpath=os.path.abspath(algo)
    else:
        bindir=os.path.dirname(os.path.abspath(__file__))
        if os.path.isfile(os.path.join(bindir,algo)):
            execpath=os.path.join(bindir,algo)
    if not execpath and algo=="scwrl":
        scwrl_pattern=re.compile("^Scwrl\d+$")
        for f in os.listdir('.'):
            if scwrl_pattern.match(f):
                execpath=os.path.abspath(f)
        if not execpath:
            for f in os.listdir(bindir):
                if scwrl_pattern.match(f):
                    execpath=os.path.abspath(os.path.join(bindir,f))
    if not execpath:
        if algo in ["scwrl","rasp"]:
            execpath=locate_execpath(algo="RASP")
        elif algo=="modeller":
            execpath=locate_modeller()
    if not execpath:
        execpath=algo
    return execpath

def sidechain(CAtxt,algo,execpath,seq,super_algo,super_execpath):
    '''reconstruct side chain'''
    PDBtxt=''
    #### make temporary folder ####
    tmp_dir="/tmp/"+os.getenv("USER")+'/sidechain'+ \
       str(random.randint(1000,9999))+'/'
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    #### parse PDB chains ####
    header_list,sequence_list=pdbtxt2seq(CAtxt,infile="",
        PERMISSIVE="ATOM",outfmt="PDB",allowX=False,SEQRES=False)
    header_list=[chain.lstrip(':') for chain in header_list]
    CAtxt_dict=dict(zip(header_list,['']*len(header_list)))

    #### parse sequence ####
    user_sequence_list=[] # sequence specified by user
    if seq:
        fp=open(seq,'rU')
        txt=fp.read()
        fp.close()
        idx=-1
        for block in txt.split('>'):
            if not block.strip():
                continue
            idx+=1
            sequence=''.join(block.splitlines()[1:])
            if len(sequence)!=len(sequence_list[idx]):
                sys.stderr.write("ERROR! sequence length for chain %s is %d in PDB but %d in fasta"%(
                    header_list[idx],len(sequence_list[idx]),len(sequence)))
                return ''
            user_sequence_list.append(sequence)

    #### split chains ####
    for line in CAtxt.splitlines():
        if line.startswith("END"):
            break
        if line.startswith("ATOM  "):# or line.startswith("HETATM"):
            chain_id=line[21].replace(' ','_')
            CAtxt_dict[chain_id]+=line+'\n'

    #### use pulchra ####
    if algo=="pulchra":
        for chain_id in header_list:
            ## write input PDB ##
            fp=open(tmp_dir+chain_id+".pdb",'w')
            fp.write(CAtxt_dict[chain_id])
            fp.close()
            ## run pulchra ##
            cmd="cd "+tmp_dir+";"+execpath+" -evp "+tmp_dir+chain_id+".pdb"
            subprocess.Popen(cmd,shell=True).communicate()
            ## read pulchra output PDB ##
            if os.path.isfile(tmp_dir+"pul_"+chain_id+".pdb"):
                fp=open(tmp_dir+"pul_"+chain_id+".pdb",'rU')
            elif os.path.isfile(tmp_dir+"rebuilt_"+chain_id+".pdb"):
                fp=open(tmp_dir+"rebuilt_"+chain_id+".pdb",'rU')
            else:
                fp=open(tmp_dir+chain_id+".rebuilt.pdb",'rU')
            for line in fp.read().splitlines():
                if line.startswith("ATOM  "):# or line.startswith("HETATM"):
                    PDBtxt+=line[:21]+chain_id.replace('_',' ')+line[22:]+'\n'
            fp.close()
            PDBtxt+="TER\n"
    elif algo in ["scwrl","rasp"]:
        fp=open(tmp_dir+"backbone.pdb",'w')
        fp.write(CAtxt)
        fp.close()
        cmd="cd "+tmp_dir+";"+execpath+" -i backbone.pdb -o sidechain.pdb"
        if seq:
            fp=open(tmp_dir+"seq.txt",'w')
            fp.write(''.join(user_sequence_list))
            fp.close()
            cmd+=" -s seq.txt"
        subprocess.Popen(cmd,shell=True).communicate()
        fp=open(tmp_dir+"sidechain.pdb",'rU')
        PDBtxt=fp.read()
        fp.close()
    elif algo=="modeller": # modeller multimeric modelling
        ## prepare template file ##
        fp=open(tmp_dir+"template.pdb",'w')
        for chain_id in header_list:
            fp.write(CAtxt_dict[chain_id]+"TER\n")
        fp.write("END\n")
        fp.close()
        
        ## prepare alignment file ##
        fp=open(tmp_dir+"alignment.ali",'w')
        if seq:
            fp.write(aln_template.substitute(
                chain1=header_list[0].replace('_',' '),
                chain2=header_list[-1].replace('_',' '),
                template_sequence='/'.join(user_sequence_list),
                target_sequence='/'.join(sequence_list)))
        else:
            fp.write(aln_template.substitute(
                chain1=header_list[0].replace('_',' '),
                chain2=header_list[-1].replace('_',' '),
                template_sequence='/'.join(sequence_list),
                target_sequence='/'.join(sequence_list)))
        fp.close()

        ## write modeller script ##
        fp=open(tmp_dir+"model-mulitchain.py",'w')
        fp.write(modeller_template.substitute(tmp_dir=tmp_dir))
        fp.close()

        ## run modeller ##
        cmd="cd "+tmp_dir+";"+execpath+" "+tmp_dir+"model-mulitchain.py"
        subprocess.Popen(cmd,shell=True).communicate()

        ## parse result ##
        target_pattern=re.compile("^target\.[\w]*\.pdb$")
        for f in os.listdir(tmp_dir):
            if target_pattern.match(f):
                model_pdb=tmp_dir+f
                if super_algo:
                    pwd=os.getcwd()
                    os.chdir(tmp_dir)
                    model_fit_pdb=superpose(model_pdb=f,
                        native_pdb="template.pdb",
                        algo=super_algo, execpath=super_execpath)
                    os.chdir(pwd)
                    model_pdb=tmp_dir+model_fit_pdb
                fp=open(model_pdb,'rU')
                PDBtxt=fp.read()
                fp.close()

    #### cleanup temporary folder ####
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    return PDBtxt

if __name__=="__main__":
    #### Default parameters for options ####
    algo='pulchra'    # pulchra, modeller, scwrl, rasp 
    execpath=''   
    super_algo=''
    super_execpath=''
    seq=''

    #### parse arguments ####
    argv=[] # input FASTA format alignment files
    for arg in sys.argv[1:]:
        if not arg.startswith('-'):
            argv.append(arg)
        elif arg.startswith('-algo='):
            algo=arg[len("-algo="):].lower()
        elif arg.startswith('-execpath='):
            execpath=arg[len("-execpath="):]
        elif arg.startswith('-super_algo='):
            super_algo=arg[len("-super_algo="):]
        elif arg.startswith('-super_execpath='):
            super_execpath=arg[len("-super_execpath="):]
        elif arg.startswith('-seq='):
            seq=arg[len("-seq="):]
        else:
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()

    if not len(argv):
        sys.stderr.write(docstring)
        exit()

    #### check if superpose module is importable ####
    importable_superpose_module=False
    if super_algo:
        try:
            from superpose import superpose
            importable_superpose_module=True
        except Exception,e:
            sys.stderr.write(str(e)+"\nERROR! Cannot import superpose module\n")
            exit()
        if not super_execpath:
            super_execpath=locate_execpath(super_algo.split('-')[0])

    #### locate executables ####
    if not execpath:
        execpath=locate_execpath(algo)

    #### perform reconstruction ####
    fp=open(argv[0],'rU')
    CAtxt=fp.read()
    fp.close()
    PDBtxt=sidechain(CAtxt,algo,execpath,seq,
        super_algo,super_execpath)
    if len(argv)>1:
        fp=open(argv[-1],'w')
        fp.write(PDBtxt)
        fp.close()
    else:
        sys.stdout.write(PDBtxt)
