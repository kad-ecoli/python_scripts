#!/usr/bin/env python
docstring='''
HomologyModelling.py [option] seq.fasta
    Comparative Modelling by Satisfaction of Spatial Restraints using MODELLER 
    guided by FASTA format gapped multiple sequence alignment "seq.fasta".

    >target.pdb
    AAAAAAAA
    >template1.pdb
    BBBBBBBB
    >template2.pdb
    CCCCCCCC

    "template1.pdb", "template2.pdb" ... are template structures.
    "target.pdb" is the target structure. 
    Prediction target must be the first entry.

options:
    -execpath=/usr/bin/mod9.15
        Path to MODELLER executable. Default is searching "mod*" in $PATH
    -md={none,low,mid,high}
        molecular dynamics refinement level. Default is none.
    -hetatm={False,True}
        whether to include ligands
    -max_models=1
        number of models to generate
    -super_algo={TMalign,TMscore,MMalign,
                 pymol-super,pymol-cealign,pymol-align}
        Algorithm for superposition of models to first template structure
        Default is not performming any superposition
    -super_execpath=/usr/bin/pymol
        Path to executable for superposition. Default is searching in $PATH
    -psipred=seq.ss2
        PSIPRED format secondary structure prediction file
    -disulfide=disulfide.txt
        list of cystein pairs forming disulfide bond
'''
import sys,os
import subprocess
import re
import shutil
from string import Template
import random
import textwrap

# Template for importable modeller script. parameters: 
# $atom_files_directory, $hetatm, $alnfile, $knowns, $sequence, 
# $max_models, $optimization, $ss_restraint, $disulfide_restraint
modeller_template=Template('''#!/usr/bin/env python
import modeller
import modeller.automodel

def comparative_modelling():
    modeller.log.minimal()
    env = modeller.environ()
    env.io.hetatm = $hetatm     # whether to include ligands (True/False)

    # list of path to PDB format template structures
    env.io.atom_files_directory = $atom_files_directory

    class MyModel(modeller.automodel.automodel):
        def special_restraints(self, aln):
            rsr = self.restraints
$ss_restraint
$disulfide_restraint

    a = MyModel(
        env, 
        alnfile =  "$alnfile",  # PIR format alignment file
        knowns =    $knowns,    # list of template structures
        sequence = "$sequence", # target sequence
    )

    $optimization  # molecular dynamics refinement

    a.starting_model= 1
    a.ending_model = $max_models # number of models to generate
    a.make() # perform the actual homology modeling

if __name__=="__main__":
    comparative_modelling()
''')

modeller_optimization_dict={ # molecular dynamics refinement code
"none":'',
"low" :'''
    a.library_schedule=modeller.automodel.autosched.very_fast # low VTFM
    a.md_level=modeller.automodel.refine.very_fast            # low MD
           ''',
"mid" :'''
    a.library_schedule=modeller.automodel.autosched.fast # thorough VTFM
    a.max_var_iterations=300
    a.md_level=modeller.automodel.refine.fast            # thorough MD
    a.repeat_optimization=2                # Repeat the whole cycle twice
           ''',
"high":'''
    a.library_schedule=modeller.automodel.autosched.slow # very thorough VTFM
    a.max_var_iterations=300
    a.md_level=modeller.automodel.refine.slow            # very thorough MD
    a.repeat_optimization=2                # Repeat the whole cycle twice and
    a.max_molpdf=1e6                       # do not stop unless obj.func.>1e6
           ''',
}

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

def FixChainID(PDBfile="pdb.pdb"):
    '''Convert all chain ID in "PDBfile" into chain A'''
    fp=open(PDBfile,'rU')
    pdb_lines=fp.read().splitlines()
    fp.close()
    
    pdb_txt=''
    for line in pdb_lines:
        if line.startswith("HETATM") or line.startswith("ATOM  "):
            pdb_txt+=line[:21]+'A'+line[22:80]+'\n'
    fp=open(PDBfile,'w')
    fp.write(pdb_txt)
    fp.close()
    return PDBfile

def parse_psipred(psipred=''):
    '''parse psipred secondary restraints into modeller command'''
    if not psipred:
        return ''
    ss_restraint='' # modeller command to add psipred restraints
    fp=open(psipred,'rU')
    psipred_txt=fp.read()
    fp.close()
    ss_pattern=re.compile("\s*(\d+)\s+[A-Za-z]\s+([HEC])")

    is_helix=False # currently reading helix
    is_strand=False # currently reading strand
    for line in ss_pattern.findall(psipred_txt):
        resi,secstruct=line # residue index and secondary structure
        resi=int(resi)
        if secstruct=='C':
            if is_helix or is_strand:
                ss_restraint+="'%d:')))\n"%(resi-1) # strand or helix ended
            is_helix=False
            is_strand=False
        elif secstruct=='H':
            if is_strand:
                ss_restraint+="'%d:')))\n"%(resi-1) # strand ended
                is_strand=False
            if is_helix==False: # start residue of helix
                ss_restraint+=' '*12+"rsr.add(secondary_structure"+ \
                    ".alpha(self.residue_range('%d:',"%resi
                is_helix=True
        elif secstruct=='E':
            if is_helix:
                ss_restraint+="'%d:')))\n"%(resi-1) # strand ended
                is_helix=False
            if is_strand==False: # start residue of helix
                ss_restraint+=' '*12+"rsr.add(secondary_structure"+ \
                    ".strand(self.residue_range('%d:',"%resi
                is_strand=True
    if is_helix or is_strand:
        ss_restraint+="'%d:')))\n"%resi # one helix ended
    return ss_restraint

def parse_disulfide(disulfide=''):
    '''parse list of cystein pairs forming disulfide bond'''
    if not disulfide:
        return ''
    fp=open(disulfide,'rU')
    disulfide_txt=fp.read()
    fp.close()
    disulfide_restraint=' '*8+"def special_patches(self,alnfile):\n"

    for line in disulfide_txt.splitlines():
        if not line.strip() or line.strip().startswith('#'):
            continue
        cys1,cys2=line.split()
        disulfide_restraint+=' '*12+"self.patch(residue_type='DISU',"+ \
        "residues=(self.residues['%s'],self.residues['%s']))\n"%(cys1,cys2)
    return disulfide_restraint

def HomologyModelling(alignment="alignment.fasta", execpath="mod9.15",
    hetatm=False, md="none",max_models=1,psipred='', disulfide=''):
    '''Homolgy Modelling by Modeller.
    return a tuple for filename of target and first template structure.
    filename of target will be empty if MODELLER models cannot be generated.

    alignment  - FASTA format input alignment
    md         - molecular dynamics refinement ("none","low", "mid", "high")
    max_models - number of models to generate
    hetatm     - whether to include ligands (True, False)'''

    # create temporary folder
    tmp_dir="/tmp/"+os.getenv("USER")+'/modeller'+md+ \
       str(random.randint(1000,9999))+'/'
    while (os.path.isdir(tmp_dir)):
        tmp_dir="/tmp/"+os.getenv("USER")+'/modeller'+md+ \
            str(random.randint(1000,9999))+'/'
    os.makedirs(tmp_dir)
    
    # read FASTA format multiple sequence alignment
    fp=open(alignment,'rU')
    fasta_txt=fp.read()
    fp.close()

    # parse input multiple sequence alignment
    target_filename=''
    target_sequence=''
    template_filename_lst=[]
    template_sequence_lst=[]
    for entry in fasta_txt.split('>'):
        if not entry.strip():
            continue
        entry_lst=entry.splitlines()
        filename=entry_lst[0]
        sequence=textwrap.fill(''.join(entry_lst[1:]).rstrip('*')+'*',75)
        if not target_filename:
            target_filename=filename
            target_sequence=sequence
        else:
            template_filename_lst.append(filename)
            template_sequence_lst.append(sequence)

    # write PIR format multiple sequence alignment
    alnfile=tmp_dir+"alignment.ali"
    fp=open(alnfile,'w')
    ali_txt=''
    for idx,template_filename in enumerate(template_filename_lst):
        template_filename_tmp=tmp_dir+str(idx)+".pdb"
        shutil.copy(template_filename,template_filename_tmp)
        FixChainID(template_filename_tmp)
        ali_txt+=">P1;%u\nstructure:%u:.:A:.:A::::\n%s\n\n"%(
            idx,idx,template_sequence_lst[idx])
    ali_txt+=">P1;target\nsequence:target:.:.:.:.::::\n%s\n\n"%(
            target_sequence)
    fp.write(ali_txt)
    fp.close()
    knowns=map(str,range(len(template_filename_lst)))

    # write homology modelling script
    modeller_script=tmp_dir+"modeller_script.py"
    fp=open(modeller_script,'w')
    modeller_parameter_dict=dict(
        hetatm=hetatm,                 # include ligand
        atom_files_directory=[tmp_dir],# PDB template folder
        alnfile=alnfile,               # alignment file
        knowns=knowns,                 # list of template structures
        sequence="target",             # target sequence
        max_models=max_models,         # number of models to generate
        optimization=modeller_optimization_dict[md], # refinement
        ss_restraint=parse_psipred(psipred), # psipred prediction
        disulfide_restraint=parse_disulfide(disulfide),
    )
    fp.write(modeller_template.substitute(modeller_parameter_dict))
    fp.close()

    # perform homology modelling script
    command="cd "+tmp_dir+"\n"+execpath+' '+modeller_script
    stdout,stderr= subprocess.Popen(command, shell=True).communicate()
    
    # grep output model files
    target_pattern=re.compile("^target\.[\w]*\.pdb$")
    target_model_lst=[f for f in os.listdir(tmp_dir) if target_pattern.match(f)]
    if not target_model_lst:
        sys.stderr.write("ERROR! Cannot generate target structures\n")
        return '',template_filename_lst[0]

    # parse output models
    target_txt_lst=[]
    for idx,target_model in enumerate(target_model_lst):
        model_idx=str(idx+1)
        target_txt_lst.append("MODEL"+' '*(9-len(model_idx))+model_idx)
        fp=open(tmp_dir+target_model,'rU')
        target_txt_lst+=[line for line in fp.read().splitlines()
            if not line.strip()=="END"]
        fp.close()
        target_txt_lst.append("ENDMDL")
    target_txt='\n'.join(target_txt_lst)+'\nEND\n'

    # write multimodel target PDB
    target_filename_tmp=tmp_dir+"target.pdb"
    fp=open(target_filename_tmp,'w')
    fp.write(target_txt)
    fp.close()
    shutil.copy(target_filename_tmp,target_filename)

    # clean up temporary folder
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    return target_filename,template_filename_lst[0]

if __name__=="__main__":
    # Default parameters for options
    execpath=''       # search MODELLER by locate_modeller()
    md="none"         # no molecular dynamics refinement
    hetatm=False      # do not include ligand
    max_models=1      # generate one model
    super_algo=''     # do not superpose model to template
    super_execpath='' # search superposition executable in $PATH
    psipred=''        # do not use secondary structure restraint
    disulfide=''      # do not define disulfide bond

    # parse arguments
    argv=[] # input FASTA format alignment files
    for arg in sys.argv[1:]:
        if not arg.startswith('-'):
            argv.append(arg)
            continue
        if arg.startswith('-execpath='):
            execpath=arg[len("-execpath="):]
        elif arg.startswith('-md='):
            md=arg[len("-md="):]
        elif arg.startswith('-hetatm='):
            hetatm=arg[len("-hetatm="):]
        elif arg.startswith('-max_models='):
            max_models=arg[len("-max_models="):]
        elif arg.startswith('-super_algo='):
            super_algo=arg[len("-super_algo="):]
        elif arg.startswith('-super_execpath='):
            super_execpath=arg[len("-super_execpath="):]
        elif arg.startswith('-psipred='):
            psipred=arg[len("-psipred="):]
        elif arg.startswith('-disulfide='):
            disulfide=arg[len("-disulfide="):]
        else:
            print >>sys.stderr, "ERROR! Unknown argument "+arg
            exit()

    if not len(argv):
        sys.stderr.write(docstring)
        exit()
    
    # locate modeller
    if not execpath:
        execpath=locate_modeller()

    # check if superpose module is importable
    importable_superpose_module=False
    if super_algo:
        try:
            from superpose import superpose
            importable_superpose_module=True
        except Exception,e:
            sys.stderr.write(str(e)+"\nERROR! Cannot import superpose module\n")
            exit()

    # run MODELLER
    for alignment in argv:
        model_pdb,native_pdb= HomologyModelling(alignment=alignment, 
            execpath=execpath, hetatm=hetatm, md=md, max_models=max_models,
            psipred=psipred, disulfide=disulfide)
        
        if not importable_superpose_module:
            continue
        # superpose model to first template
        model_fit_pdb=superpose(model_pdb=model_pdb,native_pdb=native_pdb,
                algo=super_algo, execpath=super_execpath)
        shutil.move(model_fit_pdb,model_pdb)
