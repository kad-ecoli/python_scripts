#!/usr/bin/env python
# 2016-01-05 Chengxin Zhang
docstring='''
superpose.py [options] model.pdb native.pdb
    Superpose "model.pdb" to "native.pdb". Only the first MODEL of 
    "native.pdb" will be used. All MODEL of "model.pdb" will be superposed.
    Write a new file at current directory called "model_fit.pdb", which 
    contains coordinates of "model.pdb" after superposition

options:
    -algo={TMalign,TMscore,MMalign,pymol-super,pymol-cealign,pymol-align,matrix}
        Algorithm used for superposition. 
        "matrix" means input TMalign style matrix. In this case, "native.pdb"
        is the matrix file.
    -execpath=/usr/bin/TMalign
        Path to executable. default is searching in the same folder of this
        script and then search in $PATH
    -writePDB={true,false}
        whether write PDB after superposition
    -offset=0
        add "offset" to residue index in model.pdb
'''
import sys,os
import subprocess
import re
import shutil
from string import Template
import random

# pattern for LOMETS templates
initdat_pattern=re.compile("\n\s*\d+\s+[-]{0,1}[.\d]+\s+\d+\s+\w+[\w\W]+?\n")

def sup_pymol(model_pdb,native_pdb, algo="pymol-super",
    execpath="/usr/bin/pymol",offset=0,check_nmr=True):
    '''Superpose "model_pdb" to "native_pdb" using pymol

    model_pdb  - moved pdb structure
    native_pdb - fixed pdb structure
    algo       - algorithm for supersition {"TMalign","TMscore","MMalign"}
    execpath   - path to TMalign/TMscore/MMalign executable
    offset     - obsolete parameter.
    check_nmr  - check whether "model_pdb" is multi-model NMR structure
    '''
    tmp_dir="/tmp/"+os.getenv("USER")+'/'+algo+ \
       str(random.randint(100000000,999999999))+'/'
    while(os.path.isdir(tmp_dir)):
        tmp_dir="/tmp/"+os.getenv("USER")+'/'+algo+ \
            str(random.randint(100000000,999999999))+'/'
    os.makedirs(tmp_dir)
    shutil.copy(model_pdb,tmp_dir+"model.pdb")
    shutil.copy(native_pdb,tmp_dir+"native.pdb")

    # check if native_pdb is multi-model NMR structure
    fp=open(tmp_dir+"native.pdb",'rU')
    txt=fp.read()
    fp.close()
    if "ENDMDL" in txt:
        fp=open(tmp_dir+"native.pdb",'w')
        fp.write(txt.split("ENDMDL")[0]+"ENDMDL\n")
        fp.close()

    # Check if model_pdb is multi-model NMR structure
    fp=open(tmp_dir+"model.pdb",'rU')
    txt=fp.read()
    fp.close()
    if check_nmr and ("\nMODEL " in txt or initdat_pattern.search(txt)):
        if "\nMODEL " in txt:
            MODEL_start_lst=[e.start()+1 for e in re.finditer('\nMODEL ',txt)]
            if txt.startswith("MODEL"):
                MODEL_start_lst=[0]+MODEL_start_lst
        else:
            MODEL_start_lst=[e.start()+1 for e in initdat_pattern.finditer(txt)]

        super_txt="cd $tmp_dir\nload $native,fix,1\n"
        ce_txt   ="cd $tmp_dir\nload $native,fix,1\n"
        for idx,start in enumerate(MODEL_start_lst):
            if idx==len(MODEL_start_lst)-1:
                MODEL=txt[start:]
            else:
                MODEL=txt[start:MODEL_start_lst[idx+1]]
                
            if "\nATOM  " in MODEL:
                single_model_file=tmp_dir+str(idx)
                fp=open(single_model_file,'w')
                fp.write(MODEL)
                fp.close()
                super_txt+='''
load %s,%u
$pymol_algo %u,fix
save %s,%u
'''%(single_model_file,idx,idx,single_model_file+"_fit",idx)
                ce_txt+='''
load %s,%u
$pymol_algo fix,%u
save %s,%u
'''%(single_model_file,idx,idx,single_model_file+"_fit",idx)

        pymol_template=Template(
            (ce_txt if algo=="pymol-cealign" else super_txt) \
#+'''join_states move, ('''+ \
#'|'.join(map(str,range(len(MODEL_start_lst)))) + '''), mode=0
#save model_fit.pdb,move,0
#'''
        )
        pymol_txt=pymol_template.substitute(
            tmp_dir   =tmp_dir,
            native    =tmp_dir+"native.pdb",
            pymol_algo=algo[len("pymol-"):],
        )
        
        pymol_pml=tmp_dir+"pymol.pml"
        
        fp=open(pymol_pml,'w')
        fp.write(pymol_txt)
        fp.close()

        command=execpath+" -c "+pymol_pml
        stdout,stderr= subprocess.Popen(command, shell=True,
            stdout = subprocess.PIPE).communicate()

        RMSD_pattern=re.compile("RMSD{0,1}\s*={0,1}\s*([.\d]+)")
        RMSD_lst=RMSD_pattern.findall(stdout)
        sys.stdout.write(''.join(["RMSD="+e+";"+'\n' for e in RMSD_lst]))
        
        pdb_txt=''

        for idx,start in enumerate(MODEL_start_lst):
            if idx==len(MODEL_start_lst)-1:
                MODEL=txt[start:]
            else:
                MODEL=txt[start:MODEL_start_lst[idx+1]]
                
            if not "\nATOM  " in MODEL:
                pdb_txt+=MODEL
            else:
                single_model_file=tmp_dir+str(idx)+"_fit"
                atom_idx=MODEL.find("\nATOM  ")
                MODEL_lst_head=[e for e in MODEL[:atom_idx].splitlines() \
                    if e.strip()]
                MODEL_lst_tail=[e for e in MODEL[atom_idx:].splitlines() \
                    if e.strip() and not e.startswith("ATOM  ") \
                                 and not e.startswith("HETATM")]
                fp=open(single_model_file,'rU')
                MODEL_lst_body=[e for e in fp.read().splitlines() if
                    e.startswith("ATOM  ") or e.startswith("HETATM")]
                fp.close()
                pdb_txt+='\n'.join(
                    MODEL_lst_head+MODEL_lst_body+MODEL_lst_tail)+'\n'

        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)

        return (pdb_txt,list(map(float,RMSD_lst)),[])
    
    super_template=Template('''
cd $tmp_dir
load $native,fix,1
load $model,move
$pymol_algo move,fix
save model_fit.pdb,move
quit
''')

    ce_template=Template('''
cd $tmp_dir
load $native,fix,1
load $model,move
$pymol_algo fix,move
save model_fit.pdb,move
quit
''')
    
    if algo=="pymol-cealign":
        pymol_template=ce_template
    else:
        pymol_template=super_template

    pymol_txt=pymol_template.substitute(
    tmp_dir   =tmp_dir,
    model     =tmp_dir+"model.pdb",
    native    =tmp_dir+"native.pdb",
    pymol_algo=algo[len("pymol-"):],
    )

    pymol_pml=tmp_dir+"pymol.pml"
    
    fp=open(pymol_pml,'w')
    fp.write(pymol_txt)
    fp.close()

    command=execpath+" -c "+pymol_pml
    stdout,stderr= subprocess.Popen(command, shell=True,
        stdout = subprocess.PIPE).communicate()
    RMSD_pattern=re.compile("RMSD{0,1}\s*={0,1}\s*([.\d]+)")
    RMSD=RMSD_pattern.findall(stdout)[-1]
    sys.stdout.write("RMSD="+RMSD+";")

    #fp=open(tmp_dir+"model_fit.pdb",'rU')
    #pdb_txt=fp.read()
    #fp.close()

    fp=open(tmp_dir+"model.pdb",'rU')
    MODEL=fp.read()
    fp.close()

    if not "ATOM  " in MODEL:
        pdb_txt=MODEL
    else:
        atom_idx=MODEL.find("ATOM  ")
        MODEL_lst_head=[e for e in MODEL[:atom_idx].splitlines() \
            if e.strip()]
        MODEL_lst_tail=[e for e in MODEL[atom_idx:].splitlines() \
            if e.strip() and not e.startswith("ATOM  ") \
                         and not e.startswith("HETATM")]
        fp=open(tmp_dir+"model_fit.pdb",'rU')
        MODEL_lst_body=[e for e in fp.read().splitlines() if \
            e.startswith("ATOM  ") or e.startswith("HETATM")]
        fp.close()

        pdb_txt='\n'.join(MODEL_lst_head+MODEL_lst_body+MODEL_lst_tail)+'\n'


    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    return (pdb_txt,[float(RMSD)],[])

def sup_TMalign(model_pdb,native_pdb, algo="TMalign",
    execpath="/usr/bin/TMalign",offset=0,check_nmr=True):
    '''Superpose "model_pdb" to "native_pdb" using TMalign/TMscore/MMalign

    model_pdb  - moved pdb structure
    native_pdb - fixed pdb structure
    algo       - algorithm for supersition {"TMalign","TMscore","MMalign"}
    execpath   - path to TMalign/TMscore/MMalign executable
    offset     - add "offset" to residue index of model.
    check_nmr  - check whether "model_pdb" is multi-model NMR structure
    '''
    tmp_dir="/tmp/"+os.getenv("USER")+'/'+algo+ \
       str(random.randint(100000000,999999999))+'/'
    while(os.path.isdir(tmp_dir)):
        tmp_dir="/tmp/"+os.getenv("USER")+'/'+algo+ \
            str(random.randint(100000000,999999999))+'/'
    os.makedirs(tmp_dir)
    shutil.copy(model_pdb,tmp_dir+"model.pdb")
    shutil.copy(native_pdb,tmp_dir+"native.pdb")

    # check if native_pdb is multi-model NMR structure
    fp=open(tmp_dir+"native.pdb",'rU')
    txt=fp.read()
    fp.close()
    if "ENDMDL" in txt:
        fp=open(tmp_dir+"native.pdb",'w')
        fp.write(txt.split("ENDMDL")[0]+"ENDMDL\n")
        fp.close()
    
    # read model pdb
    fp=open(tmp_dir+"model.pdb",'rU')
    txt=fp.read()
    fp.close()

    # reindex residue index
    if offset!=0:
        lines=txt.splitlines()
        txt=""
        for line in lines:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                resi=str(int(line[22:26])+offset)
                line=line[:22]+' '*(4-len(resi))+resi+line[26:]
            txt+=line+'\n'
        fp=open(tmp_dir+"model.pdb",'w')
        fp.write(txt)
        fp.close()

    # Check if model_pdb is multi-model NMR structure
    if check_nmr and ("\nMODEL " in txt or initdat_pattern.search(txt)):
        if "\nMODEL " in txt:
            MODEL_start_lst=[e.start()+1 for e in re.finditer('\nMODEL ',txt)]
            if txt.startswith("MODEL"):
                MODEL_start_lst=[0]+MODEL_start_lst
        else:
            MODEL_start_lst=[e.start()+1 for e in initdat_pattern.finditer(txt)]
        pdb_txt='' #txt[:MODEL_start_lst[0]]
        RMSD_lst=[]
        TMscore_lst=[]
        for idx,start in enumerate(MODEL_start_lst):
            if idx==len(MODEL_start_lst)-1:
                MODEL=txt[start:]
            else:
                MODEL=txt[start:MODEL_start_lst[idx+1]]
                
            if not "ATOM  " in MODEL:
                pdb_txt+=MODEL
            else:
                single_model_file=tmp_dir+str(idx)
                fp=open(single_model_file,'w')
                fp.write(MODEL)
                fp.close()

                pdb_txt_tmp,RMSD,TMscore=sup_TMalign(single_model_file,
                    tmp_dir+"native.pdb", algo,execpath,check_nmr=False)
                pdb_txt+=pdb_txt_tmp
                RMSD_lst.append(RMSD)
                TMscore_lst.append(TMscore)

                if idx<len(MODEL_start_lst)-1:
                    sys.stdout.write('\n')
                #if os.path.isfile(single_model_file):
                    #os.remove(single_model_file)
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)
        return (pdb_txt,RMSD_lst,TMscore_lst)
            

    matrix_txt=tmp_dir+"matrix.txt"
    command=execpath+' '+tmp_dir+"model.pdb"+' '+tmp_dir+"native.pdb"
    if algo=="TMalign":
        command+= " -m "+tmp_dir+"matrix.txt" \
             +" && cat "+tmp_dir+"matrix.txt"
    stdout,stderr= subprocess.Popen(command, shell=True, 
        stdout = subprocess.PIPE).communicate()
    RMSD_pattern=re.compile("RMSD\s*=\s*([.\d]+)")
    TMscore_pattern=re.compile("TM-score\s*=\s*([.\d]+)")
    try:
        RMSD=RMSD_pattern.findall(stdout)[-1]
        TMscore=TMscore_pattern.findall(stdout)[-1]
        sys.stdout.write("RMSD="+RMSD+"; TMscore="+TMscore+";")
    except:
        RMSD="0.00"
        TMscore="0.0000"
        sys.stdout.write("RMSD=0.00; TMscore=0.0000;")
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)
        return ('',[float(RMSD)],[float(TMscore)])

    pattern=re.compile("\s[123]"+"\s+[-]{0,1}[.\d]+"*4+"\s*\n")
    m=0
    t=[0,0,0]
    u=[[0,0,0] for e in range(3)]
    for m,line in enumerate(pattern.findall(stdout)):
        t[m],u[m][0],u[m][1],u[m][2]=map(float,
            re.split("\s+",line.strip())[1:])

    fp=open(tmp_dir+"model.pdb",'rU')
    pdb_lines=fp.read().splitlines()
    fp.close()

    if not m:
        return '\n'.join(pdb_lines)+'\n'
        
    
    ''' Code for rotating Chain_1 from (x,y,z) to (X,Y,Z):
    do i=1,L
      X(i)=t(1)+u(1,1)*x(i)+u(1,2)*y(i)+u(1,3)*z(i)
      Y(i)=t(2)+u(2,1)*x(i)+u(2,2)*y(i)+u(2,3)*z(i)
      Z(i)=t(3)+u(3,1)*x(i)+u(3,2)*y(i)+u(3,3)*z(i)
    enddo
    '''
    pdb_txt=''
    for idx,line in enumerate(pdb_lines):
        if not line.startswith("ATOM  ") and not line.startswith("HETATM"):
            pdb_txt+=line+'\n'
            continue
        x=float(line[30:38])
        y=float(line[38:46])
        z=float(line[46:54])
        X='%.3f'%(t[0]+u[0][0]*x+u[0][1]*y+u[0][2]*z)
        Y='%.3f'%(t[1]+u[1][0]*x+u[1][1]*y+u[1][2]*z)
        Z='%.3f'%(t[2]+u[2][0]*x+u[2][1]*y+u[2][2]*z)
        line=line[:30]+' '*(8-len(X))+X+' '*(8-len(Y))+Y+ \
                       ' '*(8-len(Z))+Z+line[54:]
        pdb_txt+=line+'\n'
        '''
COLUMNS        DATA TYPE       CONTENTS                            
-----------------------------------------------------------------------
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         Atom serial number.
13 - 16        Atom            Atom name.
17             Character       Alternate location indicator.
18 - 20        Residue name    Residue name.
22             Character       Chain identifier.
23 - 26        Integer         Residue sequence number.
27             AChar           Code for insertion of residues.
31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)       Occupancy.
61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
73 - 76        LString(4)      Segment identifier, left-justified.
77 - 78        LString(2)      Element symbol, right-justified.
79 - 80        LString(2)      Charge on the atom.
        '''

    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)

    return (pdb_txt,[float(RMSD)],[float(TMscore)])

def sup_matrix(model_pdb,matrix_file, offset=0):
    '''Superpose "model_pdb" to "native_pdb" using TMalign/TMscore/MMalign

    model_pdb  - moved pdb structure
    native_pdb - fixed pdb structure
    algo       - algorithm for supersition {"TMalign","TMscore","MMalign"}
    execpath   - path to TMalign/TMscore/MMalign executable
    offset     - add "offset" to residue index of model.
    check_nmr  - check whether "model_pdb" is multi-model NMR structure
    '''
    tmp_dir="/tmp/"+os.getenv("USER")+'/'+algo+ \
       str(random.randint(100000000,999999999))+'/'
    while(os.path.isdir(tmp_dir)):
        tmp_dir="/tmp/"+os.getenv("USER")+'/'+algo+ \
            str(random.randint(100000000,999999999))+'/'
    os.makedirs(tmp_dir)
    shutil.copy(model_pdb,tmp_dir+"model.pdb")
    shutil.copy(matrix_file,tmp_dir+"matrix.txt")

    # read model pdb
    fp=open(tmp_dir+"model.pdb",'rU')
    txt=fp.read()
    fp.close()

    # reindex residue index
    if offset!=0:
        lines=txt.splitlines()
        txt=""
        for line in lines:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                resi=str(int(line[22:26])+offset)
                line=line[:22]+' '*(4-len(resi))+resi+line[26:]
            txt+=line+'\n'
        fp=open(tmp_dir+"model.pdb",'w')
        fp.write(txt)
        fp.close()

    fp=open(tmp_dir+"matrix.txt",'rU')
    stdout=fp.read()
    fp.close()

    pattern=re.compile("^\s*[0123]"+"\s+[-]{0,1}[.\d]+"*4+"\s*$")
    m=-1
    t=[0,0,0]
    u=[[0,0,0] for e in range(3)]
    for line in stdout.splitlines():
        matchall=pattern.findall(line)
        if not matchall:
            continue
        m+=1
        t[m],u[m][0],u[m][1],u[m][2]=map(float,
            re.split("\s+",line.strip())[1:])

    fp=open(tmp_dir+"model.pdb",'rU')
    pdb_lines=fp.read().splitlines()
    fp.close()

    if m<=0:
        return '\n'.join(pdb_lines)+'\n'
        
    
    ''' Code for rotating Chain_1 from (x,y,z) to (X,Y,Z):
    do i=1,L
      X(i)=t(1)+u(1,1)*x(i)+u(1,2)*y(i)+u(1,3)*z(i)
      Y(i)=t(2)+u(2,1)*x(i)+u(2,2)*y(i)+u(2,3)*z(i)
      Z(i)=t(3)+u(3,1)*x(i)+u(3,2)*y(i)+u(3,3)*z(i)
    enddo
    '''
    pdb_txt=''
    for idx,line in enumerate(pdb_lines):
        if not line.startswith("ATOM  ") and not line.startswith("HETATM"):
            pdb_txt+=line+'\n'
            continue
        x=float(line[30:38])
        y=float(line[38:46])
        z=float(line[46:54])
        X='%.3f'%(t[0]+u[0][0]*x+u[0][1]*y+u[0][2]*z)
        Y='%.3f'%(t[1]+u[1][0]*x+u[1][1]*y+u[1][2]*z)
        Z='%.3f'%(t[2]+u[2][0]*x+u[2][1]*y+u[2][2]*z)
        line=line[:30]+' '*(8-len(X))+X+' '*(8-len(Y))+Y+ \
                       ' '*(8-len(Z))+Z+line[54:]
        pdb_txt+=line+'\n'
        '''
COLUMNS        DATA TYPE       CONTENTS                            
-----------------------------------------------------------------------
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         Atom serial number.
13 - 16        Atom            Atom name.
17             Character       Alternate location indicator.
18 - 20        Residue name    Residue name.
22             Character       Chain identifier.
23 - 26        Integer         Residue sequence number.
27             AChar           Code for insertion of residues.
31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)       Occupancy.
61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
73 - 76        LString(4)      Segment identifier, left-justified.
77 - 78        LString(2)      Element symbol, right-justified.
79 - 80        LString(2)      Charge on the atom.
        '''

    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)

    return pdb_txt

def superpose(model_pdb="mobile.pdb",native_pdb="target.pdb",
    algo="TMalign", execpath='',writePDB=True,offset=0):
    '''Superpose "model_pdb" to "native_pdb" '''
    if not execpath and algo!="matrix": # default path to executables
        if algo.lower().startswith("pymol") or algo=="super":
            execpath="pymol"
        else:
            execpath=algo

        execpath_cur=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),execpath)
        if os.path.isfile(execpath_cur):
            execpath=execpath_cur

    if algo in ("TMscore","TMalign","MMalign"):
        pdb_txt,RMSD_lst,TMscore_lst=sup_TMalign(
            model_pdb,native_pdb,algo,execpath,offset)
    elif algo.startswith("pymol-"):
        pdb_txt,RMSD_lst,TMscore_lst=sup_pymol(
            model_pdb,native_pdb,algo,execpath,offset)
    elif algo=="matrix":
        pdb_txt=sup_matrix(model_pdb,native_pdb,offset)

    model_pdb_basename=os.path.basename(model_pdb)
    ext_idx=model_pdb_basename.rfind('.')
    if ext_idx==-1:
        model_fit_pdb=model_pdb_basename+"_fit"
    else:
        model_fit_pdb=model_pdb_basename[:ext_idx]+"_fit" \
                     +model_pdb_basename[ext_idx:]
    if not writePDB:
        return ''
    fp=open(model_fit_pdb,'w')
    fp.write(pdb_txt)
    fp.close()
    return model_fit_pdb


if __name__=="__main__":
    algo="TMalign"
    execpath=''
    offset=0
    writePDB=True
    argv=[]
    for arg in sys.argv[1:]:
        if not arg.startswith('-'):
            argv.append(arg)
            continue
        if arg.startswith('-algo='):
            algo=arg[len("-algo="):]
        elif arg.startswith('-execpath='):
            execpath=arg[len("-execpath="):]
        elif arg.startswith('-writePDB='):
            writePDB=(arg[len("-writePDB="):].lower()=="true")
        elif arg.startswith("-offset="):
            offset=int(arg[len("-offset="):])
        else:
            sys.stderr.write("ERROR! Unknown argument "+arg+'\n')

    if len(argv)<2:
        sys.stderr.write(docstring)
        exit()

    native_pdb=argv[-1]
    for model_pdb in argv[:-1]:
        sys.stdout.write(superpose(
            model_pdb,native_pdb,algo,execpath,writePDB,offset)+'\n')
