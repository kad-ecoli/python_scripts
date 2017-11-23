#!/usr/bin/env python
docstring='''pdb2sarst.py pdb.pdb > sarst.fasta
    convert PDB file to SARST string using DSSP defined backbone torsion angles

options:
    -show_seq={true,false}   
        whether to show sequence

    -show_ss={8,3,1,0}    
        whether to show secondary structure.
        8 - DSSP 8 state secondary structure assignment
        3 - 3 state secondary structure assignment
        2 - show the number of helix, strand, coil, all residues 
            in the sequence header
        1 - just show whether it is random coil (0), all alpha (1), all 
            beta (2), alpha beta (3) protein in the sequence header
        0 - do not show secondary structure assignment

    -show_sarst={true,false}
        whether to calculate ramachandran code (SARST code)

    -show_chain={true,false}
        whether to show chain ID in sequence name

    -pulchra_path=./pulchra
        path to pulchra executable. By default it is guessed by
        location of this script. pulchra is used to construct full atom
        model from backbone model when the input structure contains too
        few atoms.

    -dssp_path=./dssp
        path to DSSP executable. By default it is guessed by
        location of this script

    -fold_type_cutoff=0.9
        cutoff to define alpha, beta, or alpha-beta protein
        0.9 - (default) protein is coil (c) if >90% of residues are coil
              protein is alpha(a) if >90% of non-coil residues are within helix
              protein is beta (b) if >90% of non-coil residues are within strand
              otherwise a protein is alpha-beta (ab)
'''
import sys,os
import subprocess
import random
import shutil

sarst_matrix=[
    "SSSIIIIQQQQQQQQZZZZZZZRRYYYYYYYYSSSS",
    "SSIIIIIQQQQFFFFZZZZZZZZYYYYYYYYYYSSS",
    "SIIIIIIILFFFFFFFZZZZZZZYYYYYYYYYYSSS",
    "SSIIIGGHHFFFFFFFZZZZZZZZYYYYYYYYYSSS",
    "SMIGGGGHHLLFFLFLZZZZZZZYYYYYYYYYZZZS",
    "MMGMGGHHHLLLLLLLLZZZZZZZYZYYYYYZZYSM",
    "MMMMMMMHHLLLLLLLZZZZZZZYYYZYYYZZYZZM",
    "MMMMMMMHHLLLLLLZZZZZZZZPYYYZZZZZZZMM",
    "MMMMMMMMLLLLLLLZZZZZZPPPPPYYYZYZZZZM",
    "MMMMMMMMMLLLLLZZZZZZZPPPPPZZZZZZZZZM",
    "MMMMMMMMMMLLLZZZZPZZPPPPPPZPZZZZZZZZ",
    "MMMMMMMMMNNNZZZZZZZZPPPPPPPPZZZZZZZM",
    "MMMMNNNNNNNNNZZZZZZZPPPPPPPZZPZWWWZM",
    "NNNNNNNNNNNKZZZZZZZPPPPPPPPPPPZZZZZW",
    "NNNNNNNNNNKNZZZZZZZPPPPPPPPPPWWZZZZZ",
    "WNNNNNNNNNKKZZZZZZZPPPPPPPPPWWWWZZZZ",
    "ZNNNNNNNNKKKKZZZZZZZPPPPPPPWWWWWWZZZ",
    "ZNNNNNNNNKKKKKZZZZZZZPPPPPPWWWWWWWZZ",
    "ZNNNNNNNKKKKDDDZZZZZZZPPPPWWWWWWWWWZ",
    "ZNNNNNNNKKKKDDDDZZZZZZZZPWWWWWWWWWWZ",
    "WVVNNNNNKKEDDDDDDZZZZZZZZWWWWWWWWWWZ",
    "VVVVVVNEEEBBDDDDDZZZZZZZWWWWWWWWWWWW",
    "VVVVVVVEEEEACCCCDZZZZZZZWWWWWWWWWZWV",
    "VVVVVVVEEEECCCTTTZZZZZZZWWWZWWWWZZWV",
    "VVVVVVVVEEEETTTTTTZZZZZRRWWZZZZZZZZV",
    "VVVVVVVVEEETTTTTTTZZZZZZRRZZZZZZZZZV",
    "VVVVVVVVVTTTTTTTTTZZZZZRRRRZZZZZZZZZ",
    "ZVVVVVVVVTTTTTTTTZZZZZRRRZRZZZZZZZZV",
    "VZVVVVVVVTTTTTTTZZZZZRRRRRZZZRZZZZZZ",
    "ZZZVVVVVVTTTTTTZZZZZRRRRRRRZZZZZZZZS",
    "ZZSVVVVQQQZZZZZZZZZRRRRRRRRRZZZZZZZZ",
    "SSZSSVQQQQQQZZZZZZZRRRRRRRRRZZZZZSSS",
    "SSSSSQQQQQQQQZZZZZZRRRRRRRRRRZYSSSSS",
    "SSSSSQQQQQQQQQZZZZZZRRRRRRRRYYYSSSSS",
    "SSSSSQQQQQQQQQQZZZZZZRRRRRRYYYYYSSSS",
    "SSSSIIQQQQQQQQQZZZZZZZRRRYYYYYYYSSSS",
]

def detectDSSP(dssp_path="dssp"):
    '''auto detect the location of dssp executable'''
    if not dssp_path:
        dssp_path="dssp"
        bindir=os.path.dirname(os.path.abspath(__file__))
        if os.path.isfile(os.path.join(bindir,"dssp")):
            dssp_path=os.path.join(bindir,"dssp")
        elif os.path.isfile(os.path.join(bindir,"mkdssp")):
            dssp_path=os.path.join(bindir,"mkdssp")
        elif os.path.isfile(os.path.join(bindir,"dsspcmbi")):
            dssp_path=os.path.join(bindir,"dsspcmbi")
    if os.path.isfile(dssp_path):
        dssp_path=os.path.abspath(dssp_path)
    return dssp_path

def detectPULCHRA(pulchra_path="pulchra"):
    '''auto detect the location of pulchra executable'''
    if not pulchra_path:
        pulchra_path="pulchra"
        bindir=os.path.dirname(os.path.abspath(__file__))
        if os.path.isfile(os.path.join(bindir,"pulchra")):
            pulchra_path=os.path.join(bindir,"pulchra")
    if os.path.isfile(pulchra_path):
        pulchra_path=os.path.abspath(pulchra_path)
    return pulchra_path

def runDSSP(infile="pdb.pdb",dssp_path="dssp",pulchra_path="pulchra"):
    '''run dssp executable "dssp_path" on PDB file "infile".
    if infile is C-alpha trace, use pulchra executable "pulchra_path" to
    reconstruct full atom structure and run dssp.
    '''
    cmd=dssp_path+" -i "+infile+"|grep -v '\.$'|grep -vP '\s+#'"
    p=subprocess.Popen(cmd,shell=True,
        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=p.communicate()

    fp=open(infile,'rU')
    CA_num=len([line for line in fp.read().split("\nEND")[0].splitlines(
        ) if (line.startswith("ATOM  ") or line.startswith("HETATM")
        ) and len(line)>=54 and line[12:16]==" CA " and line[16] in " A"])
    fp.close()

    #if stderr.startswith("DSSP could not be created due to an error:"):
    if len(stdout.splitlines())<0.9*CA_num:
        # extract chain ID
        chainID=' '
        fp=open(infile,'rU')
        for line in fp:
            if line.startswith("ATOM  ") and line[13:15]=="CA":
                chainID=line[21]
                break

        tmp_dir="/tmp/"+os.getenv("USER")+"/dssp2sarst"+str(random.randint(
            1000,9999))+os.path.basename(infile).split('.')[0]+'/'
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
        tmp_pdb=tmp_dir+'xxxx.pdb'
        shutil.copy(infile,tmp_pdb)
        pul_cmd=' '.join(['cd',tmp_dir,';',pulchra_path,'-epc xxxx.pdb'])

        subprocess.Popen(pul_cmd, stdout=subprocess.PIPE, shell=True,
                ).communicate()
        if os.path.isfile(tmp_dir+"xxxx.rebuilt.pdb"):
            infile=tmp_dir+"xxxx.rebuilt.pdb"
        elif os.path.isfile(tmp_dir+"rebuilt_xxxx.pdb"):
            infile=tmp_dir+"rebuilt_xxxx.pdb"
        elif os.path.isfile(tmp_dir+"pul_xxxx.pdb"):
            infile=tmp_dir+"pul_xxxx.pdb"
        else:
            sys.stderr.write(stderr)
            shutil.rmtree(tmp_dir)
            return stdout

        txt=''
        fp=open(infile,'rU')
        for line in fp:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                txt+=line[:21]+chainID+line[22:54]+"  1.00  0.00           "+line[13]+"  \n"
        fp.close()
        fp=open(infile,'w')
        fp.write(txt)
        fp.close()

        cmd=dssp_path+" -i "+infile+"|grep -v '\.$'|grep -vP '\s+#'"
        p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        stdout,stderr=p.communicate()
        shutil.rmtree(tmp_dir)
    return stdout

def parseDSSP(dssp_txt='',show_seq=True,show_ss=8,show_sarst=True, 
    fold_type_cutoff=0.9):
    '''parse cleaned dssp output. 
    output:
        dssp_dict: a dict whose key is chain ID and value
                   is a list of [sequence,ss,sarst]
        ss_type_dict: a dict whose key is chain ID and value is fold type
                   0 - random coil, 1 - all alpha
                   2 - all beta,    3 - alpha+beta alpha/beta
        ss_num_dict: a dict whose key is chain ID and value is a four element
                   list, for number of [helix, strand, coil, all] residues
    '''
    dssp_dict=dict()
    ss_type_dict=dict()
    ss_num_dict=dict()
    for line in dssp_txt.splitlines():
        if len(line)<115:
            continue
        chainID=line[11]
        AA=line[13]
        if  AA=='!': # chain break
            continue
        if not chainID in dssp_dict:
            dssp_dict[chainID]=['']*sum([show_seq,show_ss>0,show_sarst])
            ss_num_dict[chainID]=[0]*4

        if show_seq:
            dssp_dict[chainID][0]+=AA

        SS=line[16].replace(' ','C')
        if show_ss:
            if show_ss<=3:
                SS=SS.replace('I','H').replace('G','H'
                    ).replace('B','E'
                    ).replace('S','C').replace('T','C')
            ss_num_dict[chainID][0]+=(SS in 'IGH')
            ss_num_dict[chainID][1]+=(SS in 'BE')
            ss_num_dict[chainID][2]+=(SS in 'STC')
            ss_num_dict[chainID][3]+=1
            dssp_dict[chainID][show_seq]+=SS

        if show_sarst:
            PHI=float(line[103:109])
            PSI=float(line[109:115])
            if PHI!=360 and PSI!=360:
                PHI-=(PHI==180)*360
                PSI+=(PSI==-180)*360
                SARST_CODE=sarst_matrix[int((180-PSI)/10)][int((180+PHI)/10)]
            else:
                SARST_CODE='X'
            dssp_dict[chainID][-1]+=SARST_CODE
    for chainID in dssp_dict:
        ss_type_dict[chainID]=0
        if show_ss:
            helix_res_num=0.
            strand_res_num=0.
            for resn in dssp_dict[chainID][show_seq]:
                helix_res_num+=(resn in 'HIG')
                strand_res_num+=(resn in 'BE')
            L=1.*len(dssp_dict[chainID][show_seq])
            if ((L-helix_res_num-strand_res_num)/L<fold_type_cutoff):
                if 1.*helix_res_num/(
                    helix_res_num+strand_res_num)>=fold_type_cutoff:
                    ss_type_dict[chainID]=1
                elif 1.*strand_res_num/(
                    helix_res_num+strand_res_num)>=fold_type_cutoff:
                    ss_type_dict[chainID]=2
                else:
                    ss_type_dict[chainID]=3
            #if set('HIG').intersection(dssp_dict[chainID][show_seq]):
                #ss_type_dict[chainID]+=1
            #if set('BE').intersection(dssp_dict[chainID][show_seq]):
                #ss_type_dict[chainID]+=2
    return dssp_dict,ss_type_dict,ss_num_dict

if __name__=="__main__":
    show_seq=True
    show_ss=8
    show_sarst=True
    show_chain=True
    dssp_path=''
    pulchra_path=''
    fold_type_cutoff=0.9

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith('-show_seq='):
            show_seq=(arg[len("-show_seq="):].lower()=="true")
        elif arg.startswith('-show_ss='):
            show_ss=int(arg[len("-show_ss="):])
        elif arg.startswith('-show_sarst='):
            show_sarst=(arg[len("-show_sarst="):].lower()=="true")
        elif arg.startswith('-show_chain='):
            show_chain=(arg[len("-show_chain="):].lower()=="true")
        elif arg.startswith('-dssp_path='):
            dssp_path=arg[len("-dssp_path="):]
        elif arg.startswith('-pulchra_path='):
            pulchra_path=arg[len("-pulchra_path="):]
        elif arg.startswith('-fold_type_cutoff='):
            fold_type_cutoff=float(arg[len("-fold_type_cutoff="):])
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()
        else:
            argv.append(arg)
        

    if not argv:
        sys.stderr.write(docstring)
        exit()

    if not dssp_path:
        dssp_path=detectDSSP(dssp_path)
    if not pulchra_path:
        pulchra_path=detectPULCHRA(pulchra_path)

    for infile in argv:
        dssp_txt=runDSSP(infile,dssp_path,pulchra_path)
        dssp_dict,ss_type_dict,ss_num_dict=parseDSSP(dssp_txt,show_seq,
            show_ss,show_sarst, fold_type_cutoff)

        PDBID=os.path.basename(infile).split('.')[0]
        txt=''
        for chainID in dssp_dict:
            header='>'+PDBID+(':'+chainID)*show_chain
            if show_ss==1:
                txt+=header+'\t%d'%ss_type_dict[chainID]+'\n'+''.join(
                    [line+'\n' for line in dssp_dict[chainID
                    ][:show_seq]+dssp_dict[chainID][show_seq+1:]])
            elif show_ss==2:
                txt+=header+'\t'+'\t'.join(map(str,ss_num_dict[chainID])
                    )+'\n'+''.join([line+'\n' for line in dssp_dict[chainID
                    ][:show_seq]+dssp_dict[chainID][show_seq+1:]])
            else:
                txt+=header+'\n'+''.join(
                    [line+'\n' for line in dssp_dict[chainID]])
        sys.stdout.write(txt)
