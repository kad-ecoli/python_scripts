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

    -show_acc=0,1,2
        whether to print relative solvent accessibility
        0 - (default) do not print solvent accessibility
        1 - horizontally print relative solvent accessibility. 
            only first digit after decimal point
        2 - vertically print relative solvent accessibility
        3 - vertically print absolute solvent accessibility
        4 - show total relative solvent accessibility in header
        5 - show total absolute solvent accessibility in header
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

MaxASA_dict={
# data from Tien et al (2013). 
# "Maximum allowed solvent accessibilites of residues in proteins". 
# PLoS ONE. 8 (11): e80635. doi:10.1371/journal.pone.0080635. PMID 24278298.
#
# Tien (theor) ,Tien (emp), Miller 1987, Rose 1985, hhsuite, sann
    'A':[ 129.0, 121.0, 113.0, 118.1, 106.0, 115.0],
    'R':[ 274.0, 265.0, 241.0, 256.0, 248.0, 225.0],
    'N':[ 195.0, 187.0, 158.0, 165.5, 157.0, 160.0],
    'D':[ 193.0, 187.0, 151.0, 158.7, 163.0, 150.0],
    'C':[ 167.0, 148.0, 140.0, 146.1, 135.0, 135.0],
    'E':[ 223.0, 214.0, 183.0, 186.2, 194.0, 190.0],
    'Q':[ 225.0, 214.0, 189.0, 193.2, 198.0, 180.0],
    'G':[ 104.0,  97.0,  85.0,  88.1,  84.0,  75.0],
    'H':[ 224.0, 216.0, 194.0, 202.5, 184.0, 195.0],
    'I':[ 197.0, 195.0, 182.0, 181.0, 169.0, 175.0],
    'L':[ 201.0, 191.0, 180.0, 193.1, 164.0, 170.0],
    'K':[ 236.0, 230.0, 211.0, 225.8, 205.0, 200.0],
    'M':[ 224.0, 203.0, 204.0, 203.4, 188.0, 185.0],
    'F':[ 240.0, 228.0, 218.0, 222.8, 197.0, 210.0],
    'P':[ 159.0, 154.0, 143.0, 146.8, 136.0, 145.0],
    'S':[ 155.0, 143.0, 122.0, 129.8, 130.0, 115.0],
    'T':[ 172.0, 163.0, 146.0, 152.5, 142.0, 140.0],
    'W':[ 285.0, 264.0, 259.0, 266.3, 227.0, 255.0],
    'Y':[ 263.0, 255.0, 229.0, 236.8, 222.0, 230.0],
    'V':[ 174.0, 165.0, 160.0, 164.5, 142.0, 155.0],

# extended from above
    'J':[ 199.0, 193.0, 181.0,187.05, 166.5, 172.5],# (L+I)/2
    'U':[ 167.0, 148.0, 140.0, 146.1, 135.0, 135.0], # C
    'Z':[ 224.0, 214.0, 186.0, 189.7, 196.0, 185.0], # (Q+E)/2
    'B':[ 194.0, 187.0, 154.5, 162.1, 160.0, 155.0], # (D+N)/2
    'O':[ 201.0, 191.0, 180.0, 193.1, 205.0, 200.0], # K
    'X':[191.02,182.19,165.99,172.92, 180.0, 159.9], # average for all
}

#aa_scale_dict={
    #'A': 8.25, 'R': 5.53, 'N': 4.06, 'D': 5.45, 'C': 1.37,
    #'Q': 3.93, 'Q': 6.75, 'G': 7.07, 'H': 2.27, 'I': 5.96,
    #'L': 9.66, 'K': 5.84, 'M': 2.42, 'F': 3.86, 'P': 4.70,
    #'S': 6.56, 'T': 5.34, 'W': 1.08, 'Y': 2.92, 'V': 6.87,
#}

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
    fold_type_cutoff=0.9, show_acc=0):
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
    acc_dict=dict()
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
            if show_acc in [4,5]:
                acc_dict[chainID]=0
            else:
                acc_dict[chainID]=''

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

        if show_acc:
            ACC=float(line[35:38])
            if show_acc in [1,2,4]:
                ACC/=MaxASA_dict[AA][0]
            if show_acc==1:
                acc_dict[chainID]+=str(min([9,int(10*ACC)]))
            elif show_acc in [2,3]:
                acc_dict[chainID]+=str(ACC)+'\n'
            elif show_acc in [4,5]:
                acc_dict[chainID]+=ACC

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
    return dssp_dict,ss_type_dict,ss_num_dict,acc_dict

if __name__=="__main__":
    show_seq=True
    show_ss=8
    show_sarst=True
    show_chain=True
    show_acc=0
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
        elif arg.startswith('-show_acc='):
            show_acc=int(arg[len("-show_acc="):])
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
        dssp_dict,ss_type_dict,ss_num_dict,acc_dict=parseDSSP(dssp_txt,
            show_seq, show_ss,show_sarst, fold_type_cutoff, show_acc)

        PDBID=os.path.basename(infile).split('.')[0]
        txt=''
        for chainID in dssp_dict:
            header='>'+PDBID+(':'+chainID)*show_chain
            if show_ss==1:
                header+='\t%d'%ss_type_dict[chainID]
            elif show_ss==2:
                header+='\t'+'\t'.join(map(str,ss_num_dict[chainID]))

            if show_acc in [4,5]:
                header+="\t%.2f"%acc_dict[chainID]

            if show_ss in [1,2]:
                txt+=header+'\n'+''.join(
                    [line+'\n' for line in dssp_dict[chainID
                    ][:show_seq]+dssp_dict[chainID][show_seq+1:]])
            else:
                txt+=header+'\n'+''.join(
                    [line+'\n' for line in dssp_dict[chainID]])

            if show_acc in [1,2,3]:
                txt+=acc_dict[chainID]+'\n'*(show_acc==1)
        sys.stdout.write(txt)
