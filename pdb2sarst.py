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
        1 - just show whether it is random coil (0), all alpha (1), all 
            beta (2), alpha beta (3) protein in the sequence header
        0 - do not show secondary structure assignment

    -show_sarst={true,false}
        whether to calculate ramachandran code (SARST code)

    -show_chain={true,false}
        whether to show chain ID in sequence name
'''
import sys,os
import subprocess

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

def detectDSSP():
    '''auto detect the location of dssp executable'''
    dssp_exe="dssp"
    bindir=os.path.dirname(os.path.abspath(__file__))
    if os.path.isfile(os.path.join(bindir,"dssp")):
        dssp_exe=os.path.join(bindir,"dssp")
    elif os.path.isfile(os.path.join(bindir,"mkdssp")):
        dssp_exe=os.path.join(bindir,"mkdssp")
    elif os.path.isfile(os.path.join(bindir,"dsspcmbi")):
        dssp_exe=os.path.join(bindir,"dsspcmbi")
    return dssp_exe

def runDSSP(infile="pdb.pdb",dssp_exe="dssp"):
    '''run dssp executable "dssp_exe" on PDB file "infile" '''
    cmd=dssp_exe+" -i "+infile+"|grep -v '\.$'|grep -vP '\s+#'"
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    return stdout

def parseDSSP(dssp_txt='',show_seq=True,show_ss=8,show_sarst=True):
    '''parse cleaned dssp output. 
    output:
        dssp_dict: a dict whose key is chain ID and value
                   is a list of [sequence,ss,sarst]
        ss_type_dict: a dict whose key is chain ID and value is fold type
                   0 - random coil, 1 - all alpha
                   2 - all beta,    3 - alpha+beta alpha/beta
    '''
    dssp_dict=dict()
    ss_type_dict=dict()
    for line in dssp_txt.splitlines():
        if len(line)<115:
            continue
        chainID=line[11]
        AA=line[13]
        if  AA=='!': # chain break
            continue
        if not chainID in dssp_dict:
            dssp_dict[chainID]=['']*sum([show_seq,show_ss>0,show_sarst])

        if show_seq:
            dssp_dict[chainID][0]+=AA

        SS=line[16].replace(' ','C')
        if show_ss:
            if show_ss<=3:
                SS=SS.replace('I','H').replace('G','H'
                    ).replace('B','E'
                    ).replace('S','C').replace('T','C')
            dssp_dict[chainID][show_seq]+=SS

        if show_sarst:
            PHI=float(line[103:109])
            PSI=float(line[109:115])
            if PHI!=360 and PSI!=360:
                SARST_CODE=sarst_matrix[int((180-PSI)/10)][int((180+PHI)/10)]
            else:
                SARST_CODE='X'
            dssp_dict[chainID][-1]+=SARST_CODE
    for chainID in dssp_dict:
        ss_type_dict[chainID]=0
        if set('HIG').intersection(dssp_dict[chainID][show_seq]):
            ss_type_dict[chainID]+=1
        if set('BE').intersection(dssp_dict[chainID][show_seq]):
            ss_type_dict[chainID]+=2
    return dssp_dict,ss_type_dict

if __name__=="__main__":
    show_seq=True
    show_ss=8
    show_sarst=True
    show_chain=True

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
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()
        else:
            argv.append(arg)
        

    if not argv:
        sys.stderr.write(docstring)
        exit()

    dssp_exe=detectDSSP()
    for infile in argv:
        dssp_txt=runDSSP(infile,dssp_exe)
        dssp_dict,ss_type_dict=parseDSSP(dssp_txt,show_seq,show_ss,show_sarst)

        PDBID=os.path.basename(infile).split('.')[0]
        txt=''
        for chainID in dssp_dict:
            header='>'+PDBID+(':'+chainID)*show_chain
            if show_ss!=1:
                txt+=header+'\n'+''.join(
                    [line+'\n' for line in dssp_dict[chainID]])
            else:
                txt+=header+'\t%d'%ss_type_dict[chainID]+'\n'+''.join(
                    [line+'\n' for line in \
                    dssp_dict[chainID][:show_seq]+dssp_dict[chainID][show_seq+1:]])
        sys.stdout.write(txt)
