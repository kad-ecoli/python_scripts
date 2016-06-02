#!/usr/bin/env python
# 2016-04-27 Chengxin Zhang
import sys,os
import textwrap
docstring='''
pdb2fasta.py pdb.pdb > seq.fasta
    convert PDB file pdb.pdb to sequence FASTA file seq.fasta

options:
-PERMISSIVE={MSE,ATOM,HETATM} how to treat nonstandatd amino acids
    MSE   - (default), only ATOM or MSE HETATM in first MODEL is converted
    ATOM  - only allow ATOM residues
    HETATM- All all ATOM & HETATM residues with CA atom, including ligands

-outfmt={PDB,COFACTOR} how to treat multichain PDB
    PDB      - convert multiple chains PDB to multiple sequence FASTA
    COFACTOR - convert PDB to single sequence FASTA. If there are 
               multiple chains in pdb, they will be treated as one chain.
'''
code_standard = {
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
    'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
    'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
    'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
    #'MSE':'M',
    }

code_with_modified_residues = {
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M', 'ILE':'I',
    'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K', 'ARG':'R', 'SER':'S',
    'THR':'T', 'TYR':'Y', 'HIS':'H', 'CYS':'C', 'ASN':'N', 'GLN':'Q',
    'TRP':'W', 'GLY':'G',                                  'MSE':'M',
    '2AS':'D', '3AH':'H', '5HP':'E', 'ACL':'R', 'AIB':'A', 'ALM':'A', 
    'ALO':'T', 'ALY':'K', 'ARM':'R', 'ASA':'D', 'ASB':'D', 'ASK':'D', 
    'ASL':'D', 'ASQ':'D', 'AYA':'A', 'BCS':'C', 'BHD':'D', 'BMT':'T', 
    'BNN':'A', 'BUC':'C', 'BUG':'L', 'C5C':'C', 'C6C':'C', 'CCS':'C', 
    'CEA':'C', 'CHG':'A', 'CLE':'L', 'CME':'C', 'CSD':'A', 'CSO':'C', 
    'CSP':'C', 'CSS':'C', 'CSW':'C', 'CXM':'M', 'CY1':'C', 'CY3':'C', 
    'CYG':'C', 'CYM':'C', 'CYQ':'C', 'DAH':'F', 'DAL':'A', 'DAR':'R', 
    'DAS':'D', 'DCY':'C', 'DGL':'E', 'DGN':'Q', 'DHA':'A', 'DHI':'H', 
    'DIL':'I', 'DIV':'V', 'DLE':'L', 'DLY':'K', 'DNP':'A', 'DPN':'F', 
    'DPR':'P', 'DSN':'S', 'DSP':'D', 'DTH':'T', 'DTR':'W', 'DTY':'Y', 
    'DVA':'V', 'EFC':'C', 'FLA':'A', 'FME':'M', 'GGL':'E', 'GLZ':'G', 
    'GMA':'E', 'GSC':'G', 'HAC':'A', 'HAR':'R', 'HIC':'H', 'HIP':'H', 
    'HMR':'R', 'HPQ':'F', 'HSD':'H', 'HSE':'H', 'HSP':'H', 'HTR':'W', 
    'HYP':'P', 'IIL':'I', 'IYR':'Y', 'KCX':'K', 'LLY':'K', 'LTR':'W',
    'LYM':'K', 'LYZ':'K', 'MAA':'A', 'MEN':'N', 'MHS':'H', 'MIS':'S',
    'MLE':'L', 'MPQ':'G', 'MSA':'G', 'MVA':'V', 'NEM':'H', 'NEP':'H', 
    'NLE':'L', 'NLN':'L', 'NLP':'L', 'NMC':'G', 'OAS':'S', 'OCS':'C', 
    'OMT':'M', 'PAQ':'Y', 'PCA':'E', 'PEC':'C', 'PHI':'F', 'PHL':'F', 
    'PR3':'C', 'PRR':'A', 'PTR':'Y', 'SAC':'S', 'SAR':'G', 'SCH':'C', 
    'SCS':'C', 'SCY':'C', 'SEL':'S', 'SEP':'S', 'SET':'S', 'SHC':'C', 
    'SHR':'K', 'SOC':'C', 'STY':'Y', 'SVA':'S', 'TIH':'A', 'TPL':'W', 
    'TPO':'T', 'TPQ':'A', 'TRG':'K', 'TRO':'W', 'TYB':'Y', 'TYQ':'Y', 
    'TYS':'Y', 'TYY':'Y', 'AGM':'R', 'GL3':'G', 'SMC':'C', 'CGU':'E',
    'CSX':'C',
    'SEC':'U', 'PYL':'O', 'ASX':'B', 'GLX':'Z', 'LLP':'X', 'UNK':'X',
    }

def pdb2seq(infile="pdb.pdb", PERMISSIVE="MSE", outfmt="PDB"):
    '''Convert PDB to sequence.
    Return two lists, one for headers and the other for sequence.

    PERMISSIVE - whether allow non-standard residues
        ATOM:   Only allow ATOM residues
        HETATM: Allow all ATOM & HETATM residues, even if they are ligands
        MSE:   (default) Disallow any non-standard amino acid apart from MSE
    '''
    txt=open(infile,'rU').read()
    txt=txt.split("\nENDMDL")[0] # Only the first model

    aa3to1=code_with_modified_residues

    chain_dict=dict() # Each chain will be one key
    for line in txt.splitlines():
        line=line+' '*(80-len(line)) # Each line contains at least 80 char

        if line[13:15]!="CA": # Carbon Alpha Only
            continue

        residue=line[17:20] # residue name

        if   PERMISSIVE == "ATOM"   and line[0:6]!="ATOM  ":
            continue
        elif PERMISSIVE == "HETATM" and not line[0:6] in ("ATOM  ","HETATM"):
            continue
        elif PERMISSIVE == "MSE" and (not  line[0:6]=="ATOM  "  and \
             not (line[0:6]=="HETATM" and residue=="MSE")):
            continue
        
        # underscore for empty chain identifier
        chain_id=line[21].replace(' ','_')
        res_num=int(line[22:26]) # residue sequence number
        aa=aa3to1[residue] if residue in aa3to1 else 'X' # one letter AA name
        residue_tuple=(res_num,aa)

        if not chain_id in chain_dict:
            chain_dict[chain_id]=[]

        if not residue_tuple in chain_dict[chain_id]:
            chain_dict[chain_id].append(residue_tuple)

    header_list=[]
    sequence_list=[]
    PDBID=os.path.basename(infile).split('.')[0]
    #PDBID=os.path.basename(infile).upper().split('.')[0]
    for chain_id in sorted(chain_dict):
        res_num_list,sequence=zip(*chain_dict[chain_id])
        header_list.append(PDBID+':'+chain_id)
        sequence_list.append(''.join(sequence))

    if outfmt=="COFACTOR":
        if len(sequence_list)>1:
            sys.stderr.write("WARNING! Multichain PDB %s\n"%infile)
        sequence=''.join(sequence_list)
        header=PDBID+'\t'+str(len(sequence))
        sequence=textwrap.fill(''.join(sequence),60)
        header_list=[header]
        sequence_list=[sequence]
    return header_list,sequence_list

def pdb2fasta(infile="pdb.pdb", PERMISSIVE="MSE", outfmt="PDB"):
    '''Convert PDB to FASTA'''
    header_list,sequence_list=pdb2seq(infile,PERMISSIVE,outfmt)
    fasta_list=['>'+header_list[i]+'\n'+ \
                  sequence_list[i] for i in range(len(sequence_list))]
    return '\n'.join(fasta_list)

if __name__=="__main__":
    PERMISSIVE="MSE"
    outfmt="PDB"
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-outfmt="):
            outfmt=arg[len("-outfmt="):].upper()
        elif arg.startswith("-PERMISSIVE="):
            PERMISSIVE=arg[len("-PERMISSIVE="):].upper()
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)<1:
        sys.stderr.write(docstring)
    
    for pdb in argv:
        print pdb2fasta(pdb, PERMISSIVE=PERMISSIVE, outfmt=outfmt)
