#!/usr/bin/env python
# 2016-04-27 Chengxin Zhang
docstring='''
pdb2fasta.py pdb.pdb > seq.fasta
    convert PDB file pdb.pdb to sequence FASTA file seq.fasta

options:
-PERMISSIVE={MSE,ATOM,HETATM} how to treat nonstandatd amino acids
    MSE   - (default), only ATOM or MSE HETATM in first MODEL is converted
    ATOM  - only allow ATOM residues
    HETATM- All all ATOM & HETATM residues with CA atom, including ligands

-allowX={true,true} whether to allow amino acid "X"
    true  - (default) allow any amino acid type, including X
    false - only allow amino acid that can mapped to 20 standard amino
            acids ACDEFGHIKLMNPQRSTVWY

-outfmt={PDB,COFACTOR} how to treat multichain PDB
    PDB      - convert multiple chains PDB to multiple sequence FASTA
    COFACTOR - convert PDB to single sequence FASTA. If there are 
               multiple chains in pdb, they will be treated as one chain.

-SEQRES={false,true} whether to convert from "SEQRES" entries
    false - (default) always convert from "ATOM" entries
    true  - convert from "SEQRES" entries if present, 
            otherwise convert from "ATOM" entries.
            (For PDB format input only)

-mol={all,protein,rna,dna} which macromolecule type to use
    all     - use " CA " or " C3'" atoms
    protein - use " CA " atoms
    rna     - use " C3'" atoms
    dna     - use " C3'" atoms, the same as rna
'''
import sys,os
import shutil
import textwrap
from fixMSE import code_with_modified_residues
import gzip,tarfile
import random

code_standard = {
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
    'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
    'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
    'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
    #'MSE':'M',
    ' DA':'a', ' DC':'c', ' DG':'g', ' DT':'t',
    '  A':'a', '  C':'c', '  G':'g', '  T':'u',
    }

def pdbbundle2seq(tarball_name="pdb-bundle.tar.gz",PERMISSIVE="MSE",
    outfmt="PDB",allowX=True,mol="all"):
    '''convert best effort/minimum PDB bundle to sequence
    '''
    chain_id_mapping=dict()
    tarball_prefix=os.path.basename(tarball_name).split('.')[0]
    tar=tarfile.open(tarball_name,'r:gz')
    names=tar.getnames()
    PDBid=names[-1].split('-')[0]

    # parse chain id mapping
    fp=tar.extractfile(PDBid+"-chain-id-mapping.txt")
    map_txt=fp.read()
    fp.close()
    names=[] # list of PDB files in tarball
    for section in map_txt.split('\n'+PDBid+"-pdb-bundle"):
        if not ':' in section:
            continue
        idx,section=section.split(".pdb:\n")
        pdb_bundle_name=PDBid+"-pdb-bundle"+idx+".pdb"
        names.append(pdb_bundle_name)
        for line in section.splitlines():
            New_chain_ID,Original_chain_ID=line.split()
            key=pdb_bundle_name.split('.')[0]+':'+New_chain_ID
            value=tarball_prefix+':'+Original_chain_ID
            chain_id_mapping[key]=value

    header_list=[]
    sequence_list=[]
    for pdb_bundle_name in names:
        # parse text in *-pdb-bundle*.pdb
        fp=tar.extractfile(pdb_bundle_name)
        txt=fp.read()
        fp.close()
        header_list_tmp,sequence_list_tmp=pdbtxt2seq(
            txt,pdb_bundle_name,PERMISSIVE,outfmt,allowX,False,mol)
        if outfmt=="PDB":
            header_list+=[chain_id_mapping[h] for h \
                in header_list_tmp]
        sequence_list+=sequence_list_tmp
    if outfmt=="COFACTOR":
        sequence=''.join([''.join(s.splitlines()) for s in sequence_list])
        sequence=textwrap.fill(''.join(sequence),60)
        header=tarball_prefix+'\t'+str(len(sequence))
        header_list=[header]
        sequence_list=[sequence]
    return header_list,sequence_list

def pdb2seq(infile="pdb.pdb", PERMISSIVE="MSE", outfmt="PDB",
    allowX=True, SEQRES=False,mol="all"):
    '''Convert PDB to sequence.
    Return two lists, one for headers and the other for sequence.

    PERMISSIVE - whether allow non-standard residues
        ATOM:   Only allow ATOM residues
        HETATM: Allow all ATOM & HETATM residues, even if they are ligands
        MSE:   (default) Disallow any non-standard amino acid apart from MSE
    '''
    if infile.endswith(".tar.gz"): # best effort/minimum PDB bundle
        return pdbbundle2seq(infile,PERMISSIVE,outfmt,allowX,mol)
    #elif infile.endswith(".cif") or infile.endswith(".cif.gz"): # PDBx/mmCIF
        #return mmCIF2seq(infile,PERMISSIVE,outfmt,allowX)
    elif infile.endswith(".gz"):
        fp=gzip.open(infile,'rU')
    else:
        fp=open(infile,'rU')
    txt=fp.read()
    fp.close()
    return pdbtxt2seq(txt,infile,PERMISSIVE,outfmt,allowX,SEQRES,mol)

def mmCIF2seq(infile="pdb.cif", PERMISSIVE="MSE",outfmt="PDB",
    allowX=True):
    ''' Convert mmCIF/PDBx format file "infile" '''
    if infile.endswith(".gz"):
        fp=gzip.open(infile,'rU')
    else:
        fp=open(infile,'rU')
    txt=fp.read()
    fp.close()

    aa3to1=code_with_modified_residues

    chain_list=[]
    chain_dict=dict() # Each chain will be one keu
    for line in txt.splitlines():
        if not line.startswith("ATOM") or line.startswith("HETATM"):
            continue
        
        line=line.split()
        if len(line)<26:
            continue

        group_PDB,atom_num,type_symbol,label_atom_id,label_alt_id, \
        label_comp_id,label_asym_id,label_entity_id,label_seq_id,  \
        pdbx_PDB_ins_code,Cartn_x,Cartn_y,Cartn_z,occupancy,       \
        B_iso_or_equiv,Cartn_x_esd,Cartn_y_esd,Cartn_z_esd,        \
        occupancy_esd,B_iso_or_equiv_esd,pdbx_formal_charge,       \
        auth_seq_id,auth_comp_id,auth_asym_id,auth_atom_id,        \
        pdbx_PDB_model_num=line

        if int(pdbx_PDB_model_num)>1:
            break # just parse the first model
        if label_atom_id!="CA" or not label_alt_id in ['.','A']: # just CA
            continue

        if PERMISSIVE=="ATOM" and group_PDB!="ATOM":
            continue
        elif PERMISSIVE=="MSE" and group_PDB!="ATOM" and label_comp_id!="MSE":
            continue

        aa=aa3to1[label_comp_id] if label_comp_id in aa3to1 else 'X'
        if not allowX and aa=='X':
            continue

        if not label_asym_id in chain_dict:
            chain_dict[label_asym_id]=[]
            chain_list.append(label_asym_id)

        residue_tuple=(label_seq_id,aa)
        if not residue_tuple in chain_dict[label_asym_id]:
            chain_dict[label_asym_id].append(residue_tuple)

    header_list=[]
    sequence_list=[]
    PDBID=os.path.basename(infile).split('.')[0]
    for chain_id in chain_list:
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

def pdbtxt2seq(txt='',infile='pdb.pdb',PERMISSIVE="MSE",outfmt="PDB",
    allowX=True,SEQRES=False,mol="all"):
    '''Convert PDB text "txt" to sequence read from PDB file "infile"
    Return two lists, one for headers and the other for sequence.

    PERMISSIVE - whether allow non-standard residues
        ATOM:   Only allow ATOM residues
        HETATM: Allow all ATOM & HETATM residues, even if they are ligands
        MSE:   (default) Disallow any non-standard amino acid apart from MSE
    '''
    txt=txt.split("\nENDMDL")[0] # Only the first model
    if not "SEQRES" in txt:
        SEQRES=False # use "ATOM" is "SEQRES" is absent
    mol=mol.lower()

    aa3to1=code_with_modified_residues
    
    chain_list=[]
    chain_dict=dict() # Each chain will be one key
    for line in txt.splitlines():
        line=line+' '*(80-len(line)) # Each line contains at least 80 char

        if SEQRES:
            if not line[:6]=="SEQRES":
                continue
            chain_id=line[11].replace(' ','_')
            tmp_seq=[]
            for residue in line[19:].split():
                if len(residue)!=3:
                    continue # only convert amino acid
                aa=aa3to1[residue] if residue in aa3to1 else 'X'
                if allowX or aa!='X':
                    tmp_seq.append((len(tmp_seq),aa))
            if tmp_seq:
                if not chain_id in chain_dict:
                    chain_dict[chain_id]=[]
                    chain_list.append(chain_id)
                chain_dict[chain_id]+=tmp_seq
            continue

        if (mol=="protein" and line[12:16]!=" CA ") or \
           (mol in {"rna","dna"} and line[12:16]!=" C3'") or \
           (mol=="all" and line[12:16]!=" CA " and line[12:16]!=" C3'"):
            continue
        if not line[16] in [' ','A']: # remove alternative location
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
        res_num=line[22:27] # residue sequence number
        aa=aa3to1[residue] if residue in aa3to1 else 'X' # one letter AA name
        if not allowX and aa=='X':
            continue
        residue_tuple=(res_num,aa)

        if not chain_id in chain_dict:
            chain_dict[chain_id]=[]
            chain_list.append(chain_id)

        if not residue_tuple in chain_dict[chain_id]:
            chain_dict[chain_id].append(residue_tuple)

    header_list=[]
    sequence_list=[]
    PDBID=os.path.basename(infile).split('.')[0]
    #PDBID=os.path.basename(infile).upper().split('.')[0]
    for chain_id in chain_list:
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


def pdb2fasta(infile="pdb.pdb", PERMISSIVE="MSE", outfmt="PDB",
    allowX=True,SEQRES=False, mol="all"):
    '''Convert PDB to FASTA'''
    header_list,sequence_list=pdb2seq(infile,PERMISSIVE,outfmt,allowX,SEQRES,mol)
    fasta_list=['>'+header_list[i]+'\n'+ \
                  sequence_list[i] for i in range(len(sequence_list))]
    return '\n'.join(fasta_list)+'\n'

if __name__=="__main__":
    PERMISSIVE="MSE"
    outfmt="PDB"
    allowX=True
    SEQRES=False
    mol="all"
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-outfmt="):
            outfmt=arg[len("-outfmt="):].upper()
        elif arg.startswith("-PERMISSIVE="):
            PERMISSIVE=arg[len("-PERMISSIVE="):].upper()
        elif arg.startswith("-allowX="):
            allowX=(arg[len("-allowX="):].lower()=="true")
        elif arg.startswith("-SEQRES="):
            SEQRES=(arg[len("-SEQRES="):].lower()=="true")
        elif arg.startswith("-mol="):
            mol=arg[len("-mol="):].lower()
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)<1:
        sys.stderr.write(docstring)
    
    for pdb in argv:
        sys.stdout.write(pdb2fasta(
            pdb, PERMISSIVE=PERMISSIVE, outfmt=outfmt, allowX=allowX,
            SEQRES=SEQRES, mol=mol))
