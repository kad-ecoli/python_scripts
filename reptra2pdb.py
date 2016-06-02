#!/usr/bin/env python
# 2016-01-03 Chengxin Zhang
docstring='''
reptra2pdb.py seq.fasta rep1.tra1.bz2

Read I-TASSER/QUARK decoy set "rep1.tra1.bz2" and output PDB according to
amino acids sequence specified by FASTA sequence "seq.fasta"

This script generates hundreds of single model PDB files at current directory, 
named "rep1.tra1.00001.pdb" and so on. This script also generates one multimodel 
PDB file "rep1.tra1.ent"
'''
import sys,os
import bz2
import re

aa1to3 = { # 3 letter to 1 letter amino acid code conversion
    'A':'ALA', 'V':'VAL', 'F':'PHE', 'P':'PRO', 'M':'MET',
    'I':'ILE', 'L':'LEU', 'D':'ASP', 'E':'GLU', 'K':'LYS',
    'R':'ARG', 'S':'SER', 'T':'THR', 'Y':'TYR', 'H':'HIS',
    'C':'CYS', 'N':'ASN', 'Q':'GLN', 'W':'TRP', 'G':'GLY',
    #'B':'ASX', 'Z':'GLX', 'U':'SEC', 'O':'PYL', 'J':'UNK', 'X':'UNK'
    }

def readFastaOrRawSequence(infile="seq.fasta"):
    '''read a sequence from a Fasta file or a text file.'''
    br=open(infile,'rU').read()
    seq=br.split('>')[1].splitlines()[1:] if '>' in br else br.splitlines()
    return ''.join([line.strip() for line in seq])

def read_bz2(infile="rep1.tra1.bz2"):
    '''read bz2 text'''
    if infile.endswith(".bz2"):
        fp=bz2.BZ2File(infile,'rU')
    else:
        fp=open(infile,'rU')
    txt=fp.read()
    fp.close()
    return txt

def split_reptra(infile="rep1.tra1.bz2"):
    '''read I-TASSER/QUARK decoy set file "rep1.tra1.bz2"
    return a list of plain text for each decoy
    '''
    txt=read_bz2(infile)

    # pattern for start section of one decoy
    pattern=re.compile("\n\s*\d+\s+[-]{0,1}[.\d]+\s+\d+\s+\d+\n")
    # index for start sections of each decoy (except the first section)
    match_index=[e.start()+1 for e in pattern.finditer(txt)]

    decoy_txt_list=[txt[:match_index[0]]] + \
        [txt[match_index[idx]:match_index[idx+1]] \
            for idx in range(len(match_index)-1)]+ \
        [txt[match_index[-1]:]]

    return decoy_txt_list

def convert_single_decoy_to_PDB(decoy_txt='',sequence=""):
    '''read I-TASSER/QUARK decoy set file "rep1.tra1.bz2"
    return a list of plain text for each decoy
    '''
    if not decoy_txt.strip():
        sys.stderr.write("ERROR! Empty trajectory!\n")
        return ''

    pattern=re.compile("^\s+[-.\d]+\s+[-.\d]+\s+[-.\d]+$")
    pdb_txt=''
    aa_index=-1
    # add ATOM entries
    for line in decoy_txt.splitlines():
        if not pattern.match(line):
            pdb_txt+=line+'\n'
            continue
        if not line.strip():
            continue
        aa_index+=1
        serial_number=str(aa_index+1) # Atom/residue serial number
        x,y,z=[str(e) for e in re.findall("[-.\d]+",line)] # coordinates
        
        # description for PDB coordinate format can be found at
        # http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM
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

        record_name="ATOM  "
        atom_serial_number=' '*(5-len(serial_number))+serial_number
        atom_name="  CA "
        alternate_location_indicator=" "
        residue_name=aa1to3[sequence[aa_index]] if sequence else \
            "ALA" # 3 letter residue name, ALA if not specified
        chain_identifier=" A"
        residue_sequence_number=' '*(4-len(serial_number))+serial_number
        residue_insertion_code='    '
        X_coordinate=' '*(8-len(x))+x
        Y_coordinate=' '*(8-len(y))+y
        Z_coordinate=' '*(8-len(z))+z
        occupancy="  1.00"
        Bfactor="  0.00"
        segment_identifier="          "
        element_symbol=" C"
        charge="  "

        line=record_name                 + \
             atom_serial_number          + \
             atom_name                   + \
             alternate_location_indicator+ \
             residue_name                + \
             chain_identifier            + \
             residue_sequence_number     + \
             residue_insertion_code      + \
             X_coordinate                + \
             Y_coordinate                + \
             Z_coordinate                + \
             occupancy                   + \
             Bfactor                     + \
             segment_identifier          + \
             element_symbol              + \
             charge

        pdb_txt+=line+'\n'

    # add CONECT entries
    for line in pdb_txt.splitlines():
        if not line.strip() or not line.startswith("ATOM  "):
            continue
        '''
COLUMNS       DATA  TYPE      FIELD        DEFINITION
-------------------------------------------------------------------------
 1 -  6        Record name    "CONECT"
 7 - 11        Integer        serial       Atom  serial number
12 - 16        Integer        serial       Serial number of bonded atom
17 - 21        Integer        serial       Serial number of bonded atom
22 - 26        Integer        serial       Serial number of bonded atom
27 - 31        Integer        serial       Serial number of bonded atom
        '''
        serial_number=line[6:11]
        if int(serial_number)==1:
            continue
        serial_number_prev=str(int(serial_number)-1)
        serial_number_prev=' '*(5-len(serial_number_prev))+serial_number_prev

        record_name="CONECT"
        pdb_txt+=record_name       + \
                 serial_number_prev+ \
                 serial_number     + \
                 ' '*64+'\n'

    return pdb_txt

def reptra2pdb(decoy_set="rep1.tra1.bz2",sequence=""):
    '''Read I-TASSER/QUARK decoy set "rep1.tra1.bz2" and output PDB according
    to amino acids sequence'''
    decoy_txt_list=split_reptra(infile=decoy_set)

    basename=decoy_set.split(os.path.sep)[-1]
    if basename.endswith(".bz2"):
        basename=basename[:-len(".bz2")]
    
    multimodel_PDB_txt=''
    for decoy_index,decoy_txt in enumerate(decoy_txt_list):
        pdb_txt=convert_single_decoy_to_PDB(
            decoy_txt=decoy_txt,sequence=sequence)

        model_serial_number=str(decoy_index+1)
        model_serial_number=' '*(4-len(model_serial_number))+model_serial_number
        multimodel_PDB_txt+="MODEL     "+model_serial_number+' '*66+'\n' \
            +pdb_txt+"ENDMDL"+' '*74+'\n'

        suffix=str(decoy_index+1)
        suffix='0'*(5-len(suffix))+suffix
        filename=basename+'.'+suffix+".pdb"
        print filename

        fp=open(filename,'w')
        fp.write(pdb_txt)
        fp.close()
    multimodel_PDB_txt+="END"+' '*77+'\n'
    filename=basename+".ent"
    print filename
    fp=open(filename,'w')
    fp.write(multimodel_PDB_txt)
    fp.close()
    return decoy_index

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()
    elif len(sys.argv)<3:
        sys.stderr.write("WARNING! You have not specified sequence. " \
            "Use ALA for all residues\n")
        for decoy_set in sys.argv[1:]:
            reptra2pdb(decoy_set)
    else:
        sequence=readFastaOrRawSequence(sys.argv[1])
        for decoy_set in sys.argv[2:]:
            reptra2pdb(decoy_set, sequence)
