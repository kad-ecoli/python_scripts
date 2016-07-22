#!/usr/bin/env python
# 2016-02-27 Chengxin Zhang
docstring='''
split_uniprot.py 4k5y.pdb
    try to split protion of PDB structures belonging to different uniprot
    accession. return the file name of splitted files.

split_uniprot.py -suffix={index,accession} 1a22.pdb.gz
    how are splitted files named.
    index    : (default) sequentially name splitted files 
               (i.e. 1a221.pdb, 1a222.pdb)
    accession: name splitted files by uniprot accessions
               (i.e. 4k5y_P34998.pdb, 4k5y_P00720.pdb)

split_uniprot.py -exclude_accession="P0AEX9,P0AEY0,P00720" 4k5y.pdb.gz
    exclude some uniprot ID in fusion protein. default is exclude maltose
    binding protein (P0AEX9,P0AEY0) and T4 Lysozyme (P00720)

split_uniprot.py -in_place={false,true} 1ank.pdb.gz
    how to name files if a PDB only one uniprot accessions.
    false: (default) name splitted file by suffix defined by "-suffix)
           (i.e. 1ank1.pdb or 1ank_P05082.pdb)
    true: do not give any suffix to splitted file (i.e. 1ank.pdb)
'''
import sys,os
import re
import gzip
import shutil

accession_pattern=re.compile(# uniprot accession (AC)
    "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")
#DBREF_pattern=re.compile("DBREF\s\s\d\w{3}\s(\w{0,1})\s+(\d+)\s+(\d+)\s+UNP\s+(\w+)\s+(\w+_\w+)\s+(\d+)\s+(\d+)\s*")
DBREF_pattern=re.compile("DBREF[\d]{0,1}\s{1,2}\d\w{3}\s(\w{0,1})\s+(\d+)\s+(\d+)\s+UNP\s+(\w+)")

chainID_column={ # column to store chainID, 0 for no chainID
    "DBREF":12, "SEQADV":16, "SEQRES":11, "MODRES":16,# Primary Structure
    "HET":12,                             # Heterogen
    "HELIX":19, "SHEET":21,               # Secondary Structure
    "SSBOND":15, "CONECT":0,              # Connectivity
    "SITE":22,                            # Miscellaneous
    "MODEL":0, "ATOM":21, "TER":21, "HETATM":21, "ENDMDL":0, # Coordinate
    "END":0,                              # Bookkeeping
}

resSeq_column={ # columns for residue sequence number, 0 for no resSeq
    "DBREF":(14,18), "SEQADV":(18,22), "MODRES":(18,22), # Primary Structure
    "HELIX":(21,25), "SHEET":(22,26),       # Secondary Structure 
    "SSBOND":(17,21), "CONECT":0,           # Connectivity
    "SITE":(23,27),                         # Miscellaneous
    "MODEL":0, "ATOM":(22,26), "TER":0, "HETATM":(22,26), "ENDMDL":0, # Coordinate
    "END":0,                                # Bookkeeping
}


def extract_uniprot(PDB_txt,accession):
    '''extract the portion of PDB for specific uniprot "accession"
    according to DBREF from PDB format text "PDB_txt"
    '''
    DBREF_list=DBREF_pattern.findall(PDB_txt)
    
    resSeq_range_dict=dict() # key is chain ID, value is range of residue index
    for DBREF in DBREF_list:
        if accession in DBREF[3:5]:
            chain=DBREF[0]
            resSeq_range=range(int(DBREF[1]),int(DBREF[2])+1)
            if chain in resSeq_range_dict:
                resSeq_range_dict[chain]+=resSeq_range
            else:
                resSeq_range_dict[chain]=resSeq_range

    PDB_lines=PDB_txt.splitlines()
    PDB_txt=''
    atom_serial_list=[]
    for line in PDB_lines:
        section=line.split()[0]
        if not section in chainID_column:
            continue
        if not chainID_column[section] and section!="CONECT":
            PDB_txt+=line+'\n'
            continue

        chain=line[chainID_column[section]]
        if not chain in resSeq_range_dict:
            continue
        resSeq_range=resSeq_range_dict[chain]

        if section in resSeq_column and (not resSeq_column[section] or \
            int(line[resSeq_column[section][0]:resSeq_column[section][1]])
            in resSeq_range):
            if section in ("ATOM","HETATM"):
                atom_serial_list.append(line[6:11])
            elif section=="CONECT" and not line[6:11] in atom_serial_list:
                continue # CONECT is after all ATOM & HETATM
            PDB_txt+=line+'\n'
    return PDB_txt

def split_uniprot(pdb_file,suffix='index',in_place=False,exclude_accession=''):
    '''extract PDB structure belonging to different uniprot accession in 
    "pdb_file". return the filename of splitted PDB file. 
    
    options:
    suffix - index: splitted files will be named in sequential order
             accession: splitted files will be named after uniprot accession
    exclude_accession - comma separated strings of uniprot accessions to 
            exclude if more than one uniprot accession is present
    in_place - what to do if only one uniprot accession is found
            if "True", the new file will not have the suffix
    '''
    filenames='' # file names of splitted 
    if pdb_file.endswith(".gz"):
        fp=gzip.open(pdb_file,'rU')
    else:
        fp=open(pdb_file,'rU')
    txt=fp.read()
    fp.close()
    filename=os.path.basename(pdb_file).split('.')[0]

    DBREF_list=DBREF_pattern.findall(txt)
    if not DBREF_list: # the PDB does not have any uniprot accession
        split_filename=filename+".pdb"
        fp=open(split_filename,'w')
        fp.write(txt)
        fp.close()
        filenames=split_filename+"\n"
        return filenames
    accession_set=set([e[3] for e in DBREF_list])
    if len(accession_set)>1 and exclude_accession:
        accession_set-=set(exclude_accession.split(','))

    for idx,accession in enumerate(accession_set):
        PDB_txt=extract_uniprot(txt,accession)
        if suffix=="index":
            split_filename=filename+str(idx+1)+".pdb"
        else:
            split_filename=filename+'_'+accession+".pdb"
        if len(accession_set)==1 and in_place==True:
            split_filename=filename+".pdb"
        fp=open(split_filename,'w')
        fp.write(PDB_txt)
        fp.close()
        filenames+=split_filename+'\n'
    return filenames

if __name__=="__main__":
    suffix="index"
    exclude_accession="P0AEX9,P0AEY0,P00720"
    in_place=False
    
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-suffix="):
            suffix=arg[len("-suffix="):]
        elif arg.startswith("-exclude_accession="):
            exclude_accession=arg[len("-exclude_accession="):]
        elif arg.startswith("-in_place="):
            in_place=(arg[len("-in_place="):].lower()=="true")
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if not len(argv):
        sys.stderr.write(docstring)
        exit()
    
    for arg in argv:
        sys.stdout.write(split_uniprot(arg,suffix,in_place,exclude_accession))
