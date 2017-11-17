#!/usr/bin/env python
docstring='''
split_states.py pdb.pdb
    split multi-model PDB file pdb.pdb into individual models:
    pdb_0001.pdb, pdb_0002.pdb, ...
'''
import sys

def split_states(infile="pdb.pdb"):
    ''' split multi-model PDB file pdb.pdb into individual models '''
    fp=open(infile,'rU')
    txt=fp.read()
    fp.close()

    PDBtxt_list=[]
    if not "ENDMDL" in txt:
        PDBtxt_list.append(txt)
    else:
        PDBtxt_list=['\n'.join(block.lstrip().splitlines()[1:])+'\n' \
            for block in txt.split("ENDMDL") \
            if "ATOM  " in block or "HETATM" in block]

    pdb_file_list=[]
    for m in range(len(PDBtxt_list)):
        suffix=str(m+1)
        pdb_file_list.append("%s_%s.pdb"%(infile.split('.')[0],
            '0'*(4-len(suffix))+suffix))
        fp=open(pdb_file_list[m],'w')
        fp.write(PDBtxt_list[m])
        fp.close()
    return pdb_file_list

if __name__=="__main__":
    if len(sys.argv)<=1:
        sys.stderr.write(docstring)
        exit()

    for arg in sys.argv[1:]:
        for f in split_states(arg):
            sys.stdout.write(f+'\n')
