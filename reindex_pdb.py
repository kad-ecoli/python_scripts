#!/usr/bin/env python
docstring='''
reindex_pdb.py startindex infile.pdb outfile.pdb
    Rearrange residue number of "infile.pdb" so the first residue start from 
    "startindex"

reindex_pdb.py seq.fasta  infile.pdb outfile.pdb
    Rearrange residue number of "infile.pdb" according to "seq.fasta" so the 
    residue index is the same as that in "seq.fasta"

if the structure contains any missing residue, residue number gaps will be removed
'''
import Bio.PDB
import sys,os

class NonHetSelect(Bio.PDB.Select): # class to select ATOM entries
    def accept_residue(self,residue):
        return 1 if residue.id[0]==' ' else 0

def read_single_sequence(fasta):
    '''read the first sequence from fasta file "fasta"'''
    fp=open(fasta,'rU')
    block_list=[b for b in fp.read().split('\n>') if b.strip()]
    fp.close()
    return ''.join(block_list[0].splitlines()[1:])

def reindex_pdb(startindex,infile,outfile):
    '''
    Read PDB file "infile", reindex residue number, and save it to "outfile".

    If startindex is interger, renumber residue number from 1

    If startindex is a sequence file, rearrange residue number according to 
    the fasta format sequence file.
    '''
    struct = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure(infile,infile)
    model=struct[0]
    chain=[c for c in model][0]

    pdb_sequence=str(Bio.PDB.PPBuilder().build_peptides(chain)[0].get_sequence())

    if os.path.isfile(startindex):
        ref_sequence=read_single_sequence(startindex)
        startindex=ref_sequence.find(pdb_sequence)+1
    elif startindex!="clean":
        try:
            startindex=int(sys.argv[1])
        except:
            sys.stderr.write(docstring)
            return

    index=startindex
    for residue in chain:
        if residue.id[0]==' ':
            residue.id=(residue.id[0],index,residue.id[-1])
            index+=1

    io=Bio.PDB.PDBIO()
    io.set_structure(chain)
    io.save(outfile,NonHetSelect())
    return
    
if __name__=="__main__":
    if len(sys.argv)<3:
        print >>sys.stderr,docstring
        exit()

    reindex_pdb(sys.argv[1], sys.argv[2], 
        sys.argv[3] if len(sys.argv)>3 else sys.stdout)
