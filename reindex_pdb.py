#!/usr/bin/env python
'''Rearrange residue number of single chain PDB so that residues index start 
from the assigned number or according to the sequence file.'''
import Bio.SeqIO, Bio.PDB
import sys,os

docstring='''reindex_pdb.py startindex infile.pdb outfile.pdb
    Rearrange residue number of `infile.pdb` so the first residue start from 
    `startindex`

reindex_pdb.py seq.fasta  infile.pdb outfile.pdb
    Rearrange residue number of `infile.pdb` according to `seq.fasta` so the 
    residue index is the same as that in `seq.fasta`
'''

if len(sys.argv)<3:
    print >>sys.stderr,docstring
    exit()
struct = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure(
    sys.argv[2],sys.argv[2])
model=struct[0]
chain=[c for c in model][0]

pdb_sequence=str(Bio.PDB.PPBuilder().build_peptides(chain)[0].get_sequence())

if os.path.isfile(sys.argv[1]):
    # One sequence FASTA
    ref_sequence=str(Bio.SeqIO.read(open(sys.argv[1],'rU'),"fasta").seq)
    ref_sequence=ref_sequence.replace(' ','').replace('-','').upper()
    startindex=ref_sequence.find(pdb_sequence)+1
else:
    try:
        startindex=int(sys.argv[1])
    except:
        print >>sys.stderr,docstring
        exit()

index=startindex
for residue in chain:
    if residue.id[0]==' ':
        residue.id=(residue.id[0],index,residue.id[-1])
        index+=1

class NonHetSelect(Bio.PDB.Select):
    def accept_residue(self,residue):
        return 1 if residue.id[0]==' ' else 0

io=Bio.PDB.PDBIO()
io.set_structure(chain)
outfile=sys.argv[3] if len(sys.argv)>3 else sys.stdout
io.save(outfile,NonHetSelect())
