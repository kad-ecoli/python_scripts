#!/usr/bin/env python
'''Convert Fasta format to Pfam format'''
import Bio.SeqIO # parse FASTA file
import sys
if len(sys.argv)<2:
    print >>sys.stderr,'''fasta2pfam [options] alignment.fasta
Convert Fasta format file to Pfam format file
Options:
    -seqnamelen=10 maximum sequence name length
'''
    exit()

entry_list=[e for e in Bio.SeqIO.parse(sys.argv[-1],"fasta")]

seqnamelen=[int(a.lstrip("-seqnamelen=")) for a in sys.argv[1:-1] \
    if a.startswith("-seqnamelen=")]
seqnamelen=seqnamelen[0] if seqnamelen else max([len(e.id) for e in entry_list])

for e in entry_list:
    print e.id[:seqnamelen].ljust(seqnamelen)+' '+str(e.seq)
    
