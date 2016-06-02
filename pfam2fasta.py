#!/usr/bin/env python
'''Convert Pfam format to Fasta format'''
import Bio.SeqIO # parse FASTA file
import sys
if len(sys.argv)<2:
    print >>sys.stderr,'''pfam2fasta alignment.pfam
    Convert Non-interleaved Pfam format file to Fasta format file'''
    exit()

print '\n'.join(['>'+'\n'.join(line.split()) for line in \
    open(sys.argv[-1],'rU').read().splitlines()])
