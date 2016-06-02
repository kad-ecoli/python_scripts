#!/usr/bin/env python
# 2015-11-15 Chengxin Zhang
import sys,os
docstring='''
./split_fasta.py seq.txt seq.fasta

Split Multiple-Sequence-FASTA file "seq.txt" into Single-Sequence-FASTA
files, whose name shall all be "seq.fasta".
Output each FASTA file into one folder, whose name is the same as the 
header (sequence name)
'''

def split_fasta(infile="seq.txt",outfile="seq.fasta",outdir=''):
    '''Split Multiple-Sequence-FASTA file into Single-Sequence-FASTA files.
    Return a list of headers. Options:
        infile  - Multiple sequence fasta for all input sequences
        outfile - a list containing target one sequence FASTA file name
        outdir  - (default: the same as dirname of "infile")
                  path to the target folders
    '''
    infile=os.path.abspath(infile)
    if not outdir:
        outdir=os.path.dirname(infile)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outfile_list=[outfile] if isinstance(outfile,str) else outfile

    ss=[] # a list of headers

    for entry in open(infile,'rU').read().split('>'): # parse by entry
        if not entry:
            continue
        s=entry.split()[0] # only preserve the first word of header
        ss.append(s)

        sequence='\n'.join(entry.splitlines()[1:])
        if not os.path.isdir(s):
            os.makedirs(s)

        for outfile in outfile_list:
            seq_txt=open(os.path.join(outdir, s,outfile),'w')
            seq_txt.write('>'+s+'\n'+sequence+'\n')
            seq_txt.close()
    return ss

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()
    elif len(sys.argv)<3:
        outfile="seq.fasta"
    else:
        outfile=sys.argv[2]

    split_fasta(infile=sys.argv[1],
               outfile=sys.argv[2] if len(sys.argv)>2 else "seq.fasta")
