#!/usr/bin/env python
# 2015-11-15 Chengxin Zhang
docstring='''
split_fasta.py seq.txt seq.fasta

    Split Multiple-Sequence-FASTA file "seq.txt" into single-sequence FASTA
    files, whose name shall all be "seq.fasta".
    Output each FASTA file into one folder, whose name is the same as the 
    header (sequence name)

option:
    -batch_size=1
        how many sequences are there in one fasta file.
        if batch_size>1, split target sequence into mulitple-sequence FASTA
        files of "batch_size" sequences
    -exclude_list=list
        exclude entries listed in "list"
'''
import sys,os

def split_fasta(infile="seq.txt",outfile="seq.fasta",outdir='',
    batch_size=1,exclude_list=''):
    '''Split Multiple-Sequence-FASTA file into Single-Sequence-FASTA files.
    Return a list of headers. Options:
        infile  - Multiple sequence fasta for all input sequences
        outfile - a list containing target one sequence FASTA file name
        outdir  - (default: the same as dirname of "infile")
                  path to the target folders
        batch_size - number of sequence in each splitted FASTA file
        exclude_list - a file listing entries that are not included
    '''
    infile=os.path.abspath(infile)
    if not outdir:
        outdir=os.path.dirname(infile)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outfile_list=[outfile] if isinstance(outfile,str) else outfile

    exclude_set=set()
    if exclude_list and os.path.isfile(exclude_list):
        fp=open(exclude_list,'rU')
        exclude_set=set([line.split()[0] for line in \
            fp.read().splitlines() if line.strip()])
        fp.close()

    header_list=[]
    sequence_list=[]

    for entry in open(infile,'rU').read().split('>'): # parse by entry
        if not entry.strip():
            continue
        header=entry.strip().split()[0]
        if header in exclude_set:
            continue
        sequence='\n'.join(entry.splitlines()[1:])
        header_list.append(header)
        sequence_list.append(sequence)

    for idx,(header,sequence) in enumerate(zip(header_list,sequence_list)):
        if idx % batch_size == 0:
            if idx:
                fp.close()
            if not os.path.isdir(header):
                os.makedirs(header)
            fp=open(os.path.join(outdir,header,outfile),'w')
        fp.write('>'+header+'\n'+sequence+'\n')
    return header_list

if __name__=="__main__":
    batch_size=1
    exclude_list=''

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-batch_size="):
            batch_size=int(arg[len("-batch_size="):])
        elif arg.startswith("-exclude_list="):
            exclude_list=arg[len("-exclude_list="):]
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! Unknown option %s\n."%arg)
        else:
            argv.append(arg)

    if len(argv)<1:
        sys.stderr.write(docstring)
        exit()
    elif len(argv)<2:
        outfile="seq.fasta"
    else:
        outfile=argv[1]
    infile=argv[0]

    if batch_size<1:
        sys.stderr.write("ERROR! batch_size must be positive integer\n")
        exit()
    split_fasta(infile,outfile,
        batch_size=batch_size,exclude_list=exclude_list)
