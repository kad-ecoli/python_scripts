#!/usr/bin/env python
docstring='''fasta2len.py seq.fasta > seq.len
    print the length of each entry in the input fasta file
'''
import sys
def fasta2len(infile="seq.fasta"):
    '''calculat the length of each entry in the input fasta file "infile"'''
    len_txt=''
    if infile=='-':
        txt=sys.stdin.read()
    else:
        fp=open(infile,'rU')
        txt=fp.read()
        fp.close()
    for block in ('\n'+txt).split('\n>'):
        block=block.strip()
        if not block:
            continue
        lines=block.splitlines()
        header=lines[0].split()[0]
        sequence=''.join(lines[1:])
        len_txt+="%s\t%u\n"%(header,len(sequence))
    return len_txt

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    len_txt=fasta2len(infile=sys.argv[1])
    if len(sys.argv)==2:
        sys.stdout.write(len_txt)
    else:
        fp=open(sys.argv[2],'w')
        fp.write(len_txt)
        fp.close()
