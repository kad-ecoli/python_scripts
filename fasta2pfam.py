#!/usr/bin/env python
docstring='''fasta2pfam [options] alignment.fasta
    Convert Fasta format file to Pfam format file

Options:
    -seqnamelen=10 maximum sequence name length
'''
import sys

def fasta2pfam(fasta_txt="seq.fasta",seqnamelen=0):
    '''read fasta text "fasta_txt". return text in pfam format
    seqnamelen - maximum sequence name length.
    '''
    pfam_txt=''
    for block in fasta_txt.split('>'):
        if not block.strip():
            continue
        header=block.split()[0]
        sequence=''.join(block.splitlines()[1:])
        if seqnamelen:
            header=header[:seqnamelen]
        pfam_txt+=header+' '+sequence+'\n'
    return pfam_txt

if __name__=="__main__":
    seqnamelen=0
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-seqnamelen="):
            seqnamelen=int(arg[len("-seqnamelen="):])
        else:
            argv.append(arg)

    if len(argv)<1:
        sys.stderr.write(docstring)
        exit()

    for arg in argv:
        if arg=='-':
            fasta_txt=sys.stdin.read()
        else:
            fp=open(arg,'rU')
            fasta_txt=fp.read()
            fp.close()
        sys.stdout.write(fasta2pfam(fasta_txt))
