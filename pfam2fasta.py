#!/usr/bin/env python
docstring='''pfam2fasta alignment.aln
    Convert pfam/selex/clustal format to Fasta format
'''
import sys

def pfam2fasta(pfam_txt="alignment.aln"):
    '''read pfam/selex/clustal text "pfam_txt". return text in fasta format
    '''
    fasta_txt=''
    fasta_dict=dict() # key is sequence name, value is sequence
    fasta_list=[]     # list of sequence name
    for line in pfam_txt.splitlines():
        if not line.strip() or line[0] in ' #\t':
            continue
        line=line.split()
        if not len(line)==2:
            continue
        header,sequence=line
        if header in fasta_list:
            fasta_dict[header]+=sequence.strip()
        else:
            fasta_list.append(header)
            fasta_dict[header]=sequence.strip()
    for header in fasta_list:
        fasta_txt+='>%s\n%s\n'%(header,fasta_dict[header])
    return fasta_txt

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    for arg in sys.argv[1:]:
        fp=open(arg,'rU')
        pfam_txt=fp.read()
        fp.close()
        sys.stdout.write(pfam2fasta(pfam_txt))
