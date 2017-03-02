#!/usr/bin/env python
docstring='''
seq_similarity [options] alignment.fasta
Generate Sequence Similarity Matrix for FASTA alignment file

Options:
    -seqname={none,index,name}   do not label column/row, label columns 
        using sequence "index", or using "header" of sequence

    -outfmt={matrix,list,first}  output format
        "matrix" - all against all, matrix format
        "list"   - all against all, first & second column for sequence index,
                   third column for sequence similarity
        "first"   - all against the first sequence

    -countgap={false,true}       whether counting gaps as residues. e.g.
        >seq1
        AA-TT-    Identity=2/3 if countgap==false, for sequence length==3
        >seq2     Identity=2/5 if countgap==true,  for sequence length==5
        ACAT--
    
    -eliminator={space,tab}      using "tab" or "space" as text eliminator
'''

import sys
import numpy
import numpy.matlib

def seq_similarity(infile="seq.fasta",countgap="false"):
    '''Calculate Sequence Similarity matrix from FASTA alignment file'''
    fp=open(infile,'rU')
    txt=fp.read()
    fp.close()
    seq_list=[''.join(block.splitlines()[1:]
        ) for block in txt.split('>') if block.strip()]
    len_list=list(map(len,seq_list))
    if max(len_list)!=min(len_list):
        sys.stderr.write("FATAL ERROR! Unaligned sequences.")
        exit()
    
    seq_num=len(seq_list)
    res_num=len(seq_list[0])

    int_array=numpy.array([numpy.array([ord(a) for a in line]) \
        for line in seq_list])

    Code_aa='ARNDCQEGHILKMFPSTWYV-';
    aa_array=dict()
    for aa in Code_aa:
        aa_array[aa]=numpy.zeros([seq_num,res_num])
        for ii in range(seq_num):
            aa_array[aa][ii][numpy.where(int_array[ii]==ord(aa))]=1

    id_matrix=numpy.zeros([seq_num,seq_num])
    for a in Code_aa:
        id_matrix+=aa_array[a].dot(aa_array[a].T)

    res_num_matrix=float(res_num)*numpy.ones([seq_num,seq_num])
    if   countgap.lower()=="true":
        gap_matrix=aa_array['-'].dot(aa_array['-'].T)
        id_matrix=(id_matrix-gap_matrix)/(res_num_matrix-gap_matrix)
    elif countgap.lower()=="false":
        gap_matrix=res_num_matrix-(1-aa_array['-']).dot((1-aa_array['-']).T)  
        id_matrix=(id_matrix-aa_array['-'].dot(aa_array['-'].T)
            )/(res_num_matrix-gap_matrix)
    
    return id_matrix

def get_header(infile="seq.fasta"):
    '''return a list of sequence header'''
    fp=open(infile,'rU')
    txt=fp.read()
    fp.close()
    return [block.split()[0] for block in txt.split('>') if block.strip()]

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    seqname="none"
    outfmt="matrix"
    countgap="false"
    eliminator="space"
    for arg in sys.argv[1:-1]:
        if arg.startswith("-seqname="):
            seqname=arg[len("-seqname="):]
        elif arg.startswith("-outfmt="):
            outfmt=arg[len("-outfmt="):]
        elif arg.startswith("-countgap="):
            countgap=arg[len("-countgap="):]
        elif arg.startswith("-eliminator="):
            eliminator=arg[len("-eliminator="):]
        else:
            sys.stderr.write("ERROR! Unknown argument "+arg)

    id_matrix=seq_similarity(sys.argv[-1],countgap)

    if   seqname=="index" in str(sys.argv[1:-1]):
        header=[str(e+1).ljust(10)+' ' for e in range(len(id_matrix))]
    elif seqname=="name" in str(sys.argv[1:-1]):
        header=[e+' ' for e in get_header(sys.argv[-1])] if eliminator=="tab" \
            else [e[:10].ljust(10)+' ' for e in get_header(sys.argv[-1])]
    elif seqname=="none":
        header=['' for e in range(len(id_matrix))]

    txt=''
    if outfmt=="list":
        for i in range(len(id_matrix)-1):
            for j in range(i+1,len(id_matrix)):
                txt+=(header[i]+header[j])+ \
                ('%.4f'%id_matrix[i][j]+'\n' if id_matrix[i][j]!=1 else "1\n")
        if eliminator=="tab":
            txt="Name_1\tName_2\tIdentity\n"+txt
    elif outfmt=="matrix":
        for i,line in enumerate(id_matrix):
            txt+=header[i]+'  '.join(
                [('%.4f'%e if e!=1 else "1        ") for e in line])+'\n'
        if seqname!="none":
            txt='Identity   '+''.join(header)+'\n'+txt
    elif outfmt=="first":
        for line in id_matrix:
            txt+='%.4f\n'%line[0]

    if eliminator=="tab":
        txt='\n'.join(['\t'.join(line.split()) for line in txt.splitlines()])+'\n'
    sys.stdout.write(txt)

