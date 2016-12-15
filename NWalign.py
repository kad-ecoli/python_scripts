#!/usr/bin/env python
# Needleman-Wunsch dynamic programming for pairwise protein sequence alignment
# Copyright 2015 Chengxin Zhang @ Yang Zhang lab
# Inspired by NWAlign.java by Ren-Xiang Yan and NWalign.f by Yang Zhang
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

docstring='''
Pairwise sequence alignment by standard Needleman-Wunsch algorithm
    NWalign F1.fasta F2.fasta (align two sequences in fasta file)
    NWalign F1.pdb F2.pdb 1   (align two sequences in PDB file)
    NWalign F1.fasta F2.pdb 2 (align Sequence 1 in fasta and 2 in pdb)
    NWalign GKDGL EVADELVSE 3 (align two sequences typed by keyboard)
    NWalign GKDGL F.fasta 4   (align Sequence 1 by keyboard and 2 in fasta)
    NWalign GKDGL F.pdb 5     (align Sequence 1 by keyboard and 2 in pdb)
'''
import sys
gap_open=-11 # gap gapopen
gap_extn=-1 # gap gapext
Blosum62Matrix=[
# A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
[ 4,-1,-2,-2 ,0,-1,-1 ,0,-2,-1,-1,-1,-1,-2,-1 ,1 ,0,-3,-2 ,0,-2,-1 ,0,-4],#A
[-1 ,5 ,0,-2,-3 ,1 ,0,-2 ,0,-3,-2 ,2,-1,-3,-2,-1,-1,-3,-2,-3,-1 ,0,-1,-4],#R
[-2 ,0 ,6 ,1,-3 ,0 ,0 ,0 ,1,-3,-3 ,0,-2,-3,-2 ,1 ,0,-4,-2,-3 ,3 ,0,-1,-4],#N
[-2,-2 ,1 ,6,-3 ,0 ,2,-1,-1,-3,-4,-1,-3,-3,-1 ,0,-1,-4,-3,-3 ,4 ,1,-1,-4],#D
[ 0,-3,-3,-3 ,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],#C
[-1 ,1 ,0 ,0,-3 ,5 ,2,-2 ,0,-3,-2 ,1 ,0,-3,-1 ,0,-1,-2,-1,-2 ,0 ,3,-1,-4],#Q
[-1 ,0 ,0 ,2,-4 ,2 ,5,-2 ,0,-3,-3 ,1,-2,-3,-1 ,0,-1,-3,-2,-2 ,1 ,4,-1,-4],#E
[ 0,-2 ,0,-1,-3,-2,-2 ,6,-2,-4,-4,-2,-3,-3,-2 ,0,-2,-2,-3,-3,-1,-2,-1,-4],#G
[-2 ,0 ,1,-1,-3 ,0 ,0,-2 ,8,-3,-3,-1,-2,-1,-2,-1,-2,-2 ,2,-3 ,0 ,0,-1,-4],#H
[-1,-3,-3,-3,-1,-3,-3,-4,-3 ,4 ,2,-3 ,1 ,0,-3,-2,-1,-3,-1 ,3,-3,-3,-1,-4],#I
[-1,-2,-3,-4,-1,-2,-3,-4,-3 ,2 ,4,-2 ,2 ,0,-3,-2,-1,-2,-1 ,1,-4,-3,-1,-4],#L
[-1 ,2 ,0,-1,-3 ,1 ,1,-2,-1,-3,-2 ,5,-1,-3,-1 ,0,-1,-3,-2,-2 ,0 ,1,-1,-4],#K
[-1,-1,-2,-3,-1 ,0,-2,-3,-2 ,1 ,2,-1 ,5 ,0,-2,-1,-1,-1,-1 ,1,-3,-1,-1,-4],#M
[-2,-3,-3,-3,-2,-3,-3,-3,-1 ,0 ,0,-3 ,0 ,6,-4,-2,-2 ,1 ,3,-1,-3,-3,-1,-4],#F
[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4 ,7,-1,-1,-4,-3,-2,-2,-1,-2,-4],#P
[ 1,-1 ,1 ,0,-1 ,0 ,0 ,0,-1,-2,-2 ,0,-1,-2,-1 ,4 ,1,-3,-2,-2 ,0 ,0 ,0,-4],#S
[ 0,-1 ,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1 ,1 ,5,-2,-2 ,0,-1,-1 ,0,-4],#T
[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1 ,1,-4,-3,-2,11 ,2,-3,-4,-3,-2,-4],#W
[-2,-2,-2,-3,-2,-1,-2,-3 ,2,-1,-1,-2,-1 ,3,-3,-2,-2 ,2 ,7,-1,-3,-2,-1,-4],#Y
[ 0,-3,-3,-3,-1,-2,-2,-3,-3 ,3 ,1,-2 ,1,-1,-2,-2 ,0,-3,-1 ,4,-3,-2,-1,-4],#V
[-2,-1 ,3 ,4,-3 ,0 ,1,-1 ,0,-3,-4 ,0,-3,-3,-2 ,0,-1,-4,-3,-3 ,4 ,1,-1,-4],#B
[-1 ,0 ,0 ,1,-3 ,3 ,4,-2 ,0,-3,-3 ,1,-1,-3,-1 ,0,-1,-3,-2,-2 ,1 ,4,-1,-4],#Z
[ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2 ,0 ,0,-2,-1,-1,-1,-1,-1,-4],#X
[-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4 ,1]]#*
seqW="ARNDCQEGHILKMFPSTWYVBZX*" # Amino acide order in scoring matrix

def readFastaOrRawSequence(infile="F1.fasta"):
    '''read a sequence from a Fasta file or a text file.'''
    br=open(infile,'rU').read()
    seq=br.split('>')[1].splitlines()[1:] if '>' in br else br.splitlines()
    return ''.join([line.strip() for line in seq])

def readPDB(infile="F1.pdb"):
    '''read a sequence from a PDB file'''
    from pdb2fasta import pdb2seq
    header_list,sequence_list=pdb2seq(infile,
        PERMISSIVE="MSE",outfmt="PDB",allowX=True)
    return sequence_list[0]

def empty_matrix(N_rows, N_cols, fill=0):
    '''Return a matrix (list of lists) which has `N_rows` rows and `N_cols`
    columns. Each element of the matrix equals to `fill`'''
    return [[fill for e in range(N_cols)]
                  for e in range(N_rows)]

def print_matrix(matrix=[[]],f1='',f2=''):
    '''Print matrix in user friendly format. For debugging'''
    f1='^'+f1 if f1 else '^'+' '*len(matrix)
    f2=' ^'+f2
    print '\t'.join(list(f2))
    for i,row in enumerate(matrix):
        print '\t'.join([f1[i]]+[str(e) for e in row])

def aa2int(f1="GKDGL"):
    ''' Convert AA sequence to interger according to the following mapping
    A R N D C Q E G H I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
    0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
    '''
    return [[i for i in range(len(seqW)) if c==seqW[i]][0] for c in f1]

def print_alignment(sequenceA='',sequenceB=''):
    '''print pairwise sequence alignment'''
    sequenceM='' # identity string. colon for exact match.
    L_ali=0
    for i in range(len(sequenceA)):
        if sequenceA[i]==sequenceB[i]:
            sequenceM+=':'
        else:
            sequenceM+=' '
        if sequenceA[i]!='-' and sequenceB[i]!='-':
            L_ali+=1
    len2=len(sequenceB.replace('-',''))
    L_id=sequenceM.count(':')
    identity=1.*L_id/len2 # Sequence identity
    print "Aligned length:    %u"%(L_ali)
    print "Identical length:    %u"%(L_id)
    print "Sequence identity:    %.3f (=   %u/   %u)\n"%(identity,L_id,len2)
    print sequenceA
    print sequenceM
    print sequenceB
    print ''.join([str(i%10) for i in range(1,len(sequenceA)+1)])+'\n' #resi

def initialize_matrix(len1=0,len2=0,alternative=True):
    '''initialize Dynamic Programming matrices
    len1 - length of sequence 1
    len2 - length of sequence 2
    alternaive - whether to adopt alternative matrix initialization by wei zheng
    '''
    val=empty_matrix(len1+1,len2+1)
    preV=empty_matrix(len1+1,len2+1) # penalty score of horizontal long gap
    preH=empty_matrix(len1+1,len2+1) # penalty score of vertical long gap
    jpV=empty_matrix(len1+1,len2+1)
    jpH=empty_matrix(len1+1,len2+1)
    idir=empty_matrix(len1+1,len2+1,fill='')

    # fill first row/column of dynamic programming score matrix val,
    # and path matrix idir
    for i in range(1,len1+1):
        idir[i][0]='-'
        jpH[i][0]=i
    for j in range(1,len2+1):
        idir[0][j]='|'
        jpV[0][j]=j
    for i in range(1,len1+1):
        val[i][0]=gap_open+gap_extn*(i-1)
    for j in range(1,len2+1):
        val[0][j]=gap_open+gap_extn*(j-1)

    if alternative:
        init_min=-999999999 # -1-sys.max
        for i in range(1,len1+1):
            preH[i][0]=gap_open+gap_extn*(i-1)
        for j in range(1,len2+1):
            preV[0][j]=gap_open+gap_extn*(j-1)
        for i in range(len1+1):
            preV[i][0]=init_min
        for j in range(len2+1):
            preH[0][j]=init_min
    else:
        for i in range(1,len1+1):
            preV[i][0]=gap_open+gap_extn*(i-1)
        for j in range(1,len2+1):
            preH[0][j]=gap_open+gap_extn*(j-1)
    return val,preV,preH,jpV,jpH,idir

def calcualte_score_gotoh(f1="GKDGL",f2="EVADELVSE", imut=Blosum62Matrix, 
    gap_open=gap_open, gap_extn=gap_extn): # "imut" is scoring matrix
    '''Calculate the dynamic programming matrix using Gotoh algorithm.
    Return these matrices: idir, jpV,jpH
    idir - Path
        \ : match-mismatch
        - : horizontal gap (deletion)
        | : vertical gap   (insertion)
    jpV - horizontal long gap number
    jpH - vertical long gap number
    All matrices are in the size of [len(f1)+1]*[len(f2)+1]
    '''
    # H, V, jpV, jpH are variable names taken from NWAlign.java,
    # which in term derived from NWalign.f
    # Since Fortran array indexing is column-wise, as oppose to row-wise
    # array indexing in Python, the variable names do not actually reflect
    # whether they are "horizontal" or "vertical" in this script

    # initialize matrices
    len1=len(f1)
    len2=len(f2)
    score=empty_matrix(len1+1,len2+1) # score of align i in f1 to j in f2
    val,preV,preH,jpV,jpH,idir=initialize_matrix(len1,len2,alternative=True)

    # convert char sequence to int sequence (aa2int)
    seq1=[[i for i in range(len(seqW)) if c==seqW[i]][0] for c in f1]
    seq2=[[i for i in range(len(seqW)) if c==seqW[i]][0] for c in f2]
    # calculate score for local match-mismatch
    for i in range(len1):
        for j in range(len2):
            score[i+1][j+1]=imut[seq1[i]][seq2[j]]

    # fill val and idir
    for i in range(1,len1+1):
        for j in range(1,len2+1):
            # penalty of consective deletion
            preV[i][j]=max(val[i][j-1]+gap_open,preV[i][j-1]+gap_extn)
            jpV[i][j]=jpV[i][j-1]+1 if preV[i][j]==preV[i][j-1]+gap_extn else 1
            # penalty of consective insertion
            preH[i][j]=max(val[i-1][j]+gap_open,preH[i-1][j]+gap_extn)
            jpH[i][j]=jpH[i-1][j]+1 if preH[i][j]==preH[i-1][j]+gap_extn else 1

            D=val[i-1][j-1]+score[i][j]        # match-mismatch \
            V=preV[i][j]                       # deletion       -
            H=preH[i][j]                       # insertion      |
            val[i][j]=max(D,V,H)

            if val[i][j]==D:
                idir[i][j]='\\'
            if val[i][j]==V:
                idir[i][j]+='-'
            if val[i][j]==H:
                idir[i][j]+='|'

    return idir, jpV,jpH

def trace_back_gotoh(idir,jpV,jpH,f1="GKDGL",f2="EVADELVSE"):
    '''trace back dynamic programming path to dicipher pairwise alignment
    return a two-element list for the diciphered pairwise alignment'''
    len1=len(f1)
    len2=len(f2)
    # re-fill first row/column of path matrix idir for back tracing
    for i in range(1,len1+1):
        idir[i][0]='|'
    for j in range(1,len2+1):
        idir[0][j]='-'
    sequenceA=''
    sequenceB=''
    len1=len(f1)
    len2=len(f2)
    i=len1
    j=len2
    while(i+j): 
        gaplen=0
        
        if 'none' in idir[i][j]:
            pass
        elif '-'  in idir[i][j]:
            gaplen=jpV[i][j]
            sequenceA='-'*gaplen+sequenceA
            sequenceB=f2[-gaplen:]+sequenceB
            f2=f2[:-gaplen]
            j=j-gaplen
        elif '|'  in idir[i][j]:
            gaplen=jpH[i][j]
            sequenceA=f1[-gaplen:]+sequenceA
            sequenceB='-'*gaplen+sequenceB
            f1=f1[:-gaplen]
            i=i-gaplen
        elif '\\' in idir[i][j]:
            sequenceA=f1[-1]+sequenceA
            sequenceB=f2[-1]+sequenceB
            f1=f1[:-1]
            f2=f2[:-1]
            i=i-1
            j=j-1
    return sequenceA,sequenceB

def NeedlemanWunsch(f1="GKDGL",f2="EVADELVSE", imut=Blosum62Matrix,
    gap_open=gap_open, gap_extn=gap_extn): # "imut" is scoring matrix
    '''Needleman-Wunsch dynamic programming algorithm
    GapPenalty=gap_open+gap_extn*(k-1), where k is gap number'''
    f1=f1.upper().replace('-','').strip()
    f2=f2.upper().replace('-','').strip()

    # idir (DP path); jpV (Horizontal jump number); jpH (Vertical jump number)
    idir,jpV,jpH=calcualte_score_gotoh(f1, f2, imut, gap_open, gap_extn)
    # sequenceA (aligned f1); sequenceB (aligned f2)
    sequenceA,sequenceB=trace_back_gotoh(idir,jpV,jpH,f1,f2)
    print_alignment(sequenceA,sequenceB)
    return sequenceA,sequenceB

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()
    elif len(sys.argv)<4:
        input_mode='0'
    else:
        input_mode=sys.argv[3]

    # parse input sequence
    if   input_mode=='0': # align two sequences in fasta file
        f1=readFastaOrRawSequence(sys.argv[1])
        f2=readFastaOrRawSequence(sys.argv[2])
    elif input_mode=='1': # align two sequences in PDB file
        f1=readPDB(sys.argv[1])
        f2=readPDB(sys.argv[2])
    elif input_mode=='2': # align Sequence 1 in fasta and 2 in pdb
        f1=readFastaOrRawSequence(sys.argv[1])
        f2=readPDB(sys.argv[2])
    elif input_mode=='3': # align two sequences typed by keyboard
        f1=sys.argv[1].strip()
        f2=sys.argv[2].strip()
    elif input_mode=='4': # align Sequence 1 by keyboard and 2 in fasta
        f1=sys.argv[1].strip()
        f2=readFastaOrRawSequence(sys.argv[2])
    elif input_mode=='5': # align Sequence 1 by keyboard and 2 in fasta
        f1=sys.argv[1].strip()
        f2=readPDB(sys.argv[2])
    else:
        sys.stderr.write("ERROR! Unsupported file type "+input_mode+docstring)
        exit()
    
    print "\nLength of sequence 1:\t"+str(len(f1))+" ->"+sys.argv[1]
    print   "Length of sequence 2:\t"+str(len(f2))+" ->"+sys.argv[2]

    NeedlemanWunsch(f1,f2)
