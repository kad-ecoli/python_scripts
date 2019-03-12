#!/usr/bin/env python
# 2016-01-28 Chengxin Zhang
docstring='''
distance_pdb.py [options] pdb.pdb
    Calculate residue distance bins in single chain PDB file "pdb.pdb"

Options:
    -atom={CA,CB} calculate distance between "CA" for all residues, or "CA" for 
        gly and "CB" for other 19 amino acids

    -outfmt={stat,list,dist} output format:
        "list": tab-eliminated list listing residue index for contact pairs
        "dist": tab-eliminated list listing residue distances for all pairs
        "stat": statistics on number of contacts at short/medm/long/all
                range and protein length L

    -range={all,short,medium,long} sequences seperation range x
        "all":     1<=x
        "short":   6<=x<12
        "medium": 12<=x<24
        "long":   24<=x   (most useful)
        (default): 6<=x


distance_pdb.py [options] pdb.pdb distance.map
    Calculate accuracy of sorted residue distance map "distance.map"
    according to single chain PDB file "pdb.pdb"

    "distance.map" should be of ResTriplet2 format

Options:
    -cutoff=8

    -atom={CB,CA} atom with which contact is considered.
        CB - CA for GLY, CB for other 19 AA
        CA - CA for all 20 AA

    -range={all,short,medium,long}

    -offset=0    add "offset" to residue index in predicted contact map
    -infmt={res,pdb} input format of contact map
        res     - ResTriplet2
        pdb     - pdb coordinate file
'''
import sys,os
import re
import gzip

from contact_pdb import calc_res_dist,calc_res_contact,calc_contact_num
#import numpy as np
#bins=np.array([0]+list(np.arange(2,22,20./64))+[22])
bins=[0]+[20./64*i+2 for i in range(64)]+[22]

# see http://predictioncenter.org/casp12/doc/rr_help.html for contact range
# defination. Note that NeBcon/NN-BAYES uses different defination for long
# range contact
short_range_def=6 # short_range_def <= separation < medm_range_def
medm_range_def=12 # medm_range_def  <= separation < long_range_def
long_range_def=24 # long_range_def  <= separation. 25 in NeBcon/NN-BAYES

def read_pseudo_distance_map(infile="model1.pdb",atom_sele="CB",
    sep_range=str(short_range_def),offset=0):
    res_dist_list=calc_res_dist(infile,atom_sele)
    res_con_list=calc_res_contact(res_dist_list,sep_range,cutoff=bins[-1])
    resi1,resi2,p=zip(*res_con_list)
    p,resi1,resi2=zip(*sorted(zip(p,resi1,resi2)))
    res_con_list=zip(resi1,resi2,p)
    return calc_res_bin(res_con_list)

def read_distance_map(infile="dist.distri",sep_range=str(short_range_def),
    offset=0):
    '''Read NN-BAYES or CASP RR format contact map. return them in a zipped 
    list with 3 fields for each residue pair. 1st field & 2nd filed are for 
    residue indices, and 3rd field is for euclidean distance.
    '''
    resi1=[] # residue index 1 list
    resi2=[] # residue index 2 list
    p=[] # distance bin of each residue pair
    fp=sys.stdin
    if infile!='-':
        if infile.endswith(".gz"):
            fp=gzip.open(infile,'rU')
        else:
            fp=open(infile,'rU')
    lines=fp.read().strip().splitlines()
    fp.close()
    pattern=re.compile('^(\d+)\s+(\d+)\s+[.\d]+\s+[.\d]+\s+[.\d]+\s+:' \
        +"\s+(\d+)"*len(bins)+'$')
    for line in lines:
        if not line.strip(): # skip empty lines
            continue
        match_list=pattern.findall(line.strip())
        if not match_list:
            continue
        data=map(float,match_list[0][2:])
        resi_idx1=int(match_list[0][0])+offset # residue index 1
        resi_idx2=int(match_list[0][1])+offset # residue index 2
        seperation=abs(resi_idx1-resi_idx2)

        if (sep_range=="short"  and not \
            short_range_def<=seperation<medm_range_def) or \
           (sep_range=="medium" and not \
            medm_range_def<=seperation<long_range_def) or \
           (sep_range=="long"   and not long_range_def<=seperation):
            continue
        elif not sep_range in ["all","short","medium","long"] \
            and seperation<int(sep_range):
            continue

        resi1.append(resi_idx1)
        resi2.append(resi_idx2)
        p.append(len(bins))

        max_prob=0
        for b in range(len(bins)):
            if data[b]>max_prob:
                p[-1]=int(b)
                max_prob=data[b]
    return zip(resi1,resi2,p)

def compare_res_contact(res_bin_list,res_pred_list):
    '''compare residue distance map "res_dist_list" calculate from pdb to 
    predicted residue distance map "res_pred_list". return the result in a zipped 
    list with 3 fields for each pair (in the same order as res_pred_list).
    1st field & 2nd filed are for residue indices, 3rd field for whether the
    distance prediction is correct.'''
    cmp_list=[]
    for i,j,b in res_pred_list:
        cmp_list.append((i,j, 1.*((i,j,b) in res_bin_list)+ \
         .75*(((i,j,b+1) in res_bin_list) or ((i,j,b-1) in res_bin_list))+ \
         .50*(((i,j,b+2) in res_bin_list) or ((i,j,b-2) in res_bin_list))+ \
         .25*(((i,j,b+3) in res_bin_list) or ((i,j,b-3) in res_bin_list))))
    return cmp_list

def calc_lnat_acc_distance(cmp_list,con_num_dict,sep_range=str(short_range_def)):
    '''Calculate residue distance accuracy using ouput if "compare_res_contact"
    and native contact number diction "con_num_dict" '''
    top_pred=dict() # top short, medm, long, all prediction
    
    if not sep_range in ["medium","long"]:
        top_pred["short"]=[res_pair for res_pair in cmp_list if \
            short_range_def<=abs(res_pair[0]-res_pair[1])<medm_range_def
            ][:con_num_dict["short"]]

    if not sep_range in ["short","long"]:
        top_pred["medm" ]=[res_pair for res_pair in cmp_list if \
            medm_range_def<=abs(res_pair[0]-res_pair[1])<long_range_def
            ][:con_num_dict["medm"]]

    if not sep_range in ["short","medium"]:
        top_pred["long" ]=[res_pair for res_pair in cmp_list if \
            long_range_def<=abs(res_pair[0]-res_pair[1])
            ][:con_num_dict["long"]]

    if not sep_range in ["short","medium","long"]:
        top_pred["all"  ]=cmp_list[:con_num_dict["all"]]

    ACC=dict() # accuracy
    for key in top_pred:
        ACC[key]=0
        if top_pred[key]:
            ACC[key]=1.*sum([e[2] for e in top_pred[key]])/con_num_dict[key]
    return ACC,top_pred

def calc_res_bin(res_con_list):
    ''' convert real value distances to distance bins '''
    res_bin_list=[]
    for index,(i,j,dist) in enumerate(res_con_list):
        for b,lb in enumerate(bins[:-1]):
            if lb<=dist and (dist<bins[b+1] or (b==len(bins)-2 and dist<=bins[b+1])):
                res_bin_list.append((i,j,b))
                break
    return res_bin_list

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    atom_sele="CB"
    cutoff=8
    outfmt="stat"
    sep_range=str(short_range_def) # "6"
    offset=0
    infmt="rr"
    file_list=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-atom="):
            atom_sele=arg[len("-atom="):]
        elif arg.startswith("-range="):
            sep_range=arg[len("-range="):]
            if sep_range=="medm":
                sep_range="medium"
        elif arg.startswith("-outfmt="):
            outfmt=arg[len("-outfmt="):]
        elif arg.startswith("-infmt="):
            infmt=arg[len("-infmt="):]
        elif arg.startswith("-offset="):
            offset=int(arg[len("-offset="):])
        elif arg.startswith("-") and len(arg)>1:
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()
        else:
            file_list.append(arg)
    if not file_list:
        sys.stderr.write(docstring+"\nERROR! No PDB file")
        exit()
    
    res_dist_list=calc_res_dist(file_list[0],atom_sele)
    res_con_list=calc_res_contact(res_dist_list,sep_range,cutoff=bins[-1])
    res_bin_list=calc_res_bin(res_con_list)

    L=map(list,zip(*res_dist_list))
    L=L[0]+L[1]
    L=max(L)-min(L)+1

    if len(file_list)==1: # calculate residue contact
        if outfmt=="dist":
            for res_pair in res_con_list:
                sys.stdout.write("%d\t%d\t%.1f\n"%(res_pair[0],res_pair[1],res_pair[2]))
        elif outfmt=="list":
            for res_pair in res_bin_list:
                sys.stdout.write("%d\t%d\t%d\n"%(res_pair[0],res_pair[1],res_pair[2]))
        elif outfmt.startswith("stat") or outfmt.startswith("lnat"):
            con_num_dict=calc_contact_num(res_con_list,L)
            key_list=["short","medm","long","all","L"]
            sys.stderr.write('\t'.join(key_list)+'\n')
            sys.stdout.write('\t'.join([str(con_num_dict[key]
                ) for key in key_list])+'\n')


    elif len(file_list)==2: # calculate distance prediction accuracy
        if infmt!="pdb":
            res_pred_list=read_distance_map(file_list[1], sep_range,offset)
        else:
            res_pred_list=read_pseudo_distance_map(file_list[1],
                atom_sele, sep_range, offset)
        cmp_list=compare_res_contact(res_bin_list,res_pred_list)

        con_num_dict=calc_contact_num(res_con_list,L)
        ACC,top_pred=calc_lnat_acc_distance(cmp_list,con_num_dict,sep_range)
        if sep_range == "short":
            key_list=["short"]
        elif sep_range == "medium":
            key_list=["medm"]
        elif sep_range == "long":
            key_list=["long"]
        else:
            key_list=["short", "medm", "long", "all"]
        sys.stderr.write('\t'.join(key_list)+'\n')
        sys.stdout.write('\t'.join(['%.3f'%ACC[key] for key in key_list]
            )+'\n')
    else:
        sys.stderr.write(docstring+"ERROR! Too many arguments.\n")
