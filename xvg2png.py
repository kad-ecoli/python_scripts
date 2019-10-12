#!/usr/bin/env python
docstring='''
xvg2img.py gromacs.xvg gromacs.png
    convert gromacs format image gromacs.xvg to png format image gromacs.png
'''
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import numpy as np
import sys
import re

title_pat=re.compile('@\s+title\s+\"([\s\S]+?)\"\n')
yaxis_pat=re.compile('@\s+yaxis\s+label\s+\"([\s\S]+?)\"\n')
xaxis_pat=re.compile('@\s+xaxis\s+label\s+\"([\s\S]+?)\"\n')
legend_pat=re.compile('@\s+\S+\s+legend\s+\"([\s\S]+?)\"\n')

def xvg2img(infile,outfile):
    data=np.loadtxt(infile,comments=["#","@"])
    fp=open(infile,'rU')
    txt=fp.read()
    fp.close()
    
    title_list =title_pat.findall(txt)
    yaxis_list =yaxis_pat.findall(txt)
    xaxis_list =xaxis_pat.findall(txt)
    legend_list=legend_pat.findall(txt)

    plt.figure()
    for i in range(1,data.shape[1]):
        if len(legend_list)>=i:
            plt.plot(data[:,0],data[:,i],label=legend_list[i-1])
        else:
            plt.plot(data[:,0],data[:,i])

    if title_list:
        plt.title(title_list[0])
    if yaxis_list:
        plt.ylabel(yaxis_list[0])
    if xaxis_list:
        plt.xlabel(xaxis_list[0])
    if len(legend_list)>1:
        plt.legend(loc="best")
    plt.savefig(outfile,dpi=300)
    return

if __name__=="__main__":
    if len(sys.argv)<=1:
        sys.stderr.write(docstring)
        exit()
    xvg2img(sys.argv[1],sys.argv[2])
