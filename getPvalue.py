#!/usr/bin/env python
docstring='''
getPvalue.py data1.csv data2.csv
    Perform Student T-test on two datasets using two-sided Student T test
    by scipy.stats.ttest_rel or scipy.stats.ttest_ind

input files:
    data1.csv data2.csv - two lists of data points, one data point per line

options:
    -paired={true,false}
        whether to the two datasets are paired
    -equal_var={true,false}
        whether to assume same variance. only two unpaired dataset can have
        inequal variance.
    -two_tailed={true,false}
        whether to perform two tailed test. if false, null hypothesis is 
        data1<=data2, alternative is data1>data2.

output:
    t-statistics p-value
'''
import numpy as np
import scipy.stats
import sys

def getPvalue(a=[],b=[],paired=True,equal_var=True,two_tailed=True):
    '''perform student T-test on array a and b,
    paired     - whether to the two datasets are paired
    equal_var  - whether to assume same variance
    two_tailed - whether to perform two tailed test
    '''
    if paired:
        tstatistic,pvalue=scipy.stats.ttest_rel(a,b)
    else:
        tstatistic,pvalue=scipy.stats.ttest_ind(a,b,equal_var=equal_var)
    if not two_tailed:
        pvalue=pvalue*0.5
        if tstatistic<0:
            pvalue=1-pvalue
    return tstatistic,pvalue

if __name__=="__main__":
    paired=True     # related dataset
    equal_var=True  # equal variance
    two_tailed=True # two tailed test

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-paired="):
            paired=(arg[len("-paired="):].lower()=="true")
        elif arg.startswith("-equal_var="):
            equal_var=(arg[len("-equal_var="):].lower()=="true")
        elif arg.startswith("-two_tailed="):
            two_tailed=(arg[len("-two_tailed="):].lower()=="true")
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)!=2:
        sys.stderr.write(docstring)
        exit()

    tstatistic,pvalue=getPvalue(np.loadtxt(argv[0]),np.loadtxt(argv[1]),
        paired,equal_var,two_tailed)

    sys.stderr.write(str(tstatistic)+'\t'+str(pvalue)+'\n')
