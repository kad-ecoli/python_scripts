#!/usr/bin/env python
docstring='''
batch_pdb_image.py list output
    Use pymol to generate gallery view of pdb files listed by "list",
    and output png format image to "output.*.png"
'''
import sys, os, subprocess
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
from string import Template

pymol_template=Template(';'.join([
    'pymol -c -q -d "load $infile, infile',
    'viewport 640, 640',
    'zoom',
    #'center infile',
    #'zoom center, 20',
    'bg_color white',
    'hide all',
    'show cartoon',
    'spectrum',
    'set ray_shadow, 0',
    'set ray_opaque_background, off',
    'ray 640, 640', # by default, this is 640x480
    'png $outfile"']))
res_count_template=Template('grep " CA " $infile|grep "^ATOM"|wc -l')

def batch_pdb_image(infile_list,outfile):
    outfile_list=[]
    plt_count=0
    ncol=4
    nrow=5
    gallery_count=1
    fontsize=10
    plt.figure(gallery_count,figsize=(8.27,11.69))
    for infile in infile_list:
        plt_count+=1
        print(infile)
        cmd=pymol_template.substitute(
            infile=infile,
            outfile=outfile+".png",
            )
        p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        p.communicate()
        cmd=res_count_template.substitute(infile=infile)
        p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        Lch=int(p.communicate()[0])
        plt.subplot(nrow,ncol,plt_count)
        img=plt.imshow(plt.imread(outfile+".png"))
        os.system("rm %s.png"%outfile)
        plt.axis("off")
        plt.title("(%s) %s (%d)"%(
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[plt_count-1],infile,Lch),
            fontsize=fontsize)
        if plt_count==ncol*nrow:
            plt_count=0
            outfile_list.append("%s.%d.png"%(outfile,gallery_count))
            print(outfile_list[-1])
            plt.tight_layout(pad=0,h_pad=0,w_pad=0)
            plt.savefig(outfile_list[-1],dpi=300)
            plt.close()
            gallery_count+=1
            plt.figure(gallery_count,figsize=(8.27,11.69))
    if plt_count:
        outfile_list.append("%s.%d.png"%(outfile,gallery_count))
        plt.tight_layout(pad=0,h_pad=0,w_pad=0)
        plt.savefig(outfile_list[-1],dpi=300)
        print(outfile_list[-1])
    return outfile_list

if __name__=="__main__":
    if len(sys.argv)!=3:
        sys.stderr.write(docstring)
        exit()

    fp=open(sys.argv[1],'rU')
    infile_list=fp.read().splitlines()
    fp.close()

    outfile=sys.argv[2]

    outfile_list=batch_pdb_image(infile_list,outfile)
