#!/usr/bin/python
docstring='''
addchainID.py input.pdb output.pdb chainID
'''
import sys

if len(sys.argv)!=4:
    sys.stderr.write(docstring)
    exit()

infile =sys.argv[1]
outfile=sys.argv[2]
chainID=sys.argv[3]

txt=''
fp=open(infile)
for line in fp.read().splitlines():
    if line.startswith("ATOM  ") or line.startswith("HETATM"):
        txt+=line[:21]+chainID+line[22:]+'\n'
fp.close()
fp=open(outfile,'w')
fp.write(txt)
fp.close()
