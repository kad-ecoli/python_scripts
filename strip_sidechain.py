#!/usr/bin/env python
docstring='''
strip_sidechain.py fullatom.pdb > backbone.pdb
    strip side chain from full atom PDB "fullatom.pdb"

options:
    -convertAA={GLY,ALA} convert any residue to one specific amino acid type
        default is do not convert
    -amide={true,false} whether preserve amide bond atomsbackbone atoms
'''
import sys
def strip_sidechain(fullatom,convertAA=False,amide=True):
    '''strip side chain from full atom PDB file "fullatom"
    convertAA - convert all amino acid to one specific type
    amide - whether preserve amide bond atoms'''
    atom_name_set=["CA","C ","N ","O "] if amide else ["CA"]
    fp=open(fullatom,'rU')
    fullatom_lines=fp.readlines()
    fp.close()

    backbone_txt=''
    for line in fullatom_lines:
        if line.startswith("ATOM") and line[13:15] in atom_name_set:
            backbone_txt+=line if not convertAA else \
                line[:17]+convertAA+line[20:]
    return backbone_txt

if __name__=="__main__":
    convertAA=False
    amide=True
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-convertAA="):
            convertAA=arg[len("-convertAA="):]
        elif arg.startswith("-amide="):
            amide=(arg[len("-amide="):].lower()=="true")
        elif arg.startswith('-'):            
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)<1:
        sys.stderr.write(docstring)
        exit()

    for arg in argv:
        sys.stdout.write(strip_sidechain(arg,convertAA,amide))
