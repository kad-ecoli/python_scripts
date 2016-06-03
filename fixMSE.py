#!/usr/bin/env python
# 2016-03-28 Chengxin Zhang
docstring='''
fixMSE.py infile.pdb outfile.pdb
    convert MSE or other non-standard residues to ATOM record in PDB
    By default, only MSE are converted

fixMSE.py -clean=true infile.pdb outfile.pdb
    In additional to MSE conversion, only preserve "ATOM" for 20 standard 
    amino acid and "TER".
'''
import sys,os

code_MSE= {'MSE':'M'}


code_with_modified_residues = {
    'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M', 'ILE':'I',
    'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K', 'ARG':'R', 'SER':'S',
    'THR':'T', 'TYR':'Y', 'HIS':'H', 'CYS':'C', 'ASN':'N', 'GLN':'Q',
    'TRP':'W', 'GLY':'G',                                  'MSE':'M',
    '2AS':'D', '3AH':'H', '5HP':'E', 'ACL':'R', 'AIB':'A', 'ALM':'A', 
    'ALO':'T', 'ALY':'K', 'ARM':'R', 'ASA':'D', 'ASB':'D', 'ASK':'D', 
    'ASL':'D', 'ASQ':'D', 'AYA':'A', 'BCS':'C', 'BHD':'D', 'BMT':'T', 
    'BNN':'A', 'BUC':'C', 'BUG':'L', 'C5C':'C', 'C6C':'C', 'CCS':'C', 
    'CEA':'C', 'CHG':'A', 'CLE':'L', 'CME':'C', 'CSD':'A', 'CSO':'C', 
    'CSP':'C', 'CSS':'C', 'CSW':'C', 'CXM':'M', 'CY1':'C', 'CY3':'C', 
    'CYG':'C', 'CYM':'C', 'CYQ':'C', 'DAH':'F', 'DAL':'A', 'DAR':'R', 
    'DAS':'D', 'DCY':'C', 'DGL':'E', 'DGN':'Q', 'DHA':'A', 'DHI':'H', 
    'DIL':'I', 'DIV':'V', 'DLE':'L', 'DLY':'K', 'DNP':'A', 'DPN':'F', 
    'DPR':'P', 'DSN':'S', 'DSP':'D', 'DTH':'T', 'DTR':'W', 'DTY':'Y', 
    'DVA':'V', 'EFC':'C', 'FLA':'A', 'FME':'M', 'GGL':'E', 'GLZ':'G', 
    'GMA':'E', 'GSC':'G', 'HAC':'A', 'HAR':'R', 'HIC':'H', 'HIP':'H', 
    'HMR':'R', 'HPQ':'F', 'HSD':'H', 'HSE':'H', 'HSP':'H', 'HTR':'W', 
    'HYP':'P', 'IIL':'I', 'IYR':'Y', 'KCX':'K', 'LLY':'K', 'LTR':'W',
    'LYM':'K', 'LYZ':'K', 'MAA':'A', 'MEN':'N', 'MHS':'H', 'MIS':'S',
    'MLE':'L', 'MPQ':'G', 'MSA':'G', 'MVA':'V', 'NEM':'H', 'NEP':'H', 
    'NLE':'L', 'NLN':'L', 'NLP':'L', 'NMC':'G', 'OAS':'S', 'OCS':'C', 
    'OMT':'M', 'PAQ':'Y', 'PCA':'E', 'PEC':'C', 'PHI':'F', 'PHL':'F', 
    'PR3':'C', 'PRR':'A', 'PTR':'Y', 'SAC':'S', 'SAR':'G', 'SCH':'C', 
    'SCS':'C', 'SCY':'C', 'SEL':'S', 'SEP':'S', 'SET':'S', 'SHC':'C', 
    'SHR':'K', 'SOC':'C', 'STY':'Y', 'SVA':'S', 'TIH':'A', 'TPL':'W', 
    'TPO':'T', 'TPQ':'A', 'TRG':'K', 'TRO':'W', 'TYB':'Y', 'TYQ':'Y', 
    'TYS':'Y', 'TYY':'Y', 'AGM':'R', 'GL3':'G', 'SMC':'C', 'CGU':'E',
    'CSX':'C',
    'SEC':'U', 'PYL':'O', 'ASX':'B', 'GLX':'Z', 'LLP':'X', 'UNK':'X',
    }

aa1to3 = { # 3 letter to 1 letter amino acid code conversion
    'A':'ALA', 'V':'VAL', 'F':'PHE', 'P':'PRO', 'M':'MET',
    'I':'ILE', 'L':'LEU', 'D':'ASP', 'E':'GLU', 'K':'LYS',
    'R':'ARG', 'S':'SER', 'T':'THR', 'Y':'TYR', 'H':'HIS',
    'C':'CYS', 'N':'ASN', 'Q':'GLN', 'W':'TRP', 'G':'GLY',
    'B':'ASX', 'Z':'GLX', 'U':'SEC', 'O':'PYL', 'X':'UNK',
    'J':'UNK', 
    }

standard_amino_acid={
    'ALA', 'VAL', 'PHE', 'PRO', 'MET', 'ILE', 'LEU', 'ASP', 'GLU', 'LYS',
    'ARG', 'SER', 'THR', 'TYR', 'HIS', 'CYS', 'ASN', 'GLN', 'TRP', 'GLY'}

def fixMSE(infile="pdb.pdb",aa3to1=code_MSE,clean=False):
    '''fix MSE and other HETATM record in PDB file "infile"
    return the fixed text
    clean - whether only preserve "ATOM" for 20 standard amino acid and "TER".
    '''
    fp=open(infile,'rU')
    txt=fp.read()
    fp.close()
    fix_txt=fixMSE_txt(txt=txt,aa3to1=aa3to1,clean=clean)
    return fix_txt 

def fixMSE_txt(txt='',aa3to1=code_MSE,clean=False):
    '''fix MSE and other HETATM record in PDB text "txt"
    return the fixed text

    clean - whether only preserve "ATOM" for 20 standard amino acid and "TER".
    '''
    fix_txt=[]
    for line in txt.splitlines():
        if not line.strip():
            continue
        if not line.startswith("ATOM  ") and not line.startswith("HETATM"):
            if clean==False or line.startswith("TER"):
                fix_txt.append(line)
                continue
            else:
                if line.startswith("ENDMDL"): # only use return model 1
                    break
            continue
        resn=line[17:20]
        if resn in aa3to1:
            line="ATOM  "+line[6:17]+ \
                aa1to3[aa3to1[resn]]+line[20:]
        if resn=="MSE" and line[12:14]=="SE": # convert SE to S
            line=line[:12]+" SD"+line[15:]
            if len(line)>=78 and line[76:78]=="SE":
                line=line[:76]+' S'+line[78:]
        if clean==False or line[17:20] in standard_amino_acid:
            fix_txt.append(line)
    return '\n'.join(fix_txt)+'\n'

if __name__=="__main__":
    clean=False
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-clean="):
            clean=(arg[len("-clean="):].lower()=="true")
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)<2:
        print >>sys.stderr,docstring
        exit()

    fix_txt=fixMSE(argv[0],clean=clean)
    if len(argv)>1:
        fp=open(argv[1],'w')
        fp.write(fix_txt)
        fp.close()
    else:
        sys.stdout.write(fix_txt)
