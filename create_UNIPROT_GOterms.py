#!/usr/bin/env python
docstring='''
create_UNIPROT_GOterms.py  [option] gene_association.*.gz
    create UNIPROT GO terms mapping from gene association file list
    Output mapping file:

    UNIPROT_GOterms.ALL (any evidence code):
        EXP IDA IPI IMP IGI IEP ISS ISO ISA ISM IGC IBA IBD IKR IRD RCA TAS NAS IC IEA
    UNIPROT_GOterms.NONIEA (reviewed evidence code):
        EXP IDA IPI IMP IGI IEP ISS ISO ISA ISM IGC IBA IBD IKR IRD RCA TAS NAS IC
    UNIPROT_GOterms.CAFA (experimental evidence code used in CAFA assessment)
        EXP IDA IMP IGI IEP TAS IC
    UNIPROT_GOterms.{ALL,NONIEA,CAFA}.is_a (trace back all parent GO terms)

option:
-excludeGO=GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575
    (default) remove root of 3 aspects and "protein-binding"
        GO:0005515 ! protein binding
        GO:0005488 ! binding
        GO:0003674 ! molecular_function
        GO:0008150 ! biological_process
        GO:0005575 ! cellular_component
'''
import sys,os
#import urllib2
import urllib
import zlib
import re
import obo2csv # parsing GO hierachy
import gzip

obo_url="http://geneontology.org/ontology/go-basic.obo"

CAFA_set=("EXP","IDA","IMP","IGI","IEP","TAS","IC")

def wget(url=''):
    '''retrieve file from internet if not exists at current directory.
    return always basename'''
    basename=os.path.basename(url)
    if not os.path.isfile(basename):
        sys.stderr.write("fetching %s\n"%url)
        urllib.urlretrieve(url, basename)
    else:
        sys.stderr.write("%s already exists\n"%basename)
    return basename

def parse_GOA(GOA='',
    excludeGO='GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575'):
    '''parse gene associatation file "GOA". 
    Return three dict whose key is DB_Object_ID and value is GO ID: 
        ALL_dict    - any evidence code
        NONIEA_dict - reviewed evidence code
        CAFA_dict   - experimental evidence code used in CAFA

    excludeGO - comma seperated strings of GO terms to be excluded
    '''
    if isinstance(excludeGO,str):
        excludeGO=excludeGO.split(',')
    
    excludeGO_count=0 # number of uninformative GO terms excluded
    GO_count=0        # number of GO terms in the fetched database

    ALL_dict=dict()   # any evidence code
    NONIEA_dict=dict()# reviewed evidence code
    CAFA_dict=dict()  # experimental evidence code used in CAFA

    fp=gzip.open(GOA,'rU') if GOA.endswith(".gz") else open(GOA,'rU')
    GOA_lines=fp.read().splitlines()
    fp.close()
    GO_count=len(GOA_lines)
    for line_idx,line in enumerate(GOA_lines):
        if int(GO_count/100) and line_idx % int(GO_count/100) ==0:
            sys.stderr.write('*')
        if line.startswith('!'): # comments
            GO_count-=1
            continue
        line=line.split('\t')
        if line[0] in {"IntAct","RNAcentral"}: 
            GO_count-=1
            continue # do not deal with PPI complex and RNA for now

        accession=line[1] # DB_Object_ID
        GOterm=line[4]    # GO ID
        evidence=line[6]  # Evidence Code

        if GOterm in excludeGO or evidence=="ND":
            excludeGO_count+=1
            continue
        
        if line[0]=="PDB":
            accession=accession.split('_')
            accession=accession[0].lower()+accession[-1]
        
        if not accession in ALL_dict:
            ALL_dict[accession]=[GOterm]
        else:
            ALL_dict[accession].append(GOterm)
        
        if evidence!="IEA":
            if not accession in NONIEA_dict:
                NONIEA_dict[accession]=[GOterm]
            else:
                NONIEA_dict[accession].append(GOterm)
            
        if evidence in CAFA_set:
            if not accession in CAFA_dict:
                CAFA_dict[accession]=[GOterm]
            else:
                CAFA_dict[accession].append(GOterm)

    for accession in ALL_dict:
        ALL_dict[accession]=','.join(sorted(set(ALL_dict[accession])))
    for accession in NONIEA_dict:
        NONIEA_dict[accession]=','.join(sorted(set(NONIEA_dict[accession])))
    for accession in CAFA_dict:
        CAFA_dict[accession]=','.join(sorted(set(CAFA_dict[accession])))
    sys.stderr.write("\n%u GO terms for %u accession, %u GO terms excluded\n"
        %(GO_count,len(ALL_dict),excludeGO_count))
    sys.stderr.write("%u accession has non-IEA GO terms, %u has CAFA GO terms\n"
        %(len(NONIEA_dict),len(CAFA_dict)))
    return ALL_dict,NONIEA_dict,CAFA_dict

def merge_GOA_dict(GOA_dict1=dict(),GOA_dict2=dict()):
    '''merge two dict returned by parse_GOA'''
    if not GOA_dict1:
        return GOA_dict2
    GOA_dict=dict()
    for DB_Object_ID in set(GOA_dict1.keys()+GOA_dict2.keys()):
        GOA_dict[DB_Object_ID]=','.join(sorted(set(
        (GOA_dict1[DB_Object_ID].split(',') if DB_Object_ID in GOA_dict1 else [])+(
         GOA_dict2[DB_Object_ID].split(',') if DB_Object_ID in GOA_dict2 else []))))
    sys.stderr.write("%u accession\n"%len(GOA_dict))
    return GOA_dict

def create_UNIPROT_GOterms(uniprot_list=[],GOA_dict=dict(),obo_dict=dict(),excludeGO=''):
    '''create UNIPROT-GO mapping text
    uniprot_list - a list of uniprot accession
    GOA_dict     - uniprot GO mapping generated by parse_GOA
    obo_dict     - dict returned by obo2csv.parse_obo_txt to trace back 
                   parent GO terms. If empty, do not trace back
    excludeGO    - comma seperated strings of GO terms to be excluded
    '''
    if isinstance(excludeGO,str):
        excludeGO=excludeGO.split(',')
    excludeGO=set(excludeGO)

    if not obo_dict:
        UNIPROT_GOterms=[p+'    '+GOA_dict[p] for p in uniprot_list if p in GOA_dict]
    else:
        #is_a(self,Term_id='',direct=True,name=False,number=False)
        UNIPROT_GOterms=[]
        for p in uniprot_list:
            if p in GOA_dict:
                GOterms_list=GOA_dict[p].split(',')
                for Term_id in set(GOterms_list):
                    GOterms_list+=obo_dict.is_a(Term_id,
                        direct=False,name=False,number=False).split()
                GOterms_set=set(GOterms_list)-excludeGO
                UNIPROT_GOterms.append(p+'    '+','.join(sorted(GOterms_set)))
    UNIPROT_GOterms_txt='\n'.join(UNIPROT_GOterms)+'\n'
    return UNIPROT_GOterms_txt

if __name__=="__main__":
    # Default parameters for options
    excludeGO="GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575" # GO to be excluded

    # parse arguments
    argv=[] # input FASTA format alignment files
    for arg in sys.argv[1:]:
        if not arg.startswith('-'):
            argv.append(arg)
            continue
        elif arg.startswith('-excludeGO='):
            excludeGO=arg[len("-excludeGO="):].upper()
        else:
            print >>sys.stderr, "ERROR! Unknown argument "+arg
            exit()

    if not len(argv):
        sys.stderr.write(docstring)
        exit()
    
    ## extract all uniprot-GO mapping
    merge_ALL_dict=dict()
    merge_NONIEA_dict=dict()
    merge_CAFA_dict=dict()
    for GOA in argv:
        ALL_dict,NONIEA_dict,CAFA_dict=parse_GOA(GOA,excludeGO)
        merge_ALL_dict=merge_GOA_dict(merge_ALL_dict,ALL_dict)
        merge_NONIEA_dict=merge_GOA_dict(merge_NONIEA_dict,NONIEA_dict)
        merge_CAFA_dict=merge_GOA_dict(merge_CAFA_dict,CAFA_dict)

    ## list of uniprot accession annotated with GO terms
    uniprot_list_ALL=ALL_dict.keys()
    uniprot_list_NONIEA=NONIEA_dict.keys()
    uniprot_list_CAFA=CAFA_dict.keys()

    ## not tracing parent GO terms
    print "UNIPROT_GOterms.ALL"
    fp=open("UNIPROT_GOterms.ALL",'w')
    fp.write(create_UNIPROT_GOterms(
        uniprot_list_ALL,merge_ALL_dict,False,excludeGO))
    fp.close()

    print "UNIPROT_GOterms.NONIEA"
    fp=open("UNIPROT_GOterms.NONIEA",'w')
    fp.write(create_UNIPROT_GOterms(
        uniprot_list_NONIEA,merge_NONIEA_dict,False,excludeGO))
    fp.close()

    print "UNIPROT_GOterms.CAFA"
    fp=open("UNIPROT_GOterms.CAFA",'w')
    fp.write(create_UNIPROT_GOterms(
        uniprot_list_CAFA,merge_CAFA_dict,False,excludeGO))
    fp.close()

    ## tracing parent GO terms
    fp=open(wget(obo_url),'rU')
    obo_txt=fp.read()
    fp.close()
    obo_dict=obo2csv.parse_obo_txt(obo_txt)

    print "UNIPROT_GOterms.ALL.is_a"
    fp=open("UNIPROT_GOterms.ALL.is_a",'w')
    fp.write(create_UNIPROT_GOterms(
        uniprot_list_ALL,merge_ALL_dict,obo_dict,excludeGO))
    fp.close()

    print "UNIPROT_GOterms.NONIEA.is_a"
    fp=open("UNIPROT_GOterms.NONIEA.is_a",'w')
    fp.write(create_UNIPROT_GOterms(
        uniprot_list_NONIEA,merge_NONIEA_dict,obo_dict,excludeGO))
    fp.close()

    print "UNIPROT_GOterms.CAFA.is_a"
    fp=open("UNIPROT_GOterms.CAFA.is_a",'w')
    fp.write(create_UNIPROT_GOterms(
        uniprot_list_CAFA,merge_CAFA_dict,obo_dict,excludeGO))
    fp.close()

