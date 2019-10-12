#!/usr/bin/env python
docstring='''
create_UNIPROT_GOterms.py  [option] gene_association.*.gz
    create UNIPROT GO terms mapping from list of gene association files
    or SIFTS PDB-GO mapping file
    Output mapping file:

    UNIPROT_GOterms.ALL (any evidence code): IEA EXP IDA IPI IMP 
        IGI IEP ISS ISO ISA ISM IGC IBA IBD IKR IRD RCA TAS NAS IC
    UNIPROT_GOterms.NONIEA (reviewed evidence code): EXP IDA IPI IMP 
        IGI IEP ISS ISO ISA ISM IGC IBA IBD IKR IRD RCA TAS NAS IC
    UNIPROT_GOterms.CAFA (evidence code used in CAFA assessment)
        EXP IDA IMP IGI IEP TAS IC
    UNIPROT_GOterms.{ALL,NONIEA,CAFA}.{F,P,C} (mapping for one Aspect)
    UNIPROT_GOterms.{ALL,NONIEA,CAFA}.is_a (trace back all parent GO terms)
    UNIPROT_GOterms.{ALL,NONIEA,CAFA}.is_a.{F,P,C}

    If the database of entry is PDB, the entry name (DB_Object_ID) will
    be converted from 101M_A to 101mA.

    Evidence code "ND" (No Biological Data Available) are excluded

option:
-excludeGO=GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575
    GO terms to be excluded
    (default) remove root of 3 aspects and "protein binding"
        GO:0005515 ! protein binding
        GO:0005488 ! binding
        GO:0003674 ! molecular_function
        GO:0008150 ! biological_process
        GO:0005575 ! cellular_component

-DB=PDB,UniProtKB,RNAcentral,IntAct 
    comma seperated string of database of entry to use
    (default) map entry for all databases
    
-ID=list 
     A file listing protein for which GO terms are mapped.
     (default) all proteins in the association files
     If a PDB chain is in the list but not in the gene association file, its
     latest chain will be mapped from "large_split_mapping.tsv" from PDB

-pdb2uniprot=pdb2uniprot.map (valid if -ID is set and the input is PDB chain)
    map GO association from uniprot instead of from PDB for all specified
    PDB chain using pdb-uniprot mapping file "pdb2uniprot.map".
    pdb2uniprot.map is a two column tabular file. First column is PDB chain
    and second column is comma seperated list of uniprot ID

-category={CAFA,NONIEA,ALL} evidence code category to use:
    CAFA  : only UNIPROT_GOterms.CAFA* are generated
    NONIEA: only UNIPROT_GOterms.{CAFA,NONIEA}* generated
    ALL   : all UNIPROT_GOterms

-infmt={GOA,SIFTS,INTERPRO} input format of gene association file
    (default): automatically detect whether format is GOA or SIFTS
    GOA      : Gene Assoication file in GAF format
    SIFTS    : PDB GO mapping by SIFTS project
    INTERPRO : TSV format InterProScan output
'''
import sys,os
import urllib
import zlib
import re
import obo2csv # parsing GO hierachy
import gzip
from fetch import large_split_mapping,wget

CAFA_set=("EXP","IDA","IMP","IGI","IEP","TAS","IC")
obo_url="http://geneontology.org/ontology/go-basic.obo"

def parse_GOA(GOA='',
    excludeGO='GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575', 
    DB='',ID_map_dict=dict(),obo_dict=dict(),category="ALL",infmt=''):
    '''parse gene associatation file or SIFTS PDB-GO mapping file "GOA". 
    Return three dict whose key is an Aspect of GO and value is dict itself, 
    whose key is DB_Object_ID and value is GO ID: 
        ALL_dict    - any evidence code
        NONIEA_dict - reviewed evidence code
        CAFA_dict   - experimental evidence code used in CAFA

    excludeGO - comma seperated strings of GO terms to be excluded
    DB        - comma seperated strings of database of entry to use
                default is use all database
    ID_map_dict - dict for mapping old PDB ID to new PDB ID
                key is DB_Object_ID in the gene association file
                value is comma separated list of entry name is the output
                UNIPROT GO mapping file
                default take all entries without mapping
    obo_dict - obo2csv.obo class for parsing SIFTS PDB-GO mapping
    infmt    - input file format. empty for auto-detection between GOA and SIFTS
               GOA      : Gene Assoication file in GAF format
               SIFTS    : PDB GO mapping by SIFTS project
               INTERPRO : TSV format InterProScan output
    '''
    if isinstance(excludeGO,str):
        excludeGO=excludeGO.split(',')
    if DB and isinstance(DB,str):
        DB=DB.split(',')
    
    excludeGO_count=0 # number of uninformative GO terms excluded
    GO_count=0        # number of GO terms in the fetched database

    ALL_dict={'F':dict(),'P':dict(),'C':dict()} # any evidence code
    NONIEA_dict={'F':dict(),'P':dict(),'C':dict()} # reviewed evidence code
    CAFA_dict={'F':dict(),'P':dict(),'C':dict()} # experimental evidence code used in CAFA

    fp=gzip.open(GOA,'rU') if GOA.endswith(".gz") else open(GOA,'rU')
    GOA_lines=fp.read().splitlines()
    fp.close()
    GO_count=len(GOA_lines)

    if not infmt:
        infmt="GOA" # gene association file
        if GOA_lines[0].startswith('#') and GOA_lines[1].startswith("PDB"):
            infmt="SIFTS" # SIFTS PDB-GO mapping file
            GOA_lines=GOA_lines[2:]

    ### parse InterproScan TSV output ###
    if infmt=="INTERPRO":
        ## parse text
        for line in GOA_lines:
            line=line.split('\t')
            accession=line[0]
            if not "GO:" in line[-1]:
                continue # only parse line with GO annotation
            GOterms_list=[e for e in line[-1].split('|') if not e in excludeGO]
            for GOterm in GOterms_list:
                for Aspect in ALL_dict:
                    if GOterm in obo_dict[Aspect]['Term']:
                        if not accession in ALL_dict[Aspect]:
                            ALL_dict[Aspect][accession]=[GOterm]
                        elif not GOterm in ALL_dict[Aspect][accession]:
                            ALL_dict[Aspect][accession].append(GOterm)
        ## sort GO terms for each accession
        for Aspect in ALL_dict:
            for accession in ALL_dict[Aspect]:
                ALL_dict[Aspect][accession].sort()
        return ALL_dict,NONIEA_dict,CAFA_dict

    ### parse GAF or SIFTS format GO association ###
    for line_idx,line in enumerate(GOA_lines):
        if int(GO_count/100) and line_idx % int(GO_count/100) ==0:
            sys.stderr.write('*')
        if line.startswith('!'): # comments
            GO_count-=1
            continue
        line=line.split('\t')

        if infmt=="SIFTS":
            source_database="PDB"
            accession=line[0]+line[1]
            GOterm=line[5]
            evidence=line[4]
            Aspect=obo2csv.GO_namespace_to_Aspect[
                obo_dict.Term(GOterm).namespace]
        else:
            source_database=line[0] # DB
            accession=line[1] # DB_Object_ID
            GOterm=line[4]    # GO ID
            evidence=line[6]  # Evidence Code
            Aspect=line[8]    # namespace

            if source_database=="PDB": # convert PDB chain ID
                accession=accession.split('_')
                accession=accession[0].lower()+accession[-1]
        
        if DB and not source_database in DB: 
            GO_count-=1
            continue # only map entries from selected database
        
        if GOterm in excludeGO or evidence=="ND":
            excludeGO_count+=1
            continue # only map GO terms not excluded

        if ID_map_dict:
            if not accession in ID_map_dict:
                GO_count-=1
                continue # only map entries specified in ID_map_dict
            else:
                accession_list=ID_map_dict[accession].split(',')
        else:
            accession_list=[accession]

        for accession in accession_list:
            if category in ["ALL"]:
                if not accession in ALL_dict[Aspect]:
                    ALL_dict[Aspect][accession]=[GOterm]
                else:
                    ALL_dict[Aspect][accession].append(GOterm)
        
            if evidence!="IEA" and category in ["ALL","NONIEA"]:
                if not accession in NONIEA_dict[Aspect]:
                    NONIEA_dict[Aspect][accession]=[GOterm]
                else:
                    NONIEA_dict[Aspect][accession].append(GOterm)
            
            if evidence in CAFA_set:
                if not accession in CAFA_dict[Aspect]:
                    CAFA_dict[Aspect][accession]=[GOterm]
                else:
                    CAFA_dict[Aspect][accession].append(GOterm)

    ## sort GO terms for each accession
    for Aspect in ALL_dict:
        for accession in ALL_dict[Aspect]:
            ALL_dict[Aspect][accession].sort()
        for accession in NONIEA_dict[Aspect]:
            NONIEA_dict[Aspect][accession].sort()
        for accession in CAFA_dict[Aspect]:
            CAFA_dict[Aspect][accession].sort()

    ## caculate number of accessions
    uniprot_list_NONIEA=NONIEA_dict['F'].keys()+NONIEA_dict['P'].keys() \
                       +NONIEA_dict['C'].keys()
    uniprot_list_CAFA  =CAFA_dict['F'].keys()  +CAFA_dict['P'].keys() \
                       +CAFA_dict['C'].keys()
    sys.stderr.write("\n%u GO terms for %u accession, %u GO terms excluded\n"
      %(GO_count,len(get_accession_list(ALL_dict)), excludeGO_count))
    sys.stderr.write("%u accession has non-IEA GO terms, %u has CAFA GO terms\n"
      %(len(get_accession_list(NONIEA_dict)), 
        len(get_accession_list(CAFA_dict))))
    return ALL_dict,NONIEA_dict,CAFA_dict

def get_accession_list(GOA_dict=dict()):
    '''get a list of all accessions in GOA_dict'''
    uniprot_list=[]
    for Aspect in GOA_dict:
        uniprot_list+=GOA_dict[Aspect].keys()
    return sorted(set(uniprot_list))

def merge_GOA_dict(GOA_dict1=dict(),GOA_dict2=dict()):
    '''merge two dict returned by parse_GOA'''
    if not GOA_dict1:
        return GOA_dict2

    GOA_dict=dict()
    for Aspect in set(GOA_dict1.keys()+GOA_dict2.keys()):
        GOA_dict[Aspect]=dict()
        for DB_Object_ID in set(GOA_dict1[Aspect].keys()+GOA_dict2[Aspect].keys()):
            GOA_dict[Aspect][DB_Object_ID]=sorted(set(
(GOA_dict1[Aspect][DB_Object_ID] if DB_Object_ID in GOA_dict1[Aspect] else [])+ \
(GOA_dict2[Aspect][DB_Object_ID] if DB_Object_ID in GOA_dict2[Aspect] else [])))
        sys.stderr.write("%u accession for Aspect %s\n"
            %(len(GOA_dict[Aspect]),Aspect))
    return GOA_dict

def create_UNIPROT_GOterms(uniprot_list=[],GOA_dict=dict(),
    obo_dict=dict(),excludeGO=''):
    '''create UNIPROT-GO mapping text
    uniprot_list - a list of uniprot accession
    GOA_dict     - uniprot GO mapping generated by parse_GOA
    obo_dict     - dict returned by obo2csv.parse_obo_txt to trace back 
                   parent GO terms. If empty, do not trace back
    excludeGO    - comma seperated strings of GO terms to be excluded
    Aspect       - GO Aspects to map. Default is to map all 3 Aspect
    '''
    if isinstance(excludeGO,str):
        excludeGO=excludeGO.split(',')
    excludeGO=set(excludeGO)
    Aspect_list='FPC'

    UNIPROT_GOterms=[]
    UNIPROT_GOterms_Aspect={'F':[],'P':[],'C':[]}
    for p in uniprot_list:
        GOterms_list=[]
        GOterms_Aspect={'F':[],'P':[],'C':[]}
        for Aspect in Aspect_list:
            if not p in GOA_dict[Aspect]:
                continue
            GOterms_Aspect[Aspect]=GOA_dict[Aspect][p]
            if obo_dict: # trace back all parent nodes
                for Term_id in set(GOterms_Aspect[Aspect]):
                    if not Term_id in obo_dict[Aspect]["Term"]:
                        continue
                    GOterms_Aspect[Aspect]+=obo_dict.is_a(Term_id,
                        direct=False,name=False,number=False).split()
            GOterms_Aspect[Aspect]=sorted(set(GOterms_Aspect[Aspect])-excludeGO)
            if GOterms_Aspect[Aspect]:
                UNIPROT_GOterms_Aspect[Aspect].append(
                    p+'    '+','.join(GOterms_Aspect[Aspect]))
                GOterms_list+=GOterms_Aspect[Aspect]

        UNIPROT_GOterms.append(p+'    '+','.join(sorted(GOterms_list)))

    GOterms_txt='\n'.join(UNIPROT_GOterms)+'\n'
    GOterms_txt_Aspect={'F':[],'P':[],'C':[]}
    for Aspect in Aspect_list:
        GOterms_txt_Aspect[Aspect]='\n'.join(
            UNIPROT_GOterms_Aspect[Aspect])+'\n'
    return GOterms_txt,GOterms_txt_Aspect

def write_UNIPROT_GOterms_txt(GOterms_txt,GOterms_txt_Aspect, filename):
    '''write output mapping file'''
    sys.stdout.write(filename+'\n')
    fp=open(filename,'w')
    fp.write(GOterms_txt)
    fp.close()
    for Aspect in ['MF','BP','CC']:
        filenameAspect=filename+'.'+Aspect
        sys.stdout.write(filenameAspect+'\n')
        fp=open(filenameAspect,'w')
        fp.write(GOterms_txt_Aspect[Aspect[-1]])
        fp.close()
    return

def map_pdb_ID_list(ID_list=[],pdb2uniprot_list=[]):
    '''map list of obsolete PDB chain to supersede PDB chain
    return a dict whose key is old ID and value is new ID

    pdb2uniprot_list - a list containing mapping from PDB chain to 
        uniprot accession
    '''
    ID_map_dict=dict()
    if not ID_list:
        return ID_map_dict
    large_split_dict=large_split_mapping()
    for old_chain in ID_list:
        if old_chain in large_split_dict.old2new:
            new_chain=large_split_dict.old2new[old_chain]
        else:
            new_chain=old_chain
        ID_map_dict[new_chain]=old_chain

    if not pdb2uniprot_list:
        return ID_map_dict
        
    pdb2uniprot_map_dict=dict()
    old_chain_list=ID_map_dict.values()
    for line in pdb2uniprot_list:
        pdb_chain,uniprot_str=line.split()
        if pdb_chain in ID_map_dict:
            pdb_chain=ID_map_dict[pdb_chain]
        elif not pdb_chain in old_chain_list:
            continue # pdb_chain is not in ID_list
        for uniprot in uniprot_str.split(','):
            if uniprot in pdb2uniprot_map_dict:
                pdb2uniprot_map_dict[uniprot]+=','+pdb_chain
            else:
                pdb2uniprot_map_dict[uniprot]=pdb_chain
    return pdb2uniprot_map_dict

if __name__=="__main__":
    # Default parameters for options
    excludeGO="GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575" # GO to be excluded
    DB=''
    ID=''
    pdb2uniprot=''
    category='ALL'
    infmt='' # guess whether it is GOA or SIFTS

    # parse arguments
    argv=[] # input FASTA format alignment files
    for arg in sys.argv[1:]:
        if not arg.startswith('-'):
            argv.append(arg)
            continue
        elif arg.startswith('-excludeGO='):
            excludeGO=arg[len("-excludeGO="):].upper()
        elif arg.startswith('-DB='):
            DB=arg[len("-DB="):]
        elif arg.startswith('-ID='):
            ID=arg[len("-ID="):]
        elif arg.startswith('-pdb2uniprot='):
            pdb2uniprot=arg[len("-pdb2uniprot="):]
        elif arg.startswith('-category='):
            category=arg[len("-category="):].upper()
        elif arg.startswith('-infmt='):
            infmt=arg[len("-infmt="):].upper()
        else:
            sys.stderr.write("ERROR! Unknown argument "+arg+'\n')
            exit()

    if pdb2uniprot and not ID:
        sys.stderr.write("ERROR! -pdb2uniprot must be used with -ID\n")
        exit()
    if not len(argv):
        sys.stderr.write(docstring)
        exit()
    
    pdb2uniprot_list=[]
    if pdb2uniprot:
        fp=open(pdb2uniprot,'rU')
        pdb2uniprot_list=[line for line in fp.read().splitlines() if line.strip()]
        fp.close()

    ID_list=[] # list of proteins for which GO terms are mapped
    if ID:
        fp=open(ID,'rU')
        ID_list=[line for line in fp.read().splitlines() if line.strip()]
        fp.close()
    ID_map_dict=map_pdb_ID_list(ID_list,pdb2uniprot_list)

    ## parse GO hierachy
    fp=open(wget(obo_url,show_url=True),'rU')
    obo_txt=fp.read()
    fp.close()
    obo_dict=obo2csv.parse_obo_txt(obo_txt)
    
    ## extract all uniprot-GO mapping
    merge_ALL_dict=dict()
    merge_NONIEA_dict=dict()
    merge_CAFA_dict=dict()
    for GOA in argv:
        ALL_dict,NONIEA_dict,CAFA_dict=parse_GOA(GOA,excludeGO,DB, \
            ID_map_dict,obo_dict,category,infmt)
        if category in ["ALL"]:
            merge_ALL_dict=merge_GOA_dict(merge_ALL_dict,ALL_dict)
        if category in ["NONIEA","ALL"]:
            merge_NONIEA_dict=merge_GOA_dict(merge_NONIEA_dict,NONIEA_dict)
        merge_CAFA_dict=merge_GOA_dict(merge_CAFA_dict,CAFA_dict)

    ## list of uniprot accession annotated with GO terms
    uniprot_list_ALL=get_accession_list(merge_ALL_dict)
    uniprot_list_NONIEA=get_accession_list(merge_NONIEA_dict)
    uniprot_list_CAFA=get_accession_list(merge_CAFA_dict)

    ## write output mapping files
    for is_a in [False,True]: # whether trace back parent node
        suffix=".is_a"*is_a

        # all GO terms
        if category in ["ALL"]:
            GOterms_txt,GOterms_txt_Aspect=create_UNIPROT_GOterms(
                uniprot_list_ALL,merge_ALL_dict,
                obo_dict if is_a else False,excludeGO)
            if infmt=="INTERPRO":
                filename="INTERPRO_GOterms"+suffix
            else:
                filename="UNIPROT_GOterms.ALL"+suffix
            write_UNIPROT_GOterms_txt(
                GOterms_txt, GOterms_txt_Aspect, filename)

        # reviewed GO terms
        if category in ["NONIEA","ALL"] and infmt!="INTERPRO":
            GOterms_txt,GOterms_txt_Aspect=create_UNIPROT_GOterms(
                uniprot_list_NONIEA,merge_NONIEA_dict,
                obo_dict if is_a else False,excludeGO)
            filename="UNIPROT_GOterms.NONIEA"+suffix
            write_UNIPROT_GOterms_txt(
                GOterms_txt, GOterms_txt_Aspect, filename)

        # experimental GO terms used in CAFA evalulation
        if infmt!="INTERPRO":
            GOterms_txt,GOterms_txt_Aspect=create_UNIPROT_GOterms(
                uniprot_list_CAFA,merge_CAFA_dict,
                obo_dict if is_a else False,excludeGO)
            filename="UNIPROT_GOterms.CAFA"+suffix
            write_UNIPROT_GOterms_txt(
                GOterms_txt, GOterms_txt_Aspect, filename)
