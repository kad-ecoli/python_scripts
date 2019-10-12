#!/usr/bin/env python
docstring='''
batch_retrieval.py list uniprot.fasta
    batch retrieve uniprot entries listed in "list" from UniProt website
    output fetch result to "uniprot.fasta"

Options:
    -outfmt={fasta,txt,tab,xml,rdf,list,html,gaf,gpad,tsv}
        output format. default "fasta"
        If output format is "fasta", update output fetch result by adding
        new entries into output fetch result.
        Otherwise, overwrite output fetch result.

        If output format is gaf, gpad or tsv, retrieve the Gene Ontology
        for respective entries. In this case, -infmt must be ACC. 
        -db is ignored.

    -infmt={ACC+ID,ACC,ID,UPARC,NF50,NF90,NF100,GENENAME,...}
        input format. default 'ACC' (uniprot accession). See
        http://www.uniprot.org/help/programmatic_access#id_mapping_examples

    -db=uniprot_sprot.fasta.gz
        (only valid if -outfmt=fasta)
        directly extract fasta sequence from FASTA format sequence database
        only accession as will be preserved as fasta header
'''
import urllib,urllib2
import sys,os
import subprocess

url_mapping = 'http://www.uniprot.org/mapping/'
url_upload = 'http://www.uniprot.org/uploadlists/'
# set email here to help uniprot guys debug in case of problems.
email = "zcx@umich.edu"
# uniprot cannot parse very long list. If given a long list, the full query
# list will be splitted into small lists of "split_size" entries
split_size=10000

from string import Template
import requests

requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?limit=100&geneProductId=UniProtKB:"

def batch_retrival(query_list=[],infmt="ACC",outfmt="fasta",db=''):
    '''
    If "db" is empty, batch retrieve uniprot entries listed in 
    "query_list" from uniprot website.

    If "db" is a file and "outfmt" is "fasta", directly extract fasta
    sequence from "db"
    '''
    if outfmt=="fasta" and db:
        if db.endswith(".gz"):
            import gzip
            fp=gzip.open(db,'rU')
        else:
            fp=open(db,'rU')
        page=''
        read_entry=False # whether to read current entry
        for line in fp: # read line by line because trembl is huge
            if line.startswith('>'):
                read_entry=(line[1:].split()[0] in query_list or ('|' in
                    line and line.split('|')[1].split()[0] in query_list))
            if read_entry:
                page+=line
        fp.close()
        return page

    params = {
        'from':infmt,
        'to':'ACC',
        'format':outfmt, # html tab fasta txt xml rdf list
        'query':' '.join(query_list),
        'compress':'yes',
    }

    data = urllib.urlencode(params)
    #request = urllib2.Request(url_mapping, data)
    request = urllib2.Request(url_upload, data)
    request.add_header('User-Agent', 'Python %s' % email)
    response = urllib2.urlopen(request)
    page = response.read()
    return page

def get_accession_list_from_fasta(fasta_file):
    '''get a list of accessions from fasta. 
    This is faster than a pure python implement.
    '''
    if not os.path.isfile(fasta_file):
        return []
    cmd='z'*fasta_file.endswith(".gz")+ "grep '>' %s|"%fasta_file+ \
        "sed 's/^>[a-z]\{1,\}|//g'|sed 's/^>//g'|cut -d' ' -f1|cut -d'|' -f1"
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    return stdout.splitlines()

def remove_entry_not_in_fasta_file(query_list,fasta_file):
    '''remove entries not in fasta file'''
    fasta_header_list=get_accession_list_from_fasta(fasta_file)
    query_list=list(set(query_list).intersection(fasta_header_list))
    return query_list

def remove_entry_in_fasta_file(query_list,fasta_file):
    '''remove entries in fasta file'''
    fasta_header_list=get_accession_list_from_fasta(fasta_file)
    query_list=list(set(query_list).difference(fasta_header_list))
    return query_list

def request_GO(query,outfmt="gaf"):
    ''' request gpad or gaf format GO annotation '''
    page=''
    txt=''
    p=1
    while p:
        r = requests.get("%s%s&page=%d"%(requestURL,query,p),
            headers={ "Accept" : "text/"+outfmt})
        if not r.ok:
            r.raise_for_status()
        responseBody = ''.join([line+'\n' for line in str(r.text
            ).splitlines() if line.startswith("UniProtKB")])
        if responseBody==txt:
            break
        txt=str(responseBody)
        page+=responseBody
        p+=1
    return page

if __name__=="__main__":
    infmt="ACC"        # input uniprot accessions
    outfmt="fasta"     # output fasta sequences
    db=''              # fasta format sequence database

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-infmt="):
            infmt=arg[len("-infmt="):].upper()
        elif arg.startswith("-outfmt"):
            outfmt=arg[len("-outfmt="):].lower()
        elif arg.startswith("-db="):
            db=arg[len("-db="):]
            if db.startswith('~/'):
                db=os.getenv("HOME")+db[1:]
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)
    
    if len(argv)<2:
        sys.stderr.write(docstring)
        exit()
    if db and not os.path.isfile(db):
        sys.stderr.write("ERROR! No such file %s\n"%db)
        exit()

    fp=open(argv[0],'rU')
    query_list=fp.read().splitlines()
    fp.close()
    sys.stdout.write("%u query entries\n"%len(query_list))

    if outfmt in ["gaf","gpad","tsv"]:
        if infmt!="ACC":
            sys.stderr.write("ERROR! -infmt must be 'ACC' when -outfmt=gaf")
            exit()
        if db:
            sys.stderr.write("WARNING! ignore -db because -outfmt=gaf")
        fp=open(argv[1],'w')
        for query in query_list:
            page=request_GO(query,outfmt)
            fp.write(page)
            fp.flush()
    else:
        if db:
            if infmt=="ACC":
                query_list=remove_entry_not_in_fasta_file(query_list,db)
                sys.stdout.write("%u query entries in db\n"%len(query_list))

        if outfmt=="fasta":
            if infmt=="ACC":
                query_list=remove_entry_in_fasta_file(query_list,argv[1])
                sys.stdout.write("%u query entries not fetched\n"%len(query_list))
            fp=open(argv[1],'a')
        else:
            fp=open(argv[1],'w')
        for i in range(0,len(query_list),split_size):
            sys.stdout.write("Retrieving entries %u-%u\n"%(
                i+1,min(i+split_size,len(query_list))))
            while (True): # keep trying to retrieve entry until success
                try:
                    page=batch_retrival(query_list[i:i+split_size],
                        infmt,outfmt,db)
                    fp.write(page)
                    fp.flush()
                    break
                except Exception,error:
                    sys.stderr.write(str(error)+'\n')
    fp.close()
