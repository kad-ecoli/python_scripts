#!/usr/bin/env python
docstring='''
batch_retrieval.py list > uniprot-yourlist
    batch retrieve uniprot entries listed in "list" from UniProt website
    output fetch result to "uniprot-yourlist"

Options:
    -outfmt={fasta,txt,tab,xml,rdf,list,html}
        output format. default "fasta"

    -infmt={ACC+ID,ACC,ID,UPARC,NF50,NF90,NF100,GENENAME,...}
        input format. default 'ACC' (uniprot accession). See
        http://www.uniprot.org/help/programmatic_access#id_mapping_examples
'''
import urllib,urllib2
import sys

url_mapping = 'http://www.uniprot.org/mapping/'
url_upload = 'http://www.uniprot.org/uploadlists/'
# set email here to help uniprot guys debug in case of problems.
email = "zcx@umich.edu"

def batch_retrival(query_list=[],infmt="ACC",outfmt="fasta"):
    params = {
        'from':infmt,
        'to':'ACC',
        'format':outfmt, # html tab fasta txt xml rdf list
        'query':' '.join(query_list),
        'compress':'yes',
    }

    data = urllib.urlencode(params)
    request = urllib2.Request(url_mapping, data)
    request.add_header('User-Agent', 'Python %s' % email)
    response = urllib2.urlopen(request)
    page = response.read()
    return page

if __name__=="__main__":
    infmt="ACC"        # input uniprot accessions
    outfmt="fasta"     # output fasta sequences
    split_size=1000000 # split full query list into small 
                       # lists of one million entries

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-infmt="):
            infmt=arg[len("-infmt="):].upper()
        elif arg.startswith("-outfmt"):
            outfmt=arg[len("-outfmt="):].lower()
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)
    
    if len(argv)==0:
        sys.stderr.write(docstring)
        exit()

    fp=open(argv[0],'rU')
    query_list=fp.read().splitlines()
    fp.close()

    txt=''
    for i in range(0,len(query_list),split_size):
        page=batch_retrival(query_list[i:i+split_size],infmt,outfmt)
        txt+=page

    if len(argv)>1:
        fp=open(argv[1],'w')
        fp.write(txt)
        fp.close()
    else:
        sys.stdout.write(txt)
