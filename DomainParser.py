#!/usr/bin/env python
docstring='''DomainParser.py 16pkA.pdb > domain_list.txt
    run DomainParser on target PDB 16pkA.pdb to paritition it into domains. 
    Domains will be listed in "domain_list.txt" in the following format:
        filename	length	domain_number	list_of_domain
    (where list of domain is space separated list of PDB residue number ranges quoted inside parathesis)
    or, in case no domain are found:
        filename	length	1

options:
    -execpath=./domainparser2.LINUX
        path to DomainParser executable. By default it is guessed by 
        location of this script
    -dssp_path=./dssp
        path to DSSP executable. By default it is guessed by
        location of this script
    -log=DomainParser.log
        path to DomainParser output log. "-" for stdout
'''
import sys,os
import shutil
import random
import subprocess
import Bio.PDB

def DomainParser(pdb_file,execpath="domainparser2.LINUX",dssp_path="dssp"):
    '''run DomainParser executable "execpath" using DSSP executable
    "dssp_path" on pdb file "pdb_file"
    '''
    #### make temporary folder ####
    basename=os.path.basename(pdb_file)
    tmp_dir="/tmp/"+os.getenv("USER")+"/DomainParser"+ \
        str(random.randint(1000,9999))+basename.split('.')[0]+'/'
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    tmp_pdb=tmp_dir+'xxxx.pdb'

    #### parse PDB files ####
    struct = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure(pdb_file,pdb_file)
    model=struct[0]
    chain=[c for c in model][0]
    io=Bio.PDB.PDBIO()
    io.set_structure(chain)
    io.save(tmp_pdb)

    #### run DomainParser ####
    cmd=' '.join(['cd',tmp_dir,';',
        'export DSSP_PATH='+dssp_path,';',
        execpath,'xxxx'+chain.id
        ])
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()

    #### parse output ####
    output_list=stdout.split()
    target,seqlen,domain_num=output_list[:3]
    domain_list=output_list[3:]

    #### output individual domain ####
    for domain_idx,domain in enumerate(domain_list):
        domain_file=basename.split('.')[0]+str(domain_idx+1)+".pdb"
        resi_list=[]
        for segment in domain.strip('()').split(';'):
            resi_start,resi_end=segment.split('-')
            resi_list+=range(int(resi_start),int(resi_end)+1)
        
        class ResiSelect(Bio.PDB.Select): # class to select domain
            def accept_residue(self,residue):
                return 1 if residue.id[1] in resi_list else 0

        io.save(domain_file,ResiSelect())
        sys.stdout.write(domain_file+'\n')

    #### cleanup temporary folder ####
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    return '\t'.join([basename,seqlen,domain_num,' '.join(domain_list)])

if __name__=="__main__":
    execpath=os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "domainparser2.LINUX")
    dssp_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "dssp")
    log='DomainParser.log'

    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-execpath="):
            execpath=os.path.abspath(arg[len("-execpath="):])
        elif arg.startswith("-dssp_path="):
            dssp_path=os.path.abspath(arg[len("-dssp_path="):])
        elif arg.startswith("-log="):
            log=arg[len("-log="):]
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            if not os.path.isfile(arg):
                sys.stderr.write("ERROR! No such file %s\n"%arg)
                exit()
            argv.append(arg)
    
    txt=''
    for arg in argv:
        txt+=DomainParser(arg,execpath,dssp_path)+'\n'

    if log and log!='-':
        fp=open(log,'w')
        fp.write(txt)
        fp.close()
    else:
        sys.stdout.write(txt)
