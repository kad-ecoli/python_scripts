#!/usr/bin/env python
# 2016-02-27 Chengxin Zhang
docstring='''
fetch 4k5y
    Download PDB 4k5y to 4k5y.pdb

fetch 4k5yA
    Download PDB 4k5y chain A to 4k5yA.pdb

fetch 3pjpB2
    Download PDB 4k5y chain A domain 1 to 3pjpB2.pdb

fetch P34998
    Download all PDB associated with uniprot ID P34998

fetch PF00406
    Download all PDB associated with pfam ID PF00406

option:
    -outfmt={pdb,fasta,list} output format
        pdb: pdb coordinate
        fasta: fasta sequence for PDB or for uniprot accession
        list: list of available PDB for uniprot accession
    -include_model={false,true} output format
        whether to download obsolete PDB such as theoretical models
    -execpath=./domainparser2.LINUX
        path to DomainParser executable. By default it is guessed by 
        location of this script
    -dssp_path=./dssp
        path to DSSP executable. By default it is guessed by
        location of this script
    -pulchra_path=./pulchra
        path to pulchra executable. By default it is guessed by 
        location of this script. pulchra is used to construct full atom
        model from backbone model when the input structure contains too
        few atoms.
'''
import sys,os
import re
import gzip
import tarfile
import urllib
import urllib2
import shutil

from DomainParser import DomainParser

pdb_mirror="ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb"
pdb_bundle_mirror="ftp://ftp.wwpdb.org/pub/pdb/compatible/pdb_bundle/"
uniprot_mirror="http://www.uniprot.org/uniprot/"
pdbe_mirror="http://www.ebi.ac.uk/pdbe-srv/view/entry/"
rcsb_fasta_mirror="http://www.rcsb.org/pdb/files/fasta.txt?structureIdList="
rcsb_pdb_mirror="http://www.rcsb.org/pdb/files/"
pdb_pfam_mapping_mirror="http://www.rcsb.org/pdb/rest/hmmer?file=hmmer_pdb_all.txt"
obsolete_mirror="ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat"
large_split_mirror="ftp://ftp.wwpdb.org/pub/pdb/compatible/pdb_bundle/large_split_mapping.tsv"

pdb_pattern=re.compile("^\d\w{3}$")
pdb_chain_pattern=re.compile("^\d\w{4,5}$")
accession_pattern=re.compile(# uniprot accession (AC)
    "^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")
entry_pattern=re.compile("^\w{1,5}_\w{1,5}") # uniprot entry name
PDBe_entry_pattern=re.compile("www\.ebi\.ac\.uk\/pdbe\-srv\/view/entry/(\d\w{3})[\w\W]+?<td>([/\w]+)<\/td><td><a\s+href=\"\.+.\/blast\/\?about=\w+\[(\d+\-\d+)\]")
#DBREF_pattern=re.compile("DBREF\s\s\d\w{3}\s\w{0,1}\s+(\d+)\s+(\d+)\s+UNP\s+(\w+)\s+(\w+_\w+)\s+(\d+)\s+(\d+)\s*")
DBREF_pattern=re.compile("DBREF[\d]{0,1}\s{1,2}\d\w{3}\s\w{0,1}\s+(\d+)\s+(\d+)\s+UNP\s+(\w+)")
PFAM_pattern=re.compile("[Pp][Ff]\d{5}")
pdb_bundle_pattern=re.compile("\d\w{3}\-pdb\-bundle\d+\.pdb")

chainID_column={ # column to store chainID, 0 for no chainID
    "DBREF":12, "SEQADV":16, "SEQRES":11, "MODRES":16,# Primary Structure
    "HET":12,                             # Heterogen
    "HELIX":19, "SHEET":21,               # Secondary Structure
    "SSBOND":15, "CONECT":0,              # Connectivity
    "SITE":22,                            # Miscellaneous
    "MODEL":0, "ATOM":21, "TER":21, "HETATM":21, "ENDMDL":0, # Coordinate
    "END":0,                              # Bookkeeping
}

resSeq_column={ # columns for residue sequence number, 0 for no resSeq
    "DBREF":(14,18), "SEQADV":(18,22), "MODRES":(18,22), # Primary Structure
    "HELIX":(21,25), "SHEET":(22,26),       # Secondary Structure 
    "SSBOND":(17,21), "CONECT":0,           # Connectivity
    "SITE":(23,27),                         # Miscellaneous
    "MODEL":0, "ATOM":(22,26), "TER":0, "HETATM":(22,26), "ENDMDL":0, # Coordinate
    "END":0,                                # Bookkeeping
}

## class the store mapping between old PDB chain and new PDB chain
class large_split_mapping():
    '''class to store old and new large structure mapping
    key is old ID, value is new ID
    '''
    def __init__(self,large_split_mapping_file=''):
        '''initialize class
        if large_split_mapping_file is given, the mapping file is used.
        otherwise, download it from internet.
        '''
        self.new_pdb_list=[]
        self.new_chain_list=[]
        self.old2new=dict()
        self.new2old=dict()

        if not large_split_mapping_file:
            large_split_mapping_file=wget(large_split_mirror)
        fp=open(large_split_mapping_file,'rU')
        txt=fp.read()
        fp.close()

        for line in txt.splitlines():
            if line.startswith("Large") or not line.strip():
                continue
            newID,oldID_str=line.split('\t')
            if ':' in newID:
                newID=newID.replace(':','')
                oldID_str=oldID_str.replace(':','')
                self.new_chain_list.append(newID)
            else:
                self.new_pdb_list.append(newID)
                
            self.new2old[newID]=oldID_str
            for oldID in oldID_str.split(','):
                self.old2new[oldID]=newID

    def __repr__(self):
        '''print mapping file'''
        txt="Large Structure PDB ID\tSplit PDB IDs\n"
        txt+='\n'.join([ID+'\t'+self.new2old[ID] for ID in self.new_pdb_list])
        txt+="\n\nLarge Structure chain IDs\nSplit chain ID\n"
        txt+='\n'.join([ID+'\t'+self.new2old[ID] for ID in self.new_chain_list])
        txt+='\n'
        return txt

class obsolete_dat(dict):
    '''class to store obsolete PDB mapping.
    key is old ID, value is a list of new ID'''
    def __init__(self,obsolete_dat_file=''):
        '''initialize class
        if obsolete_dat_file is given, the mapping file used
        otherwise, download it from internet.
        '''
        if not obsolete_dat_file:
            obsolete_dat_file=wget(obsolete_mirror)
        fp=open(obsolete_dat_file,'rU')
        txt=fp.read()
        fp.close()
        for line in txt.splitlines()[1:]:
            line=line.split()
            if len(line)<3:
                continue
            self[line[2]]=line[3:]

    def __repr__(self):
        '''print mapping file'''
        txt='\n'.join([ID+'\t'+','.join(self[ID]) for ID in self])+'\n'
        return txt

def wget(url='',outfile='',no_err=False,show_url=False):
    '''retrieve file from internet if not exists at current directory.
    return file output filename if download successfully.
    return empty string otherwise.

    outfile - output file name. By default it is the basename
    no_err - whether supress downloading error
    show_url - whether print url at stdout
    '''
    if not outfile:
        outfile=os.path.basename(url)
    if not os.path.isfile(outfile):
        if show_url:
            sys.stderr.write("fetching %s\n"%url)
        try:
            urllib.urlretrieve(url, outfile)
        except Exception,err:
            if not err:
                sys.stderr.write(str(err)+'\n')
            return ''
    elif show_url:
        sys.stderr.write("%s already exists\n"%outfile)
    return outfile

def fetch_bundle(PDBid,no_err=False):
    '''fetch Best effort/minimal PDB format files for large structures
    return the tarball file name. 
    no_err - whether supress downloading error'''
    # As large structures (containing >62 chains and/or 99999 ATOM lines) 
    # cannot be represented in the legacy PDB file format, data are 
    # available in the PDB archive as single PDBx/mmCIF files representing
    # the entire structure, and as TAR files containing a collection of 
    # best effort/minimal files in the PDB file format to support users 
    # and software tools that rely solely on the PDB file format. 
    PDBid=PDBid.lower()
    tarball_name=PDBid+"-pdb-bundle.tar.gz"
    return wget(pdb_bundle_mirror+PDBid[1:3]+'/'+PDBid+'/'+tarball_name
        ,tarball_name,no_err)

def extract_chain_from_bundle(tarball_name,chainID_list=[]):
    '''split best effort/minimal PDB tarball into PDB files containing 
    individual chainID'''
    # Best effort/minimal PDB format files contain only HEADER, AUTHOR, 
    # JRNL, CRYST1, SCALEn, ATOM, HETATM records.
    # The following are not included: OBSLTE, TITLE, CAVEAT, COMPND,
    # SOURCE, KEYWDS, EXPDTA, REVDAT, SPRSDE, REMARKS, DBREF, SEQADV, 
    # SEQRES, MODRES, HET, HETNAM, HETSYN, FORMUL, HELIX, SHEET, SSBOND, 
    # LINK, CISPEP, SITE, ORIGXn, MTRIXn, CONECT.
    PDBid=tarball_name[:4]
    chain_id_mapping=dict()
    PDB_chain_file_list=[]
    tar=tarfile.open(tarball_name,'r:gz')

    # parse chain id mapping
    fp=tar.extractfile(PDBid+"-chain-id-mapping.txt")
    map_txt=fp.read()
    fp.close()
    for section in map_txt.split('\n'+PDBid+"-pdb-bundle"):
        if not ':' in section:
            continue
        idx,section=section.split(".pdb:\n")
        pdb_bundle_name=PDBid+"-pdb-bundle"+idx+".pdb"
        for line in section.splitlines():
            New_chain_ID,Original_chain_ID=line.split()
            chain_id_mapping[Original_chain_ID]=(
                pdb_bundle_name,New_chain_ID)

    # parse each chain
    PDB_chain_file_list=[]
    for Original_chain_ID in chainID_list:
        if not Original_chain_ID in chain_id_mapping:
            sys.stderr.write("ERROR! no chain %s in %s\n"%(
                Original_chain_ID,PDBid))
            continue
        pdb_bundle_name,New_chain_ID=chain_id_mapping[Original_chain_ID]
        PDB_chain_file=PDBid+Original_chain_ID+".pdb"
        PDB_chain_file_list.append(PDB_chain_file)
        if os.path.isfile(PDB_chain_file):
            continue
       
        # parse text in *-pdb-bundle*.pdb
        fp=tar.extractfile(pdb_bundle_name)
        pdb_lines=fp.read().splitlines()
        fp.close()

        # extract text for specific chain
        PDB_txt=''    
        for line in pdb_lines:
            section=line.split()[0]
            if section in chainID_column and (not chainID_column[section] \
                or line[chainID_column[section]]==New_chain_ID):
                PDB_txt+=line+'\n'

        # write PDB for individual chain
        fp=open(PDB_chain_file,'w')
        fp.write(PDB_txt)
        fp.close()
    return PDB_chain_file_list

def obsolete2supersede(PDBid):
    '''get superseding PDB id instead of obsolete PDBid'''
    if not "obsolete_dict" in globals():
        global obsolete_dict
        obsolete_dict=obsolete_dat()
    PDBid=PDBid.upper()
    if PDBid in obsolete_dict and obsolete_dict[PDBid]:
        return obsolete_dict[PDBid][0]
    return ''

def fetch(PDBid,include_model=False):
    '''download PDBid.pdb if standard PDB format is present
    download superseding PDB if the provided PDBid is obsolete
    download the PDBid-pdb-bundle.tar.gz if only best 
    effort/minimal PDB format files exists

    include_model - whether download obsolete PDB such as theoretical model
    '''
    PDBid=PDBid.lower()
    PDB_file=PDBid+".pdb"
    PDB_gz_file=PDB_file+".gz"
    if os.path.isfile(PDB_file):
        return PDB_file # skip already downloaded PDB file

    PDB_gz_file=wget(pdb_mirror+PDBid+".ent.gz")

    if not PDB_gz_file and include_model: # download from rcsb website
        PDB_gz_file=wget(rcsb_pdb_mirror+PDBid.upper()+".pdb.gz",PDB_gz_file)
        if not os.path.getsize(PDB_gz_file):
            os.remove(PDB_gz_file)
            PDB_gz_file=''

    if not PDB_gz_file: # download best effort/minimal PDB bundle
        tarball_name=fetch_bundle(PDBid,no_err=True)
        if tarball_name:
            return tarball_name

    if not PDB_gz_file: # download superseding PDB instead of obsolete PDB
        PDBid_supersede=obsolete2supersede(PDBid)
        if PDBid_supersede:
            sys.stderr.write("%s superseded by %s\n"%(PDBid,PDBid_supersede))
            PDB_file_supersede=fetch(PDBid_supersede,include_model)
            if PDB_file_supersede:
                return PDB_file_supersede

    if not PDB_gz_file: # cannot fetch anything: do not fetch computation model
        return ''

    fp_gz=gzip.open(PDB_gz_file,'rb')
    fp_pdb=open(PDB_file,'w')
    shutil.copyfileobj(fp_gz,fp_pdb)
    fp_pdb.close()
    fp_gz.close()
    os.unlink(PDB_gz_file)
    return PDB_file

def obsolete_chain2supersede_chain(PDB_chain):
    '''get superseding PDB chain instead of obsolete PDBid chain'''
    if not "large_split_dict" in globals(): # check if it is initialized
        global large_split_dict
        large_split_dict=large_split_mapping()
    if PDB_chain in large_split_dict.old2new:
        return large_split_dict.old2new[PDB_chain]
    
    PDBid_supersede=obsolete2supersede(PDB_chain[:4])
    if not PDBid_supersede:
        return ''
    return PDBid_supersede+PDB_chain[4:]

def fetch_chain(PDB_chain,include_model=False,
    execpath="domainparser2.LINUX",dssp_path="dssp",pulchra_path="pulchra"):
    '''download PDB and split into specific chain
    include_model - whether download obsolete PDB such as theoretical model
    '''
    PDB_file=fetch(PDB_chain[:4],include_model)
    if not PDB_file:
        return ''
    if not PDB_file[:4]==PDB_chain[:4]: # PDB_chain was obsolete
        PDB_chain=obsolete_chain2supersede_chain(PDB_chain)

    PDB_chain_file_list=extract_chain(PDB_file,PDB_chain[4:],
        execpath,dssp_path,pulchra_path)
    return PDB_chain_file_list

def extract_chain(PDB_file,chainID_list=[],
    execpath="domainparser2.LINUX",dssp_path="dssp",pulchra_path="pulchra"):
    '''split PDB_file into PDB files containing individual chainID'''
    if isinstance(chainID_list,str):
        chainID_list=[chainID_list]
    if PDB_file.endswith(".tar.gz"):
        return extract_chain_from_bundle(PDB_file,chainID_list)
    fp=open(PDB_file,'rU')
    pdb_lines=fp.read().splitlines()
    fp.close()
    PDB_chain_file_list=[]
    for chainID in chainID_list:
        PDB_chain_file=PDB_file.split('.')[0]+chainID+".pdb"
        PDB_chain_file_list.append(PDB_chain_file)
        if os.path.isfile(PDB_chain_file):
            continue
        PDB_txt=''
        atom_serial_list=[]
        for line in pdb_lines:
            section=line.split()[0]
            if section in chainID_column and (not chainID_column[section] \
                or line[chainID_column[section]]==chainID):
                if section in ("ATOM","HETATM"):
                    atom_serial_list.append(line[6:11])
                elif section=="CONECT" and not line[6:11] in atom_serial_list:
                    continue # CONECT is after all ATOM & HETATM
                PDB_txt+=line+'\n'
        if not atom_serial_list:
            if re.match("^\w+?(\d+)$",chainID):
                chainID,domainID=re.findall("^(\w+?)(\d+)$",chainID)[0]
                tmp_PDB_chain_file_list=extract_chain(PDB_file,[chainID])
                if tmp_PDB_chain_file_list:
                    log=DomainParser(tmp_PDB_chain_file_list[0],
                        execpath,dssp_path,pulchra_path)
                    continue
            sys.stderr.write("ERROR! no chain %s in %s\n"%(
                chainID,PDB_file.split('.')[0]))
            PDB_chain_file_list=PDB_chain_file_list[:-1]
            continue
        fp=open(PDB_chain_file,'w')
        fp.write(PDB_txt)
        fp.close()
    return PDB_chain_file_list

def fetch_pdb_list_for_pfam(pfamID):
    '''retrieve a list of all available PDB for a pfam ID'''
    pdb_pfam_mapping_txt="pdb_pfam_mapping.txt"
    if not os.path.isfile(pdb_pfam_mapping_txt):
        wget(pdb_pfam_mapping_mirror,pdb_pfam_mapping_txt)
    fp=open(pdb_pfam_mapping_txt,'rU')
    txt=fp.read()
    fp.close()
    pdb_pfam_mapping_pattern=re.compile(
        "\n(\d\w{3})\s(\w)\s(\d+)\s(\d+)\s"+pfamID)
    return [(e[0],e[1],e[2]+'-'+e[3]) for e in \
        pdb_pfam_mapping_pattern.findall(txt)]

def fetch_pdb_list_for_uniprot(accession,all_chain=True):
    '''retrieve a list of all available PDB chain for a uniprot accession
    all_chain - True:  (default) all available chains
              - False: only the first chain if multiple chains are 
                       present in one PDB 
    '''
    fp=urllib2.urlopen(uniprot_mirror+accession)
    txt=fp.read()
    fp.close()
    pdbe_list=sorted(set(PDBe_entry_pattern.findall(txt)))
    if not all_chain:
        for i,pdbe in enumerate(pdbe_list):
            pdbe_list[i]=(pdbe[0],pdbe[1].split('/')[0],pdbe[2])
    return pdbe_list

def extract_uniprot(PDB_chain_file,accession,PERMISSIVE=True):
    '''extract the portion of PDB for specific uniprot "accession"
    according to DBREF. 

    PERMISSIVE: what to do when the only DBREF does not match accession
        True - If only one UNP is found, copy PDB_chain_file to 
               PDB_accession_file
        False - Raise error if the only DBREF UNP does not match accession
    '''
    accession=accession.strip()
    PDB_accession_file=PDB_chain_file.split('.')[0]+'_'+accession+".pdb"
    if os.path.isfile(PDB_accession_file):
        return PDB_accession_file
    fp=open(PDB_chain_file,'rU')
    PDB_txt=fp.read()
    fp.close()
    DBREF_list=DBREF_pattern.findall(PDB_txt)
    UNP_set=list(set([DBREF[2] for DBREF in DBREF_list]))
    if not UNP_set:
        shutil.copy(PDB_chain_file,PDB_accession_file)
        return PDB_accession_file

    if PERMISSIVE and len(UNP_set)==1:
        accession=UNP_set[0]
    
    resSeq_range=[]
    for DBREF in DBREF_list:
        if accession in DBREF[2:4]:
            resSeq_range+=range(int(DBREF[0]),int(DBREF[1])+1)
    if not resSeq_range:
        print >>sys.stderr,"ERROR! Cannot parse "+PDB_accession_file
        return ''
    resSeq_range=[str(resSeq) for resSeq in set(resSeq_range)]

    PDB_lines=PDB_txt.splitlines()
    PDB_txt=''
    atom_serial_list=[]
    for line in PDB_lines:
        section=line.split()[0]
        if section in resSeq_column and (not resSeq_column[section] \
        or line[resSeq_column[section][0]:resSeq_column[section][1]
           ].strip() in resSeq_range):
            if section in ("ATOM","HETATM"):
                atom_serial_list.append(line[6:11])
            elif section=="CONECT" and not line[6:11] in atom_serial_list:
                continue # CONECT is after all ATOM & HETATM
            PDB_txt+=line+'\n'
    fp=open(PDB_accession_file,'w')
    fp.write(PDB_txt)
    fp.close()
    return PDB_accession_file

def fetch_pfam(pfamID,
    execpath="domainparser2.LINUX",dssp_path="dssp",pulchra_path="pulchra"):
    '''download all available pdb for a pfam accession'''
    pdbe_list=fetch_pdb_list_for_pfam(pfamID)
    PDB_chain_file_list=[]
    for pdbe in pdbe_list:
        PDB_chain=pdbe[0].lower()+pdbe[1]
        PDB_chain_file_list+=fetch_chain(PDB_chain,
            execpath,dssp_path,pulchra_path)
    return PDB_chain_file_list

def fetch_uniprot(accession,all_chain=False,
    execpath="domainparser2.LINUX",dssp_path="dssp",pulchra_path="pulchra"):
    '''download all available pdb for a uniprot accession
    all_chain - False: (default)only the first chain if multiple 
                       chains are present in one PDB 
              - True:  all available chains
    '''
    pdbe_list=fetch_pdb_list_for_uniprot(accession,all_chain)
    PDB_chain_file_list=[]
    for pdbe in pdbe_list:
        if all_chain:
            for chain in pdbe[1].split('/'):
                PDB_chain=pdbe[0].lower()+chain
                PDB_chain_file_list+=fetch_chain(PDB_chain,
                    execpath,dssp_path,pulchra_path)
        else:
                PDB_chain=pdbe[0].lower()+pdbe[1].split('/')[0]
                PDB_chain_file_list+=fetch_chain(PDB_chain,
                    execpath,dssp_path,pulchra_path)
    PDB_accession_file_list=[extract_uniprot(PDB_chain_file,accession
        ) for PDB_chain_file in PDB_chain_file_list]
    return PDB_accession_file_list

def fetch_uniprot_sequence(accession):
    '''retrieve fasta sequence for accession'''
    fp=urllib2.urlopen(uniprot_mirror+accession+".fasta")
    txt=fp.read()
    fp.close()
    return txt

def fetch_pdb_sequence(PDBid):
    '''retrieve SEQRES sequence for PDBid'''
    fp=urllib2.urlopen(rcsb_fasta_mirror+PDBid)
    txt=fp.read()
    fp.close()
    return txt

def fetch_pdb_chain_sequence(PDB_chain):
    '''retrieve SEQRES sequence for PDB_chain'''
    PDBid=PDB_chain[:4].upper()
    chainID=PDB_chain[4:]
    fasta=fetch_pdb_sequence(PDB_chain[:4])
    txt='\n'.join(['>'+s for s in fasta.split('>') \
        if s.startswith(PDBid+':'+chainID)])
    return txt

if __name__=="__main__":
    outfmt="pdb"
    include_model=False
    execpath=os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "domainparser2.LINUX")
    dssp_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "dssp")
    pulchra_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),
        "pulchra")
    
    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-outfmt="):
            outfmt=arg[len("-outfmt="):]
        elif arg.startswith("-include_model="):
            include_model=(arg[len("-include_model="):].lower()=="true")
        elif arg.startswith("-execpath="):
            execpath=os.path.abspath(arg[len("-execpath="):])
        elif arg.startswith("-dssp_path="):
            dssp_path=os.path.abspath(arg[len("-dssp_path="):])
        elif arg.startswith("-pulchra_path="):
            log=os.path.abspath(arg[len("-pulchra_path="):])
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown argument %s\n"%arg)
            exit()
        else:
            argv.append(arg)
    if len(argv)<1:
        sys.stderr.write(docstring)
        exit()

    if outfmt=="pdb":
        for arg in argv:
            if pdb_pattern.match(arg):
                print fetch(arg,include_model)
            elif pdb_chain_pattern.match(arg):
                print '\n'.join(fetch_chain(arg,include_model,
                    execpath,dssp_path,pulchra_path))
            elif accession_pattern.match(arg) or entry_pattern.match(arg):
                print '\n'.join(fetch_uniprot(arg,
                    execpath,dssp_path,pulchra_path))
            elif PFAM_pattern.match(arg):
                print '\n'.join(fetch_pfam(arg,
                    execpath,dssp_path,pulchra_path))
    elif outfmt=="list":
        for arg in argv:
            if accession_pattern.match(arg) or entry_pattern.match(arg):
                print '\n'.join(['\t'.join(line) for line in \
                    fetch_pdb_list_for_uniprot(arg,all_chain=True)])
            elif PFAM_pattern.match(arg):
                print '\n'.join(['\t'.join(line) for line in \
                    fetch_pdb_list_for_pfam(arg)])
    elif outfmt=="fasta":
        for arg in argv:
            if accession_pattern.match(arg) or entry_pattern.match(arg):
                sys.stdout.write(fetch_uniprot_sequence(arg))
            elif pdb_pattern.match(arg):
                sys.stdout.write(fetch_pdb_sequence(arg))
            elif pdb_chain_pattern.match(arg):
                sys.stdout.write(fetch_pdb_chain_sequence(arg))
