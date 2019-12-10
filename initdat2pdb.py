#!/usr/bin/env python
docstring='''
initdat2pdb.py [options] [seq.txt] init.dat
    Convert LOMETS threading result "init.dat" to PDB format C-alpha traces.
    Using FASTA format "seq.txt" as target sequence. If not specified,
    Recover target sequence from "init.dat". Unaligned positions in target 
    sequence will be represented by '-'

outputs:
    alignment.fasta     - multiple sequence alignment
    template[1-300].pdb - template structures
    target[1-300].pdb   - query-template alignments
    template.ent - multimodel PDB for all template structures 
    target.ent   - multimodel PDB for all target structures 

options:
    -prefix=''
        prefix for all output filenames
    -templates=10
        maximum number of template structures to parse.
        default is parsing all templates
    -good={false,true}
        whether only parse good templates (zscore>zscore0)
        false: (default) parse both good and bad templates
        true:  only parse good templates, autodetect threader
        dat:   only parse good templates from init.dat
        threader: only parse good templates from init.threader
    -template_only=LOMETS_
        only regenerate the template structures, using "LOMETS"
        as file suffix.
'''
import sys,os
import shutil
import re
import random
from fixMSE import code_with_modified_residues

zscore0={ # dict for good template zscore cutoff
    "MUS":6.1    , "QQQ":6.0    , "GGGd":12.0 , "JJJb":8.7    , 
    "RRR6":9.8   , "SPX":6.9    , "VVV":7.0   , "WWW":6       , 
    "HHP":11     , "OOO":25     , "BBB":3.2   , "RRR3":18     , 
    "IIIe":10    , "IIIj":15    , "PRC":21    , "FRM":4.8     , 
    "FF3":33     , "RAP":7      , "RAP2":6.8  , "pgen":6.3    ,   
    "mgen":5.2   , "phyre2":97.0, "hhpred":100, "hhpredo":101 , 
    "ROS":-1     , "QUA":-1     , "RQ":-1     , "RQ2":-1      ,

    "MUSTER" :5.8, "dPPAS":9.3  , "wdPPAS":9.0, "wPPAS":7.5   , 
    "wMUSTER":6.0, "dPPAS2":10.5, "PPAS":7.0  , "Env-PPAS":8.0,
    }

def read_one_sequence(infile="seq.txt"):
    '''Read one sequence from FASTA/plain-text format sequence file'''
    fp=open(infile,'rU')
    txt=fp.read().strip()
    fp.close()
    if txt.startswith('>'): # FASTA format
        lines=[e for e in txt.split('>') if e.strip()][0].splitlines()[1:]
    else: # plain text
        lines=txt.splitlines()
    sequence=''.join([line.strip() for line in lines])
    return sequence

def split_initdat(initdat="init.dat",templates=0):
    '''read LOMETS result "init.dat"
    only read the first "templates" structures if "templates" is not zeros
    return a list of plat text for each template
    '''
    fp=open(initdat,'rU')
    txt=fp.read()
    fp.close()

    # pattern for start section of one template section
    pattern=re.compile("\n\s*\d+\s+[-]{0,1}[.\d]+\s+\d+\s+\w+[\w\W]+?\n")
    # index for start sections of each decoy (except the first section)
    match_index=[e.start()+1 for e in pattern.finditer(txt)]
    if not match_index:
        return []

    initdat_txt_list=[txt[:match_index[0]]] + \
        [txt[match_index[idx]:match_index[idx+1]] for idx in \
            range(len(match_index)-1)]+[txt[match_index[-1]:]]
    if templates:
        initdat_txt_list=initdat_txt_list[:templates+1]
    return initdat_txt_list

def convert_initdat_txt(sequence='',initdat_txt='',CONECT=False):
    '''Using "sequence" as target sequence
    Convert LOMETS init.dat text to PDB format text.

    Return Values:
        template_pdb_txt- text for templates in PDB format
        target_pdb_txt  - text for unrefined targets in PDB format
        alignment       - a list for target-template aligment. target is the
                          first element. target will be guessed from 
                          initdat_txt and unaligned region will be marked
                          '-'. If "sequence" is not empty, use ungapped
                          region of "sequence" as target sequence

    Options:
        CONECT - whether to add "CONECT" entries
    '''
    template_pdb_txt=''
    target_pdb_txt=''
    template_sequence=''
    template_line=''
    for line in initdat_txt.splitlines():
        if not line.strip():
            continue
        if not line.startswith("ATOM  "):
            template_pdb_txt+=line+'\n'
            target_pdb_txt  +=line+'\n'
            continue
        # description for PDB coordinate format can be found at
        # http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM
        '''
COLUMNS        DATA TYPE       CONTENTS                            
-----------------------------------------------------------------------
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         Atom serial number.
13 - 16        Atom            Atom name.
17             Character       Alternate location indicator.
18 - 20        Residue name    Residue name.
22             Character       Chain identifier.
23 - 26        Integer         Residue sequence number.
27             AChar           Code for insertion of residues.
31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)       Occupancy.
61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
73 - 76        LString(4)      Segment identifier, left-justified.
77 - 78        LString(2)      Element symbol, right-justified.
79 - 80        LString(2)      Charge on the atom.
        '''
        
        record_name="ATOM  "
        #atom_serial_number=line[21:26]
        atom_serial_number=line[6:11]
        atom_name="  CA "
        alternate_location_indicator=" "
        residue_name=line[17:20]
        chain_identifier=" A"
        residue_sequence_number=line[22:26]
        residue_insertion_code='    '
        X_coordinate=line[30:38]
        Y_coordinate=line[38:46]
        Z_coordinate=line[46:54]

        occupancy="  1.00"
        Bfactor="  0.00"
        segment_identifier="          "
        element_symbol=" C"
        charge="  "

        residue_name_short=code_with_modified_residues[line[17:20]]

        template_segment=line[54:59]+ \
            code_with_modified_residues[line[60:63]]+"    " \
            if line[60:63] in code_with_modified_residues else "          "
        target_segment=line[21:26]+ \
            residue_name_short+"    "
            #code_with_modified_residues[line[17:20]]+"    "

        template_residue_name=line[60:63]
        if not line[60:63] in code_with_modified_residues:
            template_residue_name=residue_name
        template_resi=line[55:59]
        if not line[55:59].strip() or \
           not line[60:63] in code_with_modified_residues:
            template_resi=residue_sequence_number
        if template_line and int(template_resi)<=int(template_line[22:26]):
            template_resi=str(1+int(template_line[22:26]))
            template_resi=' '*(4-len(template_resi))+template_resi


        template_line= \
            record_name                 + \
            atom_serial_number          + \
            atom_name                   + \
            alternate_location_indicator+ \
            template_residue_name       + \
            chain_identifier            + \
            template_resi               + \
            residue_insertion_code      + \
            X_coordinate                + \
            Y_coordinate                + \
            Z_coordinate                + \
            occupancy                   + \
            Bfactor                     + \
            target_segment              + \
            element_symbol              + \
            charge

        target_line  =                    \
            record_name                 + \
            line[21:26]                 + \
            atom_name                   + \
            alternate_location_indicator+ \
            residue_name                + \
            chain_identifier            + \
            residue_sequence_number     + \
            residue_insertion_code      + \
            X_coordinate                + \
            Y_coordinate                + \
            Z_coordinate                + \
            occupancy                   + \
            Bfactor                     + \
            template_segment            + \
            element_symbol              + \
            charge

        template_pdb_txt+=template_line+'\n'
        target_pdb_txt  +=target_line  +'\n'

        resi=int(residue_sequence_number)-1
        sequence+='-'*(resi-len(sequence))
        sequence=sequence[:resi]+residue_name_short+sequence[resi+1:]

        template_sequence+='-'*(resi-len(template_sequence))
        if line[60:63] in code_with_modified_residues:
            template_sequence=template_sequence[:resi]+ \
                code_with_modified_residues[line[60:63]]+template_sequence[resi+1:]
        else:
            template_sequence=template_sequence[:resi]+ \
                code_with_modified_residues[residue_name]+template_sequence[resi+1:]

    if CONECT: # add CONECT entries
        '''
COLUMNS       DATA  TYPE      FIELD        DEFINITION
-------------------------------------------------------------------------
 1 -  6        Record name    "CONECT"
 7 - 11        Integer        serial       Atom  serial number
12 - 16        Integer        serial       Serial number of bonded atom
17 - 21        Integer        serial       Serial number of bonded atom
22 - 26        Integer        serial       Serial number of bonded atom
27 - 31        Integer        serial       Serial number of bonded atom
        '''
        template_atom_serial_number_list=[int(line[6:11]) for line in \
            template_pdb_txt.splitlines() if line.startswith("ATOM  ")]
        target_atom_serial_number_list  =[int(line[6:11]) for line in \
            target_pdb_txt.splitlines()   if line.startswith("ATOM  ")]
        for idx in range(1,len(template_atom_serial_number_list)):
            template_serial_prev=str(template_atom_serial_number_list[idx-1])
            target_serial_prev  =str(target_atom_serial_number_list[idx-1])

            template_serial     =str(template_atom_serial_number_list[idx])
            target_serial       =str(target_atom_serial_number_list[idx])

            record_name="CONECT"
            template_line=record_name         + \
                ' '*(5-len(template_serial_prev))+template_serial_prev+ \
                ' '*(5-len(template_serial     ))+template_serial     + \
                ' '*64

            target_line  =record_name         + \
                ' '*(5-len(target_serial_prev  ))+target_serial_prev  + \
                ' '*(5-len(target_serial       ))+target_serial       + \
                ' '*64
            
            template_pdb_txt+=template_line+'\n'
            target_pdb_txt  +=target_line  +'\n'

    template_sequence+='-'*(len(sequence)-len(template_sequence))
    return template_pdb_txt,target_pdb_txt,[sequence,template_sequence]

def initdat2pdb(sequence='', initdat="init.dat",
    prefix='', templates=0, good=False, template_only=''):
    '''Convert LOMETS result "initdat" into C-alpha traces

    sequence  - target sequence
    initdat   - file name for LOMETS template file
    prefix    - prefix for all output filenames
    templates - maximum number of template structures to parse.
                Parsing all templates if set to 0
    good      - whether to parse good templates only
                False: all templates both good and bad
                "dat": good templates in init.dat
                threader: good templates in init.threader
    template_only - whether to generate templates only

    return a dictinary "output_dict" for output files
    output_dict["sequence"]  - target sequence
    output_dict["alignment"] - multiple sequence alignment FASTA file
    output_dict["template"]  - template structure PDB files
    output_dict["target"]    - target structure PDB files
    '''
    output_dict={
        "sequence" : sequence,
        "alignment": prefix+"alignment.fasta",
        "template" :[prefix+"template.ent"],
        "target"   :[prefix+"target.ent"],
    }
    
    tmp_dir="/tmp/"+os.getenv("USER")+'/initdat'+prefix+ \
       str(random.randint(1000,9999))+'/'
    scratch_dir="/scratch/%s/%s"%(os.getenv("USER"),os.getenv("SLURM_JOBID"))
    if os.getenv("SLURM_JOBID") and os.path.isdir(scratch_dir):
         tmp_dir=scratch_dir+'/initdat'+prefix+ \
            str(random.randint(1000,9999))+'/'
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    initdat_tmp=tmp_dir+"init.dat"
    shutil.copy(initdat,initdat_tmp)
    
    initdat_txt_list=split_initdat(initdat=initdat_tmp,templates=templates)
    # first element of initdat_txt_list is a short summary for LOMETS run
    templates=len(initdat_txt_list)-1

    template_ent_tmp=tmp_dir+"template.ent"
    target_ent_tmp  =tmp_dir+"target.ent"
    template_ent_txt=initdat_txt_list[0]
    target_ent_txt  =initdat_txt_list[0]

    alignment_txt=''
    full_alignment=[sequence] if sequence else ['']
    alignment_fasta_tmp=tmp_dir+"alignment.fasta"

    for idx,initdat_txt in enumerate(initdat_txt_list[1:]):
        template_idx=str(idx+1)

        # first line that describes threading result
        thread_line=initdat_txt.splitlines()[0] 
        if good:
            zscore=float(thread_line.strip().split()[1])
            cutoff=zscore0[thread_line.split()[-1] if good=="dat" else good]
            if zscore<=cutoff:
                continue
        print thread_line

        MODEL_header="MODEL"+' '*(9-len(template_idx))+template_idx+'\n'
        template_ent_txt+=MODEL_header
        target_ent_txt  +=MODEL_header

        template_pdb_txt,target_pdb_txt,alignment=convert_initdat_txt(
            sequence=full_alignment[0],initdat_txt=initdat_txt,CONECT=True)

        full_alignment=[alignment[0]]+full_alignment[1:-1]+[alignment[-1]]
        alignment_txt+=">template%s.pdb\n%s\n"%(template_idx,alignment[-1])


        template_ent_txt+=template_pdb_txt
        target_ent_txt  +=target_pdb_txt

        template_ent_txt+="ENDMDL\n"
        target_ent_txt  +="ENDMDL\n"
        
        if not template_only:
            target_pdb_tmp  =tmp_dir+"target"  +template_idx+".pdb"
            fp=open(target_pdb_tmp  ,'w')
            fp.write(target_pdb_txt)
            fp.close()
            target_pdb_file  =prefix+"target"  +template_idx+".pdb"
            shutil.copy(target_pdb_tmp  ,target_pdb_file  )
            output_dict["target"  ].append(target_pdb_file  )

        template_pdb_tmp=tmp_dir+"template"+template_idx+".pdb"
        fp=open(template_pdb_tmp,'w')
        fp.write(template_pdb_txt)
        fp.close()
        template_pdb_file=prefix+"template"+template_idx+".pdb" \
            if not template_only else prefix+template_only+template_idx+".pdb"
        shutil.copy(template_pdb_tmp,template_pdb_file)
        output_dict["template"].append(template_pdb_file)

    template_ent_txt+="END\n"
    target_ent_txt  +="END\n"
    if not template_only:
        fp=open(template_ent_tmp,'w')
        fp.write(template_ent_txt)
        fp.close()
        fp=open(target_ent_tmp  ,'w')
        fp.write(target_ent_txt)
        fp.close()

        shutil.copy(template_ent_tmp,prefix+"template.ent")
        shutil.copy(target_ent_tmp  ,prefix+"target.ent")

    if not len(initdat_txt_list):
        txt=''
    else:
        if not template_only:
            txt=">target\n%s\n%s"%(alignment[0],alignment_txt)
        else:
            txt=">model.pdb\n%s\n"%alignment[0]
            for b,block in enumerate(alignment_txt.split('>')):
                if not block.strip():
                    continue
                name,sequence=block.splitlines()
                txt+=">%s%d.pdb\n%s\n"%(template_only,b,sequence)

    fp=open(alignment_fasta_tmp,'w')
    fp.write(txt)
    fp.close()
    if not template_only:
        shutil.copy(alignment_fasta_tmp,prefix+"alignment.fasta")
    else:
        shutil.copy(alignment_fasta_tmp,prefix+template_only+".fasta")

    ## clean up temporary folder ##
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    return output_dict

if __name__=="__main__":
    prefix=''
    templates=0
    good=False
    template_only=''

    argv=[]
    for arg in sys.argv[1:]:
        if not arg.startswith('-'):
            argv.append(arg)
            continue
        if arg.startswith('-prefix='):
            prefix=arg[len("-prefix="):]
        elif arg.startswith("-templates="):
            templates=int(arg[len("-templates="):])
        elif arg.startswith("-template_only="):
            template_only=arg[len("-template_only="):
                ].strip('"').strip("'")
        elif arg.startswith("-good="):
            if arg[len("-good="):].lower()!="false":
                good=arg[len("-good="):]
        else:
            print >>sys.stderr, "ERROR! Unknown argument "+arg
            exit()

    if len(argv)<1:
        sys.stderr.write(docstring)
        exit()
    
    sequence='' if len(argv)==1 else read_one_sequence(argv[0])

    if str(good).lower()=="true":
        good=argv[-1][argv[-1].rfind('.')+1:]

    initdat2pdb(sequence=sequence,initdat=argv[-1],prefix='',
        templates=templates,good=good,template_only=template_only)
