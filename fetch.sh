#!/bin/bash
# 2015-11-15 Chengxin Zhang
# fetch PDB files from internet

##### MIRROR and FORMAT START ####
    MIRROR=ftp://ftp.wwpdb.org/pub
    #MIRROR=ftp://ftp.ebi.ac.uk/pub/databases
    #MIRROR=ftp://ftp.pdbj.org/pub

    DIR=pdb/data/structures/all

    #FORMAT=XML
    #FORMAT=XML-extatom
    #FORMAT=XML-noatom
    #FORMAT=mmCIF
    #FORMAT=nmr_chemical_shifts
    #FORMAT=nmr_restraints
    #FORMAT=nmr_restraints_v2
    FORMAT=pdb
    #FORMAT=structure_factors
##### MIRROR and FORMAT END ####
##### 3rd Party software START ####
    PYMOL=pymol
##### 3rd Party software END ####

if [ `echo $1|grep -P '^[\d][\w]{3}$'` ];then  # PDB coordinate
    PDBid=$(echo $1|tr "[:upper:]" "[:lower:]")
    if [ -f "$PDBid".pdb ];then
        echo "$PDBid".pdb
    else
        FILE=pdb$PDBid.ent.gz
        URLtoCheck=$MIRROR/$DIR/$FORMAT/$FILE
        if [ `wget -O/dev/null -q $URLtoCheck && echo exists` ];then
	    wget $URLtoCheck &>/dev/null 
	    cat $FILE|gunzip > $PDBid.pdb
	    rm $FILE
            echo "$PDBid".pdb
        else
    	    echo "ERROR! Cannot fetch PDB $PDBid"
        fi
    fi
elif [ `echo $1|grep -P '^[\d][\w]{4}$'` ];then  # PDB chain
    PDBid=$(echo $1|tr "[:upper:]" "[:lower:]"|grep -ohP "^[\d][\w]{3}")
    CHAIN=$(echo $1|tr "[:lower:]" "[:upper:]"|grep -ohP "[\w]$")
    OUTPUT="$PDBid$CHAIN".pdb   # output file name
    SELE="$PDBid"_"$CHAIN"      # PYMOL selection name
    echo "$0 $PDBid"
    $0 $PDBid
    if [ -f "$PDBid".pdb ];then
        TMPDIR=/tmp/$USER/split_chain
        if [ ! -d $TMPDIR ];then # create tmp folder
            mkdir -p $TMPDIR
        fi
        echo "
load $PDBid.pdb
remove resn hoh
split_chain $PDBid
save $OUTPUT, $SELE
" > "$TMPDIR/$PDBid".pml
        $PYMOL -c $TMPDIR/$PDBid.pml
    fi
elif [ `echo $1|grep -iP '^[TR][\d]{4}$'` ];then # CASP target
    PDBid=$(echo $1|tr "[:lower:]" "[:upper:]")
    MIRROR=http://zhanglab.ccmb.med.umich.edu
    CASP_lst="8 9 10 11 12"
    DIR=decoys/casp
    FILE=$PDBid.native.pdb
    for casp in $CASP_lst;do
    	URLtoCheck=$MIRROR/$DIR$casp/$FILE
	if [ `wget -O/dev/null -q $URLtoCheck && echo exists` ];then
	    echo "CASP$casp target"
	    wget $URLtoCheck &>/dev/null 
	    mv $FILE $PDBid.pdb
            echo "$PDBid".pdb
	    exit
	fi
    done
elif [ `echo $1|grep -P '[\w]{1,10}_[\w]{1,5}'` ] || [ `echo $1|grep -P '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'` ];then # uniprot Accession number
    echo "Uniprot accession"
    MIRROR="http://www.uniprot.org/uniprot"
    URLtoCheck=$MIRROR/$1
    if [ `wget -O/dev/null -q $URLtoCheck && echo exists` ];then
	wget $URLtoCheck &>/dev/null 
	PDB_lst=$(cat $1|grep -ohP 'www.ebi.ac.uk\/pdbe-srv\/view\/entry\/[\d][\w]{3}'|grep -ohP '[\d][\w]{3}$'|sort -n|uniq)
	rm $1
        echo "$1"
	for pdb in $PDB_lst;do
	    #echo $pdb >> $1.list
	    $0 $pdb
	done
    else
    	echo "ERROR! Cannot fetch UniProt $1"
    fi
elif [ `echo $1|grep -P '^[Pp][Ff][\d]{5}'` ];then # Pfam family
    echo "PFAM"
    MIRROR=ftp://ftp.ebi.ac.uk/pub/databases/Pfam/mappings
    if [ ! -f "pdb_pfam_mapping.txt" ];then
        echo "wget $MIRROR/pdb_pfam_mapping.txt"
        wget $MIRROR/pdb_pfam_mapping.txt &>/dev/null 
    fi
    echo "pdb pfam mapping list ready"
    PDB_lst=$(cat pdb_pfam_mapping.txt|grep -i $1|grep -ohP '^[\d][\w]{3,4}\b'|sort -n|uniq)
    for pdb in $PDB_lst;do
        #echo $pdb
        $0 $pdb
    done
else
    echo "fetch 1ANK"
    echo "    Fetch PDB coodinate file 1ank.pdb from PDB FTP"
    echo "fetch T0388"
    echo "    Fetch native PDB coodinate file T0388.pdb in CASP from zhanglab"
    echo "fetch KAD_ECOLI"
    echo "    Fetch all PDB files for UniProtKB entry KAD_ECOLI"
    echo "fetch P76347"
    echo "    Fetch all PDB files for uniprot accession P76347"
    echo "fetch P76347"
    echo "    Fetch all PDB files for uniprot accession P76347"
    echo "fetch PF00001"
    echo "    Fetch all PDB files for Pfam family PF00001"
    exit
fi

if [ $2 ];then
    shift
    $0 $@
fi
if [ -f $TMPDIR ];then
    rm -rf $TMPDIR
fi
