#!/bin/bash
# grep UniProtKB accession numbers
PATTERN="\b[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\b"
if [ $# -eq 0 ];then
    cat|grep -ohP $PATTERN|sort -n |uniq
else
    grep -ohP $PATTERN $@ |sort -n |uniq
fi
