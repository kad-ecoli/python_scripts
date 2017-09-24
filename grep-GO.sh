#!/bin/bash
# grep Gene Ontology code
PATTERN="\bGO[:][\d]{7}\b"
if [ $# -eq 0 ];then
    cat|grep -ohP $PATTERN|sort -n |uniq
else
    grep -ohP $PATTERN $@ |sort -n |uniq
fi
