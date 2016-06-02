#!/bin/bash
# grep Gene Ontology code
PATTERN="GO[:][\d]{7}"
if [ $# -eq 0 ];then
    cat|grep -ohP $PATTERN|sort -n |uniq
else
    grep -ohP $PATTERN $@ |sort -n |uniq
fi
