#!/bin/bash

# If you have wget installed, this bash script will download 
# topofiles need for modeling Seaside with CSZ events

REMOTE=http://depts.washington.edu/ptha/dtopo/CSZB
DTOPODIR=./dtopofiles

FILES="CSZ_SM1.tt3
    CSZ_M1.tt3
    "
    

for FILE in $FILES
do
    if [ -f "$DTOPODIR/$FILE" ]; then
        echo "$DTOPODIR/$FILE exists, not fetching."
    else 
        wget $REMOTE/$FILE -P $DTOPODIR
    fi
done

