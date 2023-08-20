#!/bin/bash

# If you have wget installed, this bash script will download 
# topofiles need for modeling Seaside with CSZ events

REMOTE=http://depts.washington.edu/clawpack/geoclaw/topo/copes_hub/
TOPODIR=./topofiles

FILES="SeasideN_13s_mhw.asc 
    SeasideN_13s_mhw.png
    SeasideS_13s_mhw.asc
    SeasideS_13s_mhw.png
    astoria_2s_mhw.asc
    crm_vol8_3sec_cropped_for_seaside.asc
    crm_3sec_isobaths.png
    "
    

for FILE in $FILES
do
    if [ -f "$TOPODIR/$FILE" ]; then
        echo "$TOPODIR/$FILE exists, not fetching."
    else 
        wget $REMOTE/$FILE -P $TOPODIR
    fi
done

