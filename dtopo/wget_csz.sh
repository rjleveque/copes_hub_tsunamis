#!/bin/bash

# If you have wget installed, this bash script will download 
# topofiles need for modeling Seaside with CSZ events

# Go to the REMOTE directory in a web browser to see what other files
# are available.

# Note that the files downloaded have xlower = 233 and this must be shifted
# by -360 to -127 by modifying the fourth line of the file before using.

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

