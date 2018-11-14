#!/bin/bash
URL=https://snap.stanford.edu/data/finefoods.txt.gz
FNAME_ZIP=finefoods.txt.gz
FNAME=finefoods.txt
DESTDIR=$1

set -x
mkdir $DESTDIR
wget -O $DESTDIR/$FNAME_ZIP $URL
gunzip -c $DESTDIR/$FNAME_ZIP > $DESTDIR/$FNAME

./preprocess_finefoods.py < $DESTDIR/$FNAME > finefoods.csv
