#!/bin/bash
URL=http://labrosa.ee.columbia.edu/millionsong/sites/default/files/challenge/train_triplets.txt.zip
FNAME_ZIP=train_triplets.txt.zip
FNAME=train_triplets.txt
DESTDIR=$1

set -x
mkdir $DESTDIR
wget -O $DESTDIR/$FNAME_ZIP $URL
unzip $DESTDIR/$FNAME_ZIP -d $DESTDIR
./preprocess_tasteprofile.py $DESTDIR/$FNAME tasteprofile_processed.csv
