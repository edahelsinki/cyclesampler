#!/bin/bash
URL=http://www2.informatik.uni-freiburg.de/~cziegler/BX/BX-CSV-Dump.zip
FNAME_ZIP=BX-CSV-Dump.zip
FNAME=BX-Book-Ratings.csv
FNAME2=BX-Book-Ratings_fixed.csv
DESTDIR=$1

set -x
mkdir $DESTDIR
wget -O $DESTDIR/$FNAME_ZIP $URL
unzip $DESTDIR/$FNAME_ZIP -d $DESTDIR

grep -via "[]?[]" $DESTDIR/$FNAME > $DESTDIR/$FNAME2

./preprocess_bxbooks.py $DESTDIR/$FNAME2 BX-Book-Ratings-processed.csv
