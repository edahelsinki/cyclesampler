#!/bin/bash
URL=http://files.grouplens.org/datasets/hetrec2011/hetrec2011-lastfm-2k.zip
FNAME_ZIP=hetrec2011-lastfm-2k.zip
FNAME_2=user_artists.dat
DESTDIR=$1

set -x
mkdir $DESTDIR
wget -O $DESTDIR/$FNAME_ZIP $URL
unzip $DESTDIR/$FNAME_ZIP -d $DESTDIR

tr -s "\t" ";" < $DESTDIR/$FNAME_2 | tail -n +2 > lastfm_listening_count.csv
