#!/bin/bash
set -x

DESTDIR=$1

mkdir $DESTDIR

wget -O $DESTDIR/ml-20m.zip http://files.grouplens.org/datasets/movielens/ml-20m.zip
unzip $DESTDIR/ml-20m.zip -d $DESTDIR
cut -d "," -f 1,2,3 $DESTDIR/ml-20m/ratings.csv | tr -s "," ";" | tail -n +2 | awk 'BEGIN {FS = ";" ; OFS=";"} ; {print $1,$2,$3=2*$3}' > movielens_20m.csv

wget -O $DESTDIR/ml-100k.csv http://files.grouplens.org/datasets/movielens/ml-100k/u.data
cut -f 1,2,3 $DESTDIR/ml-100k.csv | tr -s "\t" ";" > movielens_100k.csv

wget -O $DESTDIR/ml-1m.zip http://files.grouplens.org/datasets/movielens/ml-1m.zip
unzip $DESTDIR/ml-1m.zip -d $DESTDIR
tr -s "::" ";" < $DESTDIR/ml-1m/ratings.dat |cut -d ";" -f 1,2,3 > movielens_1m.csv
