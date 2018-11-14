#!/usr/bin/env python2

import sys

## Read data once
def preprocess_tasteprofile(fname_in, fname_out):
    with open(fname_in, "r") as f:
        data = f.readlines()

    ## Build maps
    map_user = {}
    map_item = {}

    count_user = 1
    count_item = 1

    for line in data:
        record = line.split("\t")
        if not record[0] in map_user:
            map_user[record[0]] = count_user
            count_user += 1

        if not record[1] in map_item:
            map_item[record[1]] = count_item
            count_item += 1

    with open(fname_out, "w") as f:
        for line in data:
            record = line.split("\t")
            user   = str(map_user[record[0]])
            item   = str(map_item[record[1]])
            rating = str(record[2])
            f.writelines(";".join([user, item, rating]))

## Run from command line
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: preprocess_tasteprofle.py <filename_in> <filename_out>"

    if len(sys.argv) >= 2:
        preprocess_tasteprofile(sys.argv[1], sys.argv[2])
