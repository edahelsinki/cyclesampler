#!/usr/bin/env python2

import sys

def preprocess_bxbooks(fname_in, fname_out):
    ## Read data once
    with open(fname_in, "r") as f_in:
        data = f_in.readlines()
    del data[0]
    
    ## Build maps
    map_user = {}
    map_book = {}

    count_user = 1
    count_book = 1

    for line in data:
        record = line.split(";")
        record[1] = record[1].replace(" ", "")
        record[1] = record[1].replace('\\"', "")
        record[1] = record[1].replace('=', "")
        record = [i.replace('"', '').strip() for i in record]
        
        if not record[0] in map_user:
            map_user[record[0]] = count_user
            count_user += 1

        if not record[1] in map_book:
            map_book[record[1]] = count_book
            count_book += 1

    ## Write data to file
    with open(fname_out, "w") as f_out:
        for line in data:
            record = line.split(";")
            record[1] = record[1].replace(" ", "")
            record[1] = record[1].replace('\\"', "")
            record[1] = record[1].replace('=', "")
            record = [i.replace('"', '').strip() for i in record]

            user   = str(map_user[record[0]])
            book   = str(map_book[record[1]])
            rating = str(record[2]).replace('"', '')
            out = ";".join([user,  book, rating+"\n"])
            f_out.writelines(out)

## Run from command line
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: preprocess_bxbooks.py <filename_in> <filename_out>"

    if len(sys.argv) >= 2:
        preprocess_bxbooks(sys.argv[1], sys.argv[2])
