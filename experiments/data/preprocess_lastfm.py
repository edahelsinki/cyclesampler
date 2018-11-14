#!/usr/bin/env python3

import sys

def preprocess_lastfm(fname_in, fname_out):
    # Read data once
    with open(fname_in, "r") as f_in:
        data = f_in.readlines()
    del data[0]
    
    user_pair = [0] * len(data)

    data = [i.strip().split("\t") for i in data]

    # It is assumed that the data is sorted by the first
    # column in increasing order.
    k = 0
    for line in data:
        if float(line[0]) < float(line[1]):
            user_pair[k] = ";".join([line[0], line[1], "1"])
            k += 1
    
    user_pair = user_pair[0:k]

    # Write data to file
    with open(fname_out, "w") as f_out:
        f_out.writelines(line + "\n" for line in user_pair)
    
# Run from command line
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: preprocess_lastfm.py <filename_in> <filename_out>")
    
    if len(sys.argv) >= 2:
        preprocess_lastfm(sys.argv[1], sys.argv[2])
