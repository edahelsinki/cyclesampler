#!/usr/bin/env python2

import sys

## Build maps
map_user = {}
map_item = {}
count_user = 1
count_item = 1

user = -1
item = -1

data = {}

for row in sys.stdin:
    row = row.strip()
    row = row.split()

    if (len(row) > 0):
        if row[0].startswith( 'review/userId' ):
            user_id = row[1]
            if not user_id in map_user:
                map_user[user_id] = count_user
                count_user += 1
            user = str(map_user[user_id])
        elif row[0].startswith( 'product/productId' ):
            prod_id = row[1]
            if not prod_id in map_item:
                map_item[prod_id] = count_item
                count_item += 1
            item = str(map_item[prod_id])
        elif row[0].startswith( 'review/score' ):
            rating = str(row[1])


    if ( len(row) == 0 ):
        k = (user + "_" + item) 
        if not k in data:
            data[k] = rating
            print ";".join([user, item, rating])

