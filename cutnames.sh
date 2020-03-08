#!/bin/bash
# Extract the header row of  the first argument and save it
head -n 1 $1 | tr ',' '\n' > recs_names.txt
# Search the positions of the columns that the second argument refers to.
cols=$(grep -E -n $2 recs_names.txt | cut -f1 -d: | paste -s -d,)
# Based on the position,
# Extract the columns from the first argument and save them.
cut -f${cols} -d"," $1 > weights.csv
# Return the result
cat weights.csv
