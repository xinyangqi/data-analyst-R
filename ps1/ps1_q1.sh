#!/bin/bash
# Create a variable with the name of the csv file. 
file=recs2015_public_v4.csv
# Check if this file exists. If not, download it.
if [ ! -f "${file}" ]
then
	wget https://www.eia.gov/consumption/residential/data/2015/csv/recs2015_public_v4.csv
fi
# Extract the header row of 'recs2015_public_v4.csv' and save it
head -n 1 recs2015_public_v4.csv | tr ',' '\n' > recs_names.txt
# Search positions of id and replicate weight columns and save them
cols=$(grep -E -n "DOEID|BRRWT*" recs_names.txt | cut -f1 -d: | paste -s -d,)
# Use positions to extract id and replicate weight columns
# from the RECS data and save them
cut -f${cols} -d"," recs2015_public_v4.csv > recs_weights.csv
