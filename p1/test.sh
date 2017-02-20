#!/bin/bash

#parameters to script pbil or ga

FILES="maxsat-problems/maxsat-random/highgirth/5SAT/*.cnf"

for f in $FILES
do
	if [ "$1" == "pbil" ]; then
		java EvolAlg $f 100 0.1 0.075 0.02 0.05 1000 p > out.txt
	elif [ "$1" == "ga" ]; then
		java EvolAlg $f 100 bs uc 0.7 0.01 1000 g > out.txt
	else
		echo "Algorithm specified incorrectly. Choose pbil or ga"
		exit 1
	fi
done

if [ ! -f out.txt ]; then
	echo "Output File Not Found!"
	exit 1
fi

if [ ! -f sol.txt ]; then
	echo "Solution File Not Found!"
	exit 1
fi

if cmp -s out.txt sol.txt; then
	echo "The files match"
else
	echo "The files do not match"
fi