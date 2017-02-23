#!/bin/bash

#parameters to script pbil or ga

FILES="sample-problems/*.cnf"

for f in $FILES
do
	if [ "$1" == "pbil" ]; then
		java EvolAlg $f 100 0.1 0.075 0.02 0.05 1000 p >> out.txt
		echo -e -n "\n" >> out.txt
	elif [ "$1" == "ga" ]; then
		java EvolAlg $f 100 bs uc 0.7 0.01 1000 g >> out.txt
		echo -e -n "\n" >> out.txt
	else
		echo "Algorithm specified incorrectly. Choose pbil or ga"
		exit 1
	fi
done