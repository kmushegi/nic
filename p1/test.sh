#!/bin/bash

make p > out.txt #p or g to be supplied as a parameter

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