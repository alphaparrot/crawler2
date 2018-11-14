#!/bin/bash

#Move output to destination hopper
cp output/*.nc ../output/$1/

#Clean up
rm -rf pickup* output STDOUT*

#run the central release script
python release.py
