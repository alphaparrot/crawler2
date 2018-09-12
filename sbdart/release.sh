#!/bin/bash

#Move output to destination hopper
cp sbout.* $1/

#Clean up
rm -rf sbout.*
rm -rf sbdart-*

