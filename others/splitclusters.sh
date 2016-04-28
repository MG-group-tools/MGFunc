#!/bin/bash

while read line; do

python2.7 /home/projects8/pr_53035/people/asli/Enzymes/Koala/scripts/change_split.py fetch $1 $line.faa $line

done < $2
