#!/bin/bash

while read line; do

python2.7 /home/projects8/pr_53035/people/asli/Enzymes/Koala/scripts/change_split.py changeheader $1/$line.faa $1/$line.new.faa $1/$line.headertable

done < $2
