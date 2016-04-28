#!/bin/bash

for i in $(ls -1 CL_$1*.aln); do

python2.7 ~/tools/thomaspy/bb_phylo.py -c -t paup $i

done
