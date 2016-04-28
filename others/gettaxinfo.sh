#!/bin/bas

uniprotids=$1


while read line; do 
echo "$line " >> $1.taxinfo
../getgene_trie.py --get Get_OX $line >> $1.taxinfo
done < $uniprotids


exit
