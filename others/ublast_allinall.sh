#!/bin/tcsh

set EVALUE=0.00001
/home/people/thomas/projects/Metagenomics/tools/MOCAT/bin/usearch -ublast $1 -db $2 -evalue $EVALUE -blast6out $1:r.allinall_ublast_nomask


#/home/projects8/pr_53035/people/asli/Enzymes/Koala/allinallblast/HKMN.cdhit.udb
#usearch -makeudb_ublast db.fasta -output db.udb
