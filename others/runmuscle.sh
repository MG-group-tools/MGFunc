#!/bin/bash

mkdir -p largeseq
for i in $(ls -1 CL_$1*.new.faa); do
count=`grep -c ">" $i`
if [ $count -gt 1000 ] ; then
   if [ $count -gt 1500 ] ; then
      mv $i largeseq
   else

      if [ ! -f $i.aln ] ; then
      muscle -in $i -out $i.mx3.aln -clwstrict -maxiters 3 -quiet -log $i.musclelog
      fi
   fi
   
   
elif [ $count -gt 500 ] ; then

if [ ! -f $i.aln ] ; then
muscle -in $i -out $i.mx8.aln -clwstrict -maxiters 8 -quiet -log $i.musclelog
fi

else
if [ ! -f $i.aln ] ; then
muscle -in $i -out $i.mx16.aln -clwstrict -quiet -log $i.musclelog
fi

fi
done
