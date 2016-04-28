#!/usr/bin/env python
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen (asli@cbs.dtu.dk)
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
import sys, gzip
import re, string
import argparse
import os


efile= sys.argv[1]
eh=open(efile,"r")
line= eh.readline()
clid = line.split(" ")[2].strip("\n")
clgos = []
while 1: 

   line= eh.readline()
   if not line: 
      with open(efile + ".new","a") as new:
         for i in clgos:
             new.write(clid + "\t" + i)	    	 
      break
      
   if line[0] =="#":
      with open(efile + ".new","a") as new:
         for i in clgos:
             new.write(clid + "\t" + i)
      clid = line.split(" ")[2].strip("\n") 
      clgos = []  
   elif len(line) == 1: 
      continue
   else:
      clgos.append(line)

eh.close()    
    
