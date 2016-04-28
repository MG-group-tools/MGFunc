#!/usr/bin/env python

import os
import sys
import re
from Bio import SeqIO
from datetime import datetime as dt
import time
import argparse
helpstr = '''
description: This is a small utility script for splitting of sequences for separate clusters. Clusters.fasta file must be already processed and the headers should be
accordingly. 
'''

epi="Author: Asli I. Ozen (asli@cbs.dtu.dk)" 


class Split:
    
      def __init__(self):
          self.start = time.time()							
	  d_ = dt.today()
	  self.timestarted = d_.strftime("%d-%m-%Y %H:%M:%S")
	  self.parseArgs()										  

      def printer(self, string):
          if self.opts.verbose: 
	     print string
      def parseArgs(self):	
	  self.parser = argparse.ArgumentParser(description=helpstr, epilog= epi)
	  self.parser.add_argument("-i", metavar ="infasta", type=str, help="Input FASTA FILE (required)",nargs=1,required=True)
	  #self.parser.add_argument("-n", metavar ="newfasta",type=str, help="output FASTA FILE extention name (default CL_ID.fsa)",nargs=1)
	  self.parser.add_argument("-cl", metavar = "clusters", type=str, help="input list of cluster ids to be separated from the fasta file of all clusters (required for -s option)",nargs=1)
	  #self.parser.add_argument("-t", metavar = "headertable", type=str, help="output header conversion table FILE name",nargs=1)
	  self.parser.add_argument("-o", metavar ="ofold", type=str, help="Output folder name")
	  self.parser.add_argument("-v", "--verbose", action="store_true" , help="increase output verbosity, prints details, warnings etc.")
	  return 1

 
      def splitclusters(self,fsafile, newfasta, cl_id):

	  handle = open(fsafile, "rU") 
	  output_handle = open(newfasta, "w")
	  for record in SeqIO.parse(handle, "fasta"):     
              header = record.description
	      clusterid=header.split(":")[0]
	      geneid=header.split(":")[1]
	      if clusterid == cl_id:
		 SeqIO.write(record, output_handle, "fasta")

	  handle.close()
	  output_handle.close() 
      
      def mainthing(self):
      
          #self.opts=self.parser.parse_args(sysargs)
	  self.fsa = self.opts.i[0]
	  self.outfold = self.opts.o
          
	  clidlist= open(self.opts.cl[0], "r").readlines()
	  for line in clidlist :
	      clfsa = self.outfold + line.strip("\n") + ".fsa"
	      self.splitclusters(self.fsa, clfsa, line.strip("\n"))
  
          timeused = (time.time() - self.start) / 60
          self.printer("### Time used: "+str(round(timeused*60)) + " seconds ("+ str(float(timeused)) + " min)\n")

if __name__ == '__main__':
	
     try:    
         obj = Split()        
         obj.opts=obj.parser.parse_args(sys.argv[1:])  # test: ["-s", "-i", "clusters_size4.semi.koalagenes.fasta", "-cl", "test.clid", "-v"]
         if obj.opts.verbose: 
            print "\n### " + sys.argv[0] + " initialized at " + obj.timestarted  
            print "### OPTIONS : " + str(obj.opts)
         obj.mainthing()

          #  obj.parser.print_help()
     except Exception,e: 
         print str(e)
         import traceback
         traceback.print_exc()
	
#	
#	except  StandardError, err:
#		print "there is an error" , err
#	except:
#		print "An unknown error occurred.\n"    
#      
#      
#
#

    
