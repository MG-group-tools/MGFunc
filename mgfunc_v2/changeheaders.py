#!/usr/bin/env python

import os
import sys
import re
from Bio import SeqIO
from datetime import datetime as dt
import time
import argparse
helpstr = '''
description: This is a small utility script for changing of headers of the fasta files to unique and short names for the consequent alignment. 
changeheader can be used with any fasta file, separately.
'''

epi="Author: Asli I. Ozen (asli@cbs.dtu.dk)" 


class Change:
    
      def __init__(self):
          self.start = time.time()							
	  d_ = dt.today()
	  self.timestarted = d_.strftime("%d-%m-%Y %H:%M:%S")
	  self.parseArgs()										  
      def printer(self,string): 
          if self.opts.verbose:
             print string


      def parseArgs(self):	
	  self.parser = argparse.ArgumentParser(description=helpstr, epilog= epi)
	  self.parser.add_argument("-i", metavar ="infasta", type=str, help="Input FASTA FILE (required)",nargs=1,required=True)
	  self.parser.add_argument("-n", metavar ="newfasta",type=str, help="output name for new FASTA FILE (default = inputfile.new.fsa)",nargs=1)
	  self.parser.add_argument("-ht", metavar = "headertable", type=str, help="output header conversion table FILE name (default = inputfile.headertable)",nargs=1)
	  self.parser.add_argument("-v", "--verbose", action="store_true" , help="increase output verbosity, prints details, warnings etc.")
	  return 1

      def changeheader(self,fsafile, newfasta, headertable):

	  handle = open(fsafile, "rU") 
	  output_handle1 = open(newfasta, "w")
	  output_handle2 = open(headertable, "w")
	  names = {}
	  a = 0
	  for record in SeqIO.parse(handle, "fasta"):
	     a += 1
	     if a < 100000:
        	 newheader = 'Seq' + str(a)

	     if a >= 100000 and a < 1000000 :
		 newheader = 'Seq' + str(a)

	     if a >= 1000000:
		 newheader = 'M' + str(a)

	     output_handle2.write(newheader + "\t" + record.description + "\n")
	     record.description = ''
	     record.name = ''     
	     record.id = newheader
	     SeqIO.write(record, output_handle1, "fasta")



	  handle.close()
	  output_handle1.close()
	  output_handle2.close()


       
      def mainthing(self):
      
          #self.opts=self.parser.parse_args(sysargs)
	  self.fsa = self.opts.i[0]
	  if self.opts.n: 
	     self.newfsa = self.opts.n[0]
	  else: 
	     tmp = self.fsa.split(".")[0:-1]
	     self.newfsa = ".".join(tmp) + ".new.fsa"
	     #print self.newfsa
	  
	  if self.opts.ht: 
	     self.ht = self.opts.ht[0]
	  else: 
	     tmp = self.fsa.split(".")[0:-1]
	     self.ht = ".".join(tmp) + ".headertable"
	  
	  self.changeheader(self.fsa, self.newfsa, self.ht)
	  
          timeused = (time.time() - self.start) / 60
          if self.opts.verbose: #prints end-time
             self.printer("### Time used: "+str(round(timeused*60)) + " seconds ("+str(round(timeused)) + " min)\n")

if __name__ == '__main__':
	
     try:    
         obj = Change()        
         obj.opts=obj.parser.parse_args(sys.argv[1:])  # test: ["-s", "-i", "clusters_size4.semi.koalagenes.fasta", "-cl", "test.clid", "-v"]
         if obj.opts.verbose: 
            obj.printer("\n### " + sys.argv[0] + " initialized at " + obj.timestarted  )
            obj.printer("### OPTIONS : " + str(obj.opts))
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

    
