#!/usr/bin/python
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen (asli@cbs.dtu.dk)\n
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
import os
import sys,gzip,zlib
#import re 
import subprocess
import glob
from Bio import SeqIO
from optparse import OptionParser
import argparse
import time
from datetime import datetime as dt
from os.path import basename
import re
import pickle


helpstr = '''
description: This script creates an index file for each sequence in a fasta file,
where the header of the sequence, start and end position indexes
of the sequence are printed out in a tab separeted format and python 
dictionary object is saved. This script is designed to help fasten the 
process when working with large datasets. 
'''

epi="Author: Asli I. Ozen (asli@cbs.dtu.dk)" 


class Fasta2Index:

    def __init__(self):
	self.timing=""
	self.starttime = time.time()
	d_ = dt.today()
        self.timestarted = d_.strftime("%d-%m-%Y %H:%M:%S")
        self.parseArgs()

    def parseArgs(self):
        self.parser = argparse.ArgumentParser(description= helpstr, epilog=epi)
	self.parser.add_argument("-p", metavar="path", default="./", type=str, help="path to the fasta sequence files")
	self.parser.add_argument("-f", metavar="fastafile",type=str, help="fasta file of the genome",required=True)
	self.parser.add_argument("-o", metavar="outputname", type=str, help="name to the output index files")
	self.parser.add_argument("-op", metavar="outpath", default="./", type=str, help="path to the output index files")
	self.parser.add_argument("-v", "--verbose", action="store_true" , help="increase output verbosity")

	#parser.add_option('-n', '--ncpu', dest="ncpu", default=10, type="int", help="Number of CPU's to use.")
        return 1
	
	
    def createDBindex(self):
        #path to fasta files
	#self.path = opts.P
	self.fasta = self.opts.f
	if self.opts.o:
	   self.fastaout = self.opts.o
	else: 	
	   self.fastaout = basename(self.fasta).split(".")[0]
	if self.opts.op[-1] == "/":
	   self.opath = self.opts.op
	else:
	   self.opath = self.opts.op + "/"
	
 
	if os.path.exists(self.opath + self.fastaout + ".db_index.p"): 
		self.db_index= pickle.load( open(self.opath + self.fastaout + ".db_index.p", "rb" ) )
	else:
	        self.db_index = {}
	        fid = open(self.fasta)
	         
	        indexfile = self.fastaout + ".indexed"
	  
	        ind = open(self.opath + indexfile, "w")
		lng = open(self.opath + self.fastaout + ".length", "w")
	        line = fid.readline()
		
	        header= line.strip("\n").split(" ")[0]                #### NOTE: assume anything after space is irrelevant
                start = fid.tell()
		
		try: 
		   sampleid=header.split("_")[0]
		except Exception,e: 
	            print str(e) , "Your gene catalog headers are wrong!! It should look like this: >SampleID_geneid. \n Yours look like this: " , header , "\n"
		    sys.exit()
	   
                seqlength= 0
		
	        while 1:
	           line = fid.readline()
		   
	           #print line.strip("\n"), fid.tell()
	           
		   if not line or not len(line):
	              stop = fid.tell()
	              self.db_index[header[1:]] = [start, stop]
	              ind.write("%s\t%d\t%d\n" %(header[1:], start, stop))
		      lng.write("%s\t%d\n" %(header[1:], seqlength)) 
	              break

	           if line[0] == ">":	    
                	stop = fid.tell() - len(line) 
                	self.db_index[header[1:]] = [start, stop]
                	ind.write("%s\t%d\t%d\n" %(header[1:], start, stop)) 
			lng.write("%s\t%d\n" %(header[1:], seqlength))
                	header= line.strip("\n").split(" ")[0]
                	start = fid.tell()
			seqlength= 0
	           else:
		     linelen = len(line.strip("\n"))
		     seqlength += linelen

	        lng.close()
	        fid.close()   
                pdbindex= open(self.opath + self.fastaout + ".db_index.p", "wb" )
	        pickle.dump(self.db_index, pdbindex)
	        pdbindex.close()

        self.timing= (time.time() - self.starttime) /60
	if self.opts.verbose: #prints end-time
              print "### Time used for running: "+str(round(self.timing*60)) + " seconds ("+str(round(self.timing)) + " min)"
	      timeended= dt.today().strftime("%d-%m-%Y %H:%M:%S")
	      #print "### " + sys.argv[0] + " finished at :", timeended
	      
	      
if __name__ == '__main__':
     
     try:        
        obj=Fasta2Index()		
	obj.opts = obj.parser.parse_args(sys.argv[1:])
	if obj.opts.verbose: 
	   print "\n### " + sys.argv[0] + " initialized at " + obj.timestarted  
	   print "### OPTIONS : " + str(obj.opts)	
	obj.createDBindex()
		
     #except IOError as i:
	# print "I/O error({0}): {1}".format(i.errno, i.strerror)   
     except Exception,e: 
	print str(e)
        import traceback
        traceback.print_exc()
###############
# INPUT LIST
# Input Fasta file to be indexed
# 
#
# OUTPUT LIST
# tab separated list, column 1 fasta header, column2-3 start and end positionindex for sequence 
# outputname: fastafile.indexed
#
#
#
# TODO: 
# Add regex for spaces or odd characters in headers and add warnings
#
#
#
#
#
#
#
#
#
#
