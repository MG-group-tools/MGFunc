#!/usr/bin/env python
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen (asli@cbs.dtu.dk)
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
import sys, gzip
import re, string
import optparse
import argparse
import os
from Bio.Blast import NCBIStandalone
from operator import itemgetter, attrgetter
sys.path.append('/home/projects5/pr_53035/people/asli/bin/lib/python2.7/site-packages')
import networkx as nx
#import kyotocabinet as kc
from os.path import basename
import time
from datetime import datetime as dt

usage= sys.argv[0] + " [-h] [-c C] [-d D] [-b B] [-bh] [-i] [-v]\n"

helpstr= '''
description: 
This script identifies the uniprot hits for each cluster.
If the option -bh is given, only outputs best hit for a gene in a cluster.\n
First an index file for the genes_vs_database.blastdecide file is created in order 
to fasten the process.
'''

epi="Author: Asli I. Ozen (asli@cbs.dtu.dk)" 

class GetUniprotInfo:
    def __init__(self):
        self.timing=""
	self.dbstarttime = time.time()
        self.parseArgs()
	
    def parseArgs(self):
        self.parser = argparse.ArgumentParser(description=helpstr, epilog=epi)
        #parser.add_option("-p", "--path", dest="P", default="./", type=str, help="path to the blast output files")
        #parser.add_argument("-o", "--opath", dest="OP", default="../parsed/", type="str", help="path to the parse output files")    	
        
        #self.parser.add_argument("-d", "--uniprotdb", type=str, help="uniprot database tab separated format",required=True)
	self.parser.add_argument("-i", metavar="clusterfile", type=str, help="cluster file",nargs=1,required=True)
        self.parser.add_argument("-b", metavar="blastdecide", type=str, help="blastdecide file",required=True)
	self.parser.add_argument("-bh", "--besthit", action = "store_true", help="If chosen only best uniprot hits will be output")
      	self.parser.add_argument("-i", "--index", action = "store_true", help="If chosen, blastdecide file is indexed first.")        
        self.parser.add_argument("-v", "--verbose", action="store_true" , help="increase output verbosity")
        return 1

	
    def AlignClustal(self, infile, outfile, options = '-quicktree', cleaner = None, do_remove_replicates = None):
        "method designed to do the actual alignment"
        id = stem(infile)
        
        com = "%s -infile=%s -outfile=%s %s" % (self.clustal, infile, outfile, options)
        #com = "%s %s %s > %s " % (self.clustal, infile, options, outfile)
        Message('\t' + com)
        os.system(com)

        
        if cleaner:
            Message('Cleaning alignment')
            cleaner.Clean(outfile, id)
            Message('Done cleaning')
        else:
            Message('skipping cleaning alignment part')

        if not exists(outfile) and exists(outfile +'.gz'):
            outfile = ifzipped_unzip(outfile +'.gz')
            
        if do_remove_replicates and myexists(outfile):
            from aln import AlnReader
            Message('Removing replicates from alignment')
            a = AlnReader(outfile)
            a.db_index = self.db_index
            a.remove_replicates(outfile, threshold = do_remove_replicates,
                                self_reference = self.config.get('self_reference'),keep=id)
				
				
				
    def AlignMuscle(self, options):
        "method designed to do the actual alignment with Muscle"
        base = ".".join(infile.split(".")[0:-1])
        
	while 1: 
        
	
	     try: 
	     
	        infile = self.jobs.get(True,1)
		outfile = base + ".aln"
	        com = "muscle -in %s -out %s -clwstrict %s" % (infile, outfile, options)
                os.system(com)

        
    #    if cleaner:
#            Message('Cleaning alignment')
#            cleaner.Clean(outfile, id)
#            Message('Done cleaning')
#        else:
#            Message('skipping cleaning alignment part')
#
#        if not exists(outfile) and exists(outfile +'.gz'):
#            outfile = ifzipped_unzip(outfile +'.gz')
#            
#        if do_remove_replicates and myexists(outfile):
#            from aln import AlnReader
#            Message('Removing replicates from alignment')
#            a = AlnReader(outfile)
#            a.db_index = self.db_index
#            a.remove_replicates(outfile, threshold = do_remove_replicates,
#                                self_reference = self.config.get('self_reference'),keep=id)

    def AlignMafft(self, infile, outfile, options = '', cleaner = None, do_remove_replicates = None):
        "method designed to do the actual alignment with Mafft"
        id = stem(infile)
        options = options.replace('-quicktree','')
        if not options:
            try:
                options = self.muscle_options
            except:
                options = ''

        
        com = "cat %s | readseq -p| mafft --clustalout - | sed 's/CLUSTAL (-like) formatted alignment by MAFFT FFT-NS-2 (v6.240)/CLUSTAL W (1.81) multiple sequence alignment/g' > %s" % (infile, outfile)
        #com = "%s %s %s > %s " % (self.clustal, infile, options, outfile)
        Message('\t' + com)
        os.system(com)

        
        if cleaner:
            Message('Cleaning alignment')
            cleaner.Clean(outfile, id)
            Message('Done cleaning')
        else:
            Message('skipping cleaning alignment part')

        if not exists(outfile) and exists(outfile +'.gz'):
            outfile = ifzipped_unzip(outfile +'.gz')
            
        if do_remove_replicates and myexists(outfile):
            from aln import AlnReader
            Message('Removing replicates from alignment')
            a = AlnReader(outfile)
            a.db_index = self.db_index
            a.remove_replicates(outfile, threshold = do_remove_replicates,
                                self_reference = self.config.get('self_reference'),keep=id)

    
    
    
    def UseMafft(self): self.Align = self.AlignMafft
    
    def UseMuscle(self): self.target = self.AlignMuscle
    
    def UseMuscleQuick(self):
        self.muscle_options = ' -diags1'
        self.Align = self.AlignMuscle
    def UseMuscleQuickest(self):
        self.muscle_options = ' -maxiters 1 -diags1 -sv -distance1 kbit20_3'
        self.Align = self.AlignMuscle



   





	self.getuniprotresults()
        self.timing= (time.time() - self.dbstarttime) /60
	if self.opts.verbose: #prints end-time
              print "### Time used: "+str(round(self.timing*60)) + " seconds ("+str(round(self.timing)) + " min)\n"


if __name__ == '__main__':


