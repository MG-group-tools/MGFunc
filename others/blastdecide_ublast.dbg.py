#!/usr/bin/env python2.7
# asli@cbs.dtu.dk, thomas@cbs.dtu.dk
import sys, gzip
import re, string
import argparse
import os
from Bio.Blast import NCBIStandalone
from operator import itemgetter, attrgetter
from datetime import datetime as dt
import time

sys.path.append('/home/projects5/pr_53035/people/asli/bin/lib/python2.7/site-packages')

helpstr = '''
This script parses blast/ublast results and filters them based on the given cut-offs. 
Blast results should be in -m 0 format or tab separated -m 6 format. With ublast, the results should be 
obtained with -blast6out option. 
Author: Asli I. Ozen ,(asli@cbs.dtu.dk)
License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
'''


                

class BlastDecision:
    def __init__(self):
	self.start = time.time()
        d_ = dt.today()
        self.timestarted = d_.strftime("%Y-%m-%d_%H%M")
	self.parseArgs()
	
    def parseArgs(self):
        self.parser = argparse.ArgumentParser(description= helpstr, prog='PROG', conflict_handler='resolve')
	self.parser.add_argument('-bd', help="pre-made blastdecide output file")
	self.parser.add_argument('-bt', help="blast tabular output file")
	self.parser.add_argument('-bo', help="blast -m 0 output file and its format (type 'blast' or 'psiblast')", nargs=2)
	self.parser.add_argument('-l', help="Query lengths file (if tabular blast output is given)")
	self.parser.add_argument('-cb',  type=str, dest="", help="save new blastdecide or not (type 'savenew' or 'nosave') ",nargs=1)
	self.parser.add_argument('-s',  default= "50", dest="", help="minimum similarity cutoff")
	self.parser.add_argument('-qc', default= "50", dest="", help="minimum query coverage cutoff")
	#self.parser.add_argument('-tc',  dest="", help="minimum target coverage cutoff")    	
	self.parser.add_argument('-e',  default= "1e-10" , dest="", help="evalue cutoff")    	
	self.parser.add_argument('-v',  action="store_true" , dest="", help="increase output verbosity")   

  
    def read_lengths(self):
            fl= open(self.lenfile,"rU")
	    self.lendict={}
	    for line in fl:
		   #print line
		   query = line.split("\t")[0]
		   query_name = query.split(" ")[0].strip(">")
		   length= int(line.split("\t")[1].strip("\n"))
		   self.lendict[query_name]=length
	    fl.close()           
                               
    def ReadBlast(self, file, iszipped = 0, is_psiblast=None):
        self.output= open(str(file)+".blastdecide", "w")
	self.selfhits=[]
        if is_psiblast:
            print >> sys.stderr, 'Parsing PSI-Blast'
            self.parser = NCBIStandalone.PSIBlastParser()
        else:
            self.parser = NCBIStandalone.BlastParser()
        if file[-3:] == '.gz' or iszipped:
            handle = gzip.open(file)
        else:
            handle = open(file)
        
        self.iter = NCBIStandalone.Iterator(handle = handle, parser = self.parser)
        self.blastDict = {}
            
        while 1:
            try:
                rec = self.iter.next()
		if not rec: break
            except:
                sys.stderr.write('Can\'t iterate on blast records anymore. Abort.\n')
                import traceback
                traceback.print_exc()
                return 'Error parsing %s' % file
            
            
	   
            self.query = rec.query.split(" ")[0] ##  blast_record.query.split(" ")[0]
            self.length = rec.query_letters
	    
            if self.length < self.min_size:
                if self.verbose: print 'Does not meet the minimum length %d' % self.min_size
                break

            if is_psiblast: rec = rec.rounds[-1]
                        
            # each alignment is one potential hit
            for n, alignment in enumerate(rec.alignments):
                # to make it easy, skip all alignments with multiple HSPS
                
                hsp = alignment.hsps[0]
		alnlength=hsp.align_length
                hit = alignment.title
		#targetlength = alignment.length
		#m = re.search("sp\|([A-Z0-9]+)\|([A-Z0-9_]+) ?(.+)?", alignment.title)

                m = re.search("sp\|(.+?)\|(.+?) (.+)?", alignment.title)
                if m: # pyphynr blast result
                    hit_sp_ac = m.group(1)
                    hit_sp_id = m.group(2)
                    hit_sp_note = m.group(3)
                elif alignment.title[0] == '>': # result from qadditional blast databases
                    hit_sp_ac = None
                    hit_sp_id = alignment.title[1:].split()[0]
                    hit_sp_note = None
                else:
                    hit_sp_ac = None
                    hit_sp_id = None
                    hit_sp_note = None


                # fix annoying dots in ids
                if hit_sp_ac: hit_sp_ac = hit_sp_ac.replace('.','_')
                if hit_sp_id: hit_sp_id = hit_sp_id.replace('.','_')
                

                #if not hit_sp_id: print 'XXXXXXX', alignment.title
                if self.verbose: print hit_sp_id,
                similarity = hsp.positives[0]/float(hsp.positives[1])*100
                if  float(hsp.expect) < int(self.HSP_max_evalue):
		    if  float(similarity) >  int(self.HSP_minimal_positives):
                	coverage = hsp.positives[1]/float(self.length)*100
                	if  float(coverage) > int(self.HSP_minimal_coverage):
                            #targetcoverage = hsp.positives[1]/float(targetlength)*100
                            #if  float(targetcoverage) > int(self.HSP_minimal_targetcov):
                	       #self.compatibles.append((hit_sp_ac, hit))
                	       #hitlist = [hit_sp_id, n+1 , hsp.positives[0]/float(hsp.positives[1])*100, hsp.positives[1]/float(self.length)*100, hsp.positives[1]/float(targetlength)*100, hsp.score, hsp.expect]
                	       hitlist = [hit_sp_id, hsp.positives[0]/float(hsp.positives[1])*100, hsp.positives[1]/float(self.length)*100, hsp.score, hsp.expect]
                               if self.cB: self.createblastDict(query,hitlist)
			       self.output.write("%s\t" % (self.query)),		    
			       for element in hitlist:		    
 	                	   self.output.write("%s\t" % element),	            
			       self.output.write("\n")
        self.output.close()	
	handle.close()
	return None
    

    def ReadBlastresultsTab(self, filename):
        
	if filename[-3:] == '.gz':
	    fh = gzip.open(filename)
	else:
	    fh= open(filename,"rU")

        #hitsdict={}   
	#hitlist = [hit_sp_id, n+1 , hsp.positives[0]/float(hsp.positives[1])*100, hsp.positives[1]/float(self.length)*100, hsp.score, hsp.expect]
	self.blastDict={}
	self.selfhits=[]
	self.read_lengths()	
	output= open(basename(filename) + ".blastdecide", "w")    
        #lines=fh.readlines()
        for line in fh:
	    line = line.strip("\n")
	    if len(line.split("\t")) > 2:
	       query = line.split("\t")[0]
	       #print query						   ###################################################################################
	       query_name = query.split(" ")[0]	      			   # OLD WAY									     #
	       hit_sp_id = line.split("\t")[1]				   # command = "grep -m 1 \"^>" +  query_name + "\" " + "HKMN.cdhit.faa.id_length"   #
	       percent_id = float(line.split("\t")[2])			   # fout2 = os.popen(command)  						     #
	       aln_len=float(line.split("\t")[3])             		   # out2 = fout2.readline().strip("\n")					     #
	       query_length=self.lendict[query_name]	   		   # if out2:									     #
	       coverage = 100*int(aln_len)/float(query_length)  	   #   query_length=out2.split("\t")[1] 					     #
	       bitscore = float(line.split("\t")[11])			   ###################################################################################
	       evalue = float(line.split("\t")[10])
	       if str(query_name) == str(hit_sp_id):
	          #print "SameSameSame"        
                  self.selfhits.append(query)
	       else:
	          if  float(evalue) < int(self.HSP_max_evalue):
		      if  float(percent_id) >  int(self.HSP_minimal_positives):		     
                	  if  float(coverage) > int(self.HSP_minimal_coverage):
				 hitlist=[hit_sp_id, percent_id, coverage, bitscore, evalue]   
				 if self.cB: self.createblastDict(query,hitlist)
				 output.write("%s\t" % (query_name)),		    
				 for element in hitlist:		    
 	                	     output.write("%s\t" % element),	            
				 output.write("\n")
 
	output.close() 					     				     
	fh.close()


    def ReadBlastdecide(self, filename):
        #hitsdict={}   
	#hitlist = [hit_sp_id, n+1 , hsp.positives[0]/float(hsp.positives[1])*100, hsp.positives[1]/float(self.length)*100, hsp.score, hsp.expect]
	fh= open(filename,"rU")	    
        lines=fh.readlines()
	output= open(basename(filename) + ".new.blastdecide", "w") 
        for line in lines:
	    line = line.strip("\n")
	    if len(line.split("\t")) > 2:
	       query = line.split("\t")[0]
	       hit_sp_id = line.split("\t")[1]
	       #n=float(line.split("\t")[2])
	       percent_id = float(line.split("\t")[2])
	       coverage = float(line.split("\t")[3])
	       #targetcoverage = float(line.split("\t")[5])	     
	       bitscore = float(line.split("\t")[4])
	       evalue = float(line.split("\t")[5])
	       if str(query_name) == str(hit_sp_id):
	          #print "SameSameSame"        
                  self.selfhits.append(query)
	       else:
	          if  float(evalue) < int(self.HSP_max_evalue):
		      if  float(percent_id) >  int(self.HSP_minimal_positives):		     
                	  if  float(coverage) > int(self.HSP_minimal_coverage):
				 hitlist=[hit_sp_id, percent_id, coverage, bitscore, evalue]   
				 if self.cB == 'savenew':
				    self.createblastDict(query,hitlist)
				    self.writeoutput(output,query,hitlist)
				 else:
				    self.createblastDict(query,hitlist)

	output.close() 
	
	if self.cB != 'savenew' and os.path.getsize(basename(filename) + ".new.blastdecide") == 0: 
	      os.system("rm " + basename(filename) + ".new.blastdecide")					     				     
	
	fh.close()

    def writeoutputt(self, oh, query, hitlist):
    	
    	oh.write("%s\t" % (query_name))
    	for element in hitlist: 	 
    	    output.write("%s\t" % element),	  
    	oh.write("\n")


    def createblastDict(self, query, hitlist):
    	 self.blastDict={}
    	 self.selfhits=[]
    	 hit_sp_id=hitlist[0]
    	 if str(query) is not str(hit_sp_id):	   
    	    #hitlist=[hit_sp_id, n, percent_id, coverage,targetcoverage, bitscore,evalue]  
    	    #hitlist=[hit_sp_id, percent_id, coverage, bitscore,evalue] 
    	    if query in self.blastDict:
    	       self.blastDict[query].append(hitlist)
    	    else:		     
    	       self.blastDict[query] = [hitlist]				
    	 									 
   										 
    def mainthing(self):							 
    	self.HSP_minimal_positives = self.opts.s  			
    	self.HSP_minimal_coverage = self.opts.qc				 
    	#self.HSP_minimal_targetcov = self.opts.tc    				 
    	self.HSP_minimal_coverage_length = 20					 
    	self.lenfile= self.opts.l						 
    	self.HSP_max_evalue = self.opts.e					 
    	self.verbose = self.opts.v						 
    	self.min_size = 0							 
    	self.cB = self.opts.cb
   
    	if self.opts.bd:
    	    blastdecide=self.opts.bd
    	    self.ReadBlastdecide(blastdecide)

    	elif self.opts.bt:
    	    blasttab =  self.opts.bt
    	    self.ReadBlastresultsTabb(blasttab)
    	else:
    	    try: 
    	       blastfile = self.opts.bo[0]
    	       typ = self.opts.bo[1]
    	       if typ == "psiblast":	   
    		  self.ReadBlast(blastfile,  is_psiblast=True)
    	       else: 
    		  self.ReadBlast(blastfile)
    	    except: 
    	       raise IOError('If you dont have Pre-made blastdecide or ublast-tab results, you should provide a normal blast output (-m0)')


    	timeused = (time.time() - self.start) / 60
    	if self.v: #prints end-time
    	    print "Time used: "+str(round(timeused*60)) + " seconds ("+str(round(timeused)) + " min)\n"


if __name__ == '__main__':

    try:    
        obj = BlastDecision()        
        obj.opts=obj.parser.parse_args(sys.argv[1:])
	obj.mainthing()
	if obj.opts.v: 
	   obj.parser.print_help()
    except Exception,e: 
         print str(e)
   #     print sys.stderr
         sys.exit()

###############
# INPUT LIST
# blast output in tab format & query lengths file : genecatalogue_vs_uniprot.blasttab OR genecatalogue_vs_genecatalogue.blasttab & genecatalogue.lengths
# blast output in -m 0 format : genecatalogue_vs_uniprot.blastout OR genecatalogue_vs_genecatalogue.blastout
# pre-made blastdecide file : genecatalogue_vs_uniprot.blastdecide
# 
# OUTPUT LIST
# new blastdecide file based on given parameters : genecatalogue_vs_uniprot.blastdecide
# if premade blastdecide is given, the blastDict is generated   : obj.blastDict
#
# OPTIONS LIST
# '-bd', '--blastdecide', help="pre-made blastdecide output file"
# '-bt', '--blasttab', help="blast tabular output file"
# '-bo', '--blastout', help="blast -m 0 output file and its type ('blast' or 'psiblast')"
# '-l', '--lengths', help="Query lengths file"
# '-s', '--similarity', default= "50", help="minimum similarity cutoff"
# '-qc', '--querycoverage',default= "50", help="minimum query coverage cutoff"
# '-tc', '--targetcoverage', help="minimum target coverage cutoff"
# '-e', '--maxevalue', default= "1e-10" , help="evalue cutoff"    
#
#
