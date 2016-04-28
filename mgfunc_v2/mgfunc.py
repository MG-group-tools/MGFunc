#!/usr/bin/env python
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen (asli@cbs.dtu.dk)
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
import sys, gzip
import re, string
import argparse
import os
from operator import itemgetter, attrgetter
import pickle
#sys.path.append('/home/projects5/pr_53035/people/asli/bin/lib/python2.7/site-packages')
#import networkx as nx
from os.path import basename
import commands, tempfile
import time
from datetime import datetime as dt
import fasta2index
import traceback
import glob
import Queue
import threading


#import Utils
#import swiss2fasta
#import swiss2tab
#import bldecide
#import clustergenes
#import clusterFilter
#import parseclusterinfo
#import cluster2fasta
#import splitclusters
#import changeneaders
#import runtrees
#import treenames2

usage = sys.argv[0] + " [-h] -c <configfile> []" 
helpstr = '''
Description: MGFunc is a functional and taxonomical annotation tool for large metagenomic data sets.
'''

epi="Author: Asli I. Ozen (asli@cbs.dtu.dk), Kosai Iskold (kosai@cbs.dtu.dk)" 


class MGFunc:

      def __init__(self): 
          self.start = time.time()
          d_ = dt.today()
          self.timestarted = d_.strftime("%d-%m-%Y %H:%M:%S")
	  self.timestamp = "_".join(self.timestarted.split(" "))
          self.parseArgs()	
	  self.inputo = {}  
	  with open(self.timestamp + ".runlog" , "w") as runlog:
	     runlog.write("RUN LOG " + self.timestamp)
	  
	  with open(self.timestamp + ".runerr" , "w") as runerr:
	     runerr.write("ERROR LOG " + self.timestamp)  
	  self.SECTIONS = ['beginning','formatswiss','ublast_uniprot','ublast_genecatalog','blastdecide_uniprot','blastdecide_genecatalog',
	      'clustergenes','clusterFilter','parseclusterinfo',
              'enrichment', 'cluster2fasta', 'splitclusters', 'runalignment','trees', 'treenames']
             
      def printer(self,string): 
          if self.opts.verbose:
             with open(self.timestamp + ".runlog", "a") as myfile:
	          myfile.write(string)
		  myfile.write("\n")

      def printer2(self,string): 
          if self.opts.verbose:
             with open(self.timestamp + ".runstdout", "a") as myfile:
	          myfile.write(string)	     
      def invokeexit(self,sectionname,specificerror):
          with open(self.timestamp + ".runerr", "a") as exfile:
               exfile.write(specificerror  + "\n")
               exfile.write("Problem with cofiguration in section:  " + sectionname + "\n")
	       exfile.write("Execution Halted!\n")
	  sys.exit()
	  
      def printerr(self,specificerror):
          with open(self.timestamp + ".runerr", "a") as errfile:
	       errfile.write("Exception error \n")
               errfile.write(specificerror  + "\n")

      def makeoptslog(self):
      
          optslog = open(self.timestamp + ".optslog", "w")
	  optslog.write("ALL OPTIONS USED IN THIS RUN: %s \n" % self.timestarted)
          self.P_obj.makedict()
	  a=0
	  for section in self.SECTIONS:
	      
	      optslog.write("###### SECTION %s: %s #######\n" % (str(a),section))
	      for varlist in self.P_obj.d[section]:	          	          
	          if len(varlist)>1:
		    optslog.write("%s %s \n" % (varlist[0],varlist[1]))
		  #else:
 		    #optslog.write("%s : True\n" % varlist[0])
              optslog.write("\n")
	      a+=1
	  optslog.close()	 
	           
#####################################################################################################################################
####    															 ####
####    						    INITIATE								 ####
####    															 ####
#####################################################################################################################################

     
      def parseArgs(self):
          self.parser = argparse.ArgumentParser(description= helpstr, epilog=epi)
	  self.parser.add_argument("-g", metavar="genecatalog",type=str, help="Fasta file of the gene catalog") #KOSAI# removed required=True
	  self.parser.add_argument("-gl", metavar="gclist",type=str, help="A file with list of fasta files of the gene catalog, each file represents a sample,\
	                                                                   must be in the same folder with the gene catalog")	
	  self.parser.add_argument("-d", metavar="datfile", type=str, help="Specify name of uniprot.dat file. Include path, if it is not in current directory.", required = True)	                                                                     
	  self.parser.add_argument("-df", metavar="datfasta", type=str, help= "Specify name of uniprot.fasta file. Include path, if it is not in current directory.\
	                                                                      If it is not found and names are not specified in configuration, \
								              swiss2fasta.py will run automatically." )		  
	  self.parser.add_argument("-dt", metavar="dattab", type=str, help= "Specify name of uniprot.dat.tab file. Include path, if it is not in current directory.\
	                                                                      If it is not found and names are not specified in configuration, \
								              swiss2tab.py will run automatically." )		  
	  self.parser.add_argument("-c", metavar="config", default="MGFunc.ini", type=str, help="configuration FILE. Make sure that the file extention is .ini")
	  self.parser.add_argument("-p", metavar="procs", default=3, type=int, help="maximum processors")

	  self.parser.add_argument("-v", "--verbose", action="store_true" , help="increase output verbosity")

      def parseconfig(self):
          if self.opts.c:
	    if os.path.exists(self.opts.c):
	       self.configfile = self.opts.c
	    else: 	     
	       serr="WHALE WHALE WHALE..Looks like someone is missing the config file. "
	       self.invokeexit("parseconfig", serr)
          
          from Utils import ParseConfig
	  self.P_obj =  ParseConfig(self.configfile)
	  self.P_obj.parser()
	  self.config = self.P_obj.cparser
	  self.p = self.opts.p
	  self.printer("### Config dictionary created\n")




#####################################################################################################################################
####																 ####
####					     INPUT  CHECK/PREPARE              						         #### DONE
####																 ####
#####################################################################################################################################

      def startmain(self):
          from Utils import StartMain
          start = StartMain(self.genecatalog, self.workdir)
	  self.startdir = start.workdir
	  start.prepareinput()
 
 
      def genecatcheck(self):  
          if self.opts.g: 
	     self.typ = "single"              
	     self.genecatalog = self.opts.g
	     self.filebase = basename(self.genecatalog).split(".")[0]
	     self.genecatpath= "/".join(os.path.abspath(self.genecatalog).split("/")[:-1]) + "/"
	     return self.typ	     
	  if self.opts.gl:
	     self.typ = "list"
	     self.genecatalog=[]
	     self.genecatpath= "/".join(os.path.abspath(self.opts.gl).split("/")[:-1]) + "/"
	     self.filebase = basename(self.opts.gl).split(".")[0]
	     
	     for line in file(self.opts.gl):	         
	         self.genecatalog.append(self.genecatpath + line.strip("\n"))	#KOSAI# Added added line.strip()+	 
	     #print self.genecatalog
	     return self.typ
	  else: 
	     msg = " How about an input Gene catalog file? Use options -g or -gl"
	     self.invokeexit("Input options -g/-gl", msg)
	  	 
#####################################################################################################################################
####    															 ####
####    						    BEGINNING								 #### DONE
####    															 ####
#####################################################################################################################################      

      
      def beginning(self):
	  self.s0 = 'beginning'
	  self.inputo[self.s0] = {}
          self.printer(self.s0)
	  if self.config.has_section(self.s0): 
                  if self.config.get(self.s0, 'workdir')[-1] == "/":
	             self.workdir = self.config.get(self.s0, 'workdir')  	     
	          else:
	             self.workdir = self.config.get(self.s0, 'workdir') + str("/")
	          os.system("mkdir -p " + self.workdir)
	          
	          
	          if self.config.get(self.s0, 'usearchpath')[-1] == "/":
	             self.usearchpath = self.config.get(self.s0, 'usearchpath')
	          else:
	             self.usearchpath = self.config.get(self.s0, 'usearchpath') + str("/")
	          
	          if self.config.get(self.s0, 'pkgdir')[-1] == "/":
	             self.pkgdir = self.config.get(self.s0, 'pkgdir')
	          else:
	             self.pkgdir = self.config.get(self.s0, 'pkgdir') + str("/")
		  self.godb = self.config.get(self.s0, 'gene_ontology_obo')
		  if os.path.exists("./" + basename(self.godb)): 
		     self.printer("Found obo file")
		  else:
		     os.system("ln -s " + self.godb + " " + basename(self.godb)) 
		  
		  if self.config.has_option(self.s0, 'preclusterfile') : self.preclusterfile = self.config.get(self.s0, 'preclusterfile') 
		  
		       
	  else: 
	     msg = "Make sure there is a [beginning] section in the .ini file"
	     self.invokeexit(self.s0, msg)

             
#####################################################################################################################################
####    															 ####
####    						    FORMAT UNIPROT DATABASE FILES					 ####   DONE
####    															 ####
#####################################################################################################################################      


      def SwissFasta(self,optslist): 
          import swiss2fasta
	  try:
	     SF_obj = swiss2fasta.main()
	     SF_obj.printer = self.printer
	     self.printer("\n### swiss2fasta initialized at "+ SF_obj.timestarted + "\n")
	     self.printer("### OPTIONS: "+str(SF_obj.args)+"\n")
	     SF_obj.mainthing()
	  except IOError as i:
	  	print "I/O error({0}): {1}".format(i.errno, i.strerror)
	  except Exception,e:
	  	print str(e)
          	traceback.print_exc()
          

      def SwissParse(self,optslist):
          import swiss2tab
	  try:
              SP_obj = swiss2tab.main()
              SP_obj.args = SP_obj.parser.parse_args(optslist)
	      SP_obj.printer = self.printer
              self.printer("\n### swiss2tab initialized at "+ SP_obj.timestarted + "\n")
              self.printer("### OPTIONS: "+str(SP_obj.args)+"\n")
              SP_obj.mainthing()
	  except IOError as i:
              print "I/O error({0}): {1}".format(i.errno, i.strerror)
	  except Exception,e:
              print str(e)
              traceback.print_exc()
              sys.exit()
	 ##################################### UNIPROT INPUT CHECK ############################################

      def uniprotcheck(self):
          self.s1 = 'formatswiss'
	  self.inputo[self.s1] = {}
          '''
	  Function checks if -makefasta and -maketab options are True/False. 
	  If False, function first checks the location of the uniprot.dat file for .tab and .fasta extensions 
	  of the Uniprot file. If the files are missing, and the options are False, gives an error and exits. 
	  '''
	  self.printer("Started Uniprot check")
	  self.convert2fasta = False
	  self.convert2tab = False
	  self.uniprotdat = self.opts.d
	  if basename(self.uniprotdat)[-3:] == ".gz":	     
	     ubase = basename(self.uniprotdat)[:-3]
	  else: 
	     ubase = basename(self.uniprotdat)
	  testloc = "/".join(os.path.abspath(self.uniprotdat).split("/")[:-1]) + "/"
	  
	  
	  
	  if self.opts.dt:
	     self.uniprottab = self.opts.dt
	  else: 	     
	     if os.path.exists(testloc + ubase + ".tab"):
	        self.uniprottab = testloc + ubase + ".tab"
		self.printer("Tab file found: " + self.uniprottab)
	     else: 
	        
	        if self.config.get(self.s1, '-maketab') == 'False':
		   ermsg = "NO Uniprot tab file found and Maketab option is False. Check input or config\n Expected file:" + testloc + ubase + ".tab\n"		   
		   self.invokeexit(self.s1,ermsg) 
		   
		if self.config.get(self.s1, '-maketab') == 'True':
	           self.convert2tab = True
          
	  
          if os.path.exists(testloc + ubase + ".tab.indexed"):
	        self.uniprottabindex = testloc + ubase + ".tab.indexed"
		self.printer("Tab Index file found: " + self.uniprottab)
          else: 
	        ermsg ="Tab index not found, run swiss2tab and keep the index file(.indexed) in the same location with the tab file\n" 
		errmsg2 = "Expected file: " + testloc + ubase + ".tab.indexed\n"
	        self.invokeexit(self.s1,ermsg + errmsg2)
		
		
	  if self.opts.df:
	     self.uniprotfasta = self.opts.df
	  else: 	     
	     if os.path.exists(testloc + ubase + ".fasta"):
	        self.uniprotfasta = testloc + ubase + ".fasta"
		self.printer("Fasta file found: " + self.uniprotfasta)

	     else: 
	        if self.config.get(self.s1, '-makefasta') == 'False':
		   ermsg = "NO Uniprot FASTA file found and Makefasta option is False. Check input or config\n"
		   self.invokeexit(self.s1,ermsg)
		   
		if self.config.get(self.s1, '-makefasta') == 'True':
	           self.convert2fasta = True

	  if os.path.exists(testloc + ubase + ".fasta.indexed"):
               self.printer("Fasta index found: " + self.uniprotfasta + ".indexed")
	  else:
               ermsg ="Fasta index not found, run fasta2index and keep the index file(.indexed) in the same location with the fasta file\n" 
	       self.invokeexit(self.s1,ermsg)

      def formatswiss(self):         
	  
	  self.printer(self.s1)
	  if self.config.has_section(self.s1): 
             if self.config.get(self.s1, '-i') == 'DEFAULT':	
	     
                   ###################UNIPROT FASTA OPTIONS ##########################################	         
      
		   if self.convert2fasta:
		      if self.config.get(self.s1, '-o') == 'DEFAULT':
			 self.uniprotfasta = ".".join(self.uniprotdat.split(".")[:-1]) + ".fasta"    #uniprottabbase.fasta
			 self.P_obj.updateConfig(self.s1, '-i', self.uniprotdat)
			 self.P_obj.updateConfig(self.s1, '-o', self.uniprotfasta)
	              else: 
		         self.uniprotfasta = self.config.get(self.s1, '-o') + ".fasta"
			 self.P_obj.updateConfig(self.s1, '-i', self.uniprotdat)
			 self.P_obj.updateConfig(self.s1, '-o', self.uniprotfasta)
			 
		      s1fopts = self.P_obj.converttolist(self.s1)
		      
		      try:		
			 self.SwissFasta(s1fopts)  # self.GC_vs_UN_ublast	
		      except Exception,e:
			 print str(e)
			 traceback.print_exc()	   

                   ###################UNIPROT TAB OPTIONS ##########################################   # SHOULD WE INDEX TAB file?
	           if self.convert2tab: 
		      if self.config.get(self.s1, '-o') == 'DEFAULT':
			 self.uniprottab = self.uniprotdat + ".tab"
			 self.P_obj.updateConfig(self.s1, '-i', self.uniprotdat)
			 self.P_obj.updateConfig(self.s1, '-o', self.uniprottab )

                      else: 
			 self.uniprottab = self.config.get(self.s1, '-o') + ".tab"			 
			 self.P_obj.updateConfig(self.s1, '-i', self.uniprotdat)
			 self.P_obj.updateConfig(self.s1, '-o', self.uniprottab)
		      
		      
		      s1topts = self.P_obj.converttolist(self.s1)
                      
		      ############## RUN FUNCTION ####################
		      if self.config.get(self.s7,'-run').upper() == "TRUE":
			 try:		
			    self.SwissParse(s1topts)  
			 except Exception,e:
			    print str(e)
			    traceback.print_exc()	   
 
      
#####################################################################################################################################
####    															 ####
####    						    RUN UBLAST 								 ####  DONE
####    															 ####
#####################################################################################################################################      
    
      def ublast(self, fasta, db, out):
          #print fasta, db, out
	  #self.UBevalue = 0.00001
	  ublastcmd=self.usearchpath + "usearch -ublast "+ fasta +  " -db " + db + " -evalue " + str(self.ubevalue) + " -blast6out " + out
	  if self.uprocs: 	             
	     ublastcmd  =  ublastcmd  + " -threads " + self.uprocs 
	  
	  os.system(ublastcmd)
          
	  return 1
	  
      def ublast_makedb(self, fasta, udbname):
             udbcmd= self.usearchpath + "usearch -makeudb_ublast " + fasta + " -output " + udbname 		
	     os.system(udbcmd)
             return 1




      def ublast_genes_vs_uniprot(self, samplefasta, savedir, savedb):
          #self.printer("Running ublast genes catalog against uniprot")	
	  self.uniprotdb_basename = ".".join(basename(self.uniprotfasta).split(".")[0:-1])  
	  if self.makedatabase_u: 	     
	     self.uniprotdb_udb = savedb  + self.uniprotdb_basename + ".udb"  # "/udb_files
	     #print self.uniprotdb_udb
	     self.ublast_makedb(self.uniprotfasta, self.uniprotdb_udb)

	  filebase = basename(samplefasta).split(".")[0]
	  GC_vs_UN_ublast = savedir + filebase + ".genesvsdb.ublast"  # NAMES could be shortened / results saved to /blast instead
	  self.ublast(samplefasta, self.uniprotdb_udb, GC_vs_UN_ublast)
	  self.GC_vs_UN[samplefasta] = GC_vs_UN_ublast
          return 1
	  
      def ublast_uniprot(self):      
	  self.s2 = 'ublast_uniprot'
	  self.printer(self.s2)
	  self.inputo[self.s2] = {}
          savedir = self.workdir + "blast/"
	  savedb = self.workdir + "databases/"
	  os.system("mkdir -p " + savedir)
	  os.system("mkdir -p " + savedb)
  	  
	  if self.config.has_section(self.s2):
	     self.unirun = self.config.get(self.s2,'-run') 		
	     if self.unirun.upper() == "TRUE":
		self.ubevalue = self.config.get(self.s2, 'ubevalue')
		if self.config.has_option(self.s2, '-t'): self.uprocs = self.config.get(self.s2, '-t')
		self.makeudb = self.config.get(self.s2, 'makeudb')
		if self.makeudb.upper() == 'TRUE':
		   self.makedatabase_u = True
		else:
		   self.makedatabase_u = False
		   if not self.config.get(self.s2, 'uniprotudb') or self.config.get(self.s2, 'uniprotudb').upper() == "DEFAULT":		      
		      msg = "If -makeudb option is False, you need to specify a .udb file name"
		      self.invokeexit(self.s2,msg)
		   else:
		      self.uniprotdb_udb = self.config.get(self.s2, 'uniprotudb') 
                
		############## RUN FUNCTION ####################
	        #if self.unirun.upper() == "TRUE":
	     	try:												       
	     	   self.GC_vs_UN = {}										       
	     	   if type(self.genecatalog) == list:								       
	     	      for sample in self.genecatalog:								       
	     		  self.ublast_genes_vs_uniprot(sample, savedir, savedb) 				       
	     	   else:											       
	     	      self.ublast_genes_vs_uniprot(self.genecatalog, savedir, savedb)  # self.GC_vs_UN_ublast	       
	     													       
	     													       
	     	except Exception,e:										       
	     	   print str(e) 										       
	     	   traceback.print_exc()									       


      def ublast_genecat_list(self, samplefasta, savedir, savedb, sample_vs_all):
	     filebase = basename(samplefasta).split(".")[0]  
	     for sample_j in self.genecatalog:
	     	 sbase = ".".join(basename(sample_j).split(".")[0:-1])
	     	 sjudb = savedb  + sbase + ".udb"    	       
	     	 sk_vs_sj_ublast= savedir + filebase +  "_vs_" + sbase + ".ublast"
		 self.ublast(samplefasta, sjudb, sk_vs_sj_ublast)
                 os.system("cat " +  sk_vs_sj_ublast + " >> " + sample_vs_all)
		 os.remove(sk_vs_sj_ublast)
             return 1 #KOSAI# indendation
      


      def ublast_genecatalog(self):
         
          self.s3 = 'ublast_genecatalog'
	  self.inputo[self.s3] = {}
	  savedir = self.workdir + "blast/"
	  savedb = self.workdir + "databases/"
	  os.system("mkdir -p " + savedir)
	  os.system("mkdir -p " + savedb)
	  
	  
	  self.printer(self.s3)
	  if self.config.has_section(self.s3):
	     self.gcrun = self.config.get(self.s3,'-run') 		
	     if self.gcrun.upper() == "TRUE":
		self.ubevalue = self.config.get(self.s3, 'ubevalue')
		if self.config.has_option(self.s3, '-t'): self.uprocs = self.config.get(self.s3, '-t')
		self.makegdb = self.config.get(self.s3, 'makegdb')
		if self.makegdb.upper() == 'TRUE':
		   self.makedatabase_g = True
		else:
		   self.makedatabase_g = False
		   
		   if not self.config.get(self.s3, 'genecatudb') or self.config.get(self.s3, 'genecatudb').upper() == "DEFAULT":
		      checkdb = self.workdir + "databases/" + ".".join(basename(self.genecatalog).split(".")[0:-1]) + ".udb"
		      if os.path.exists(checkdb):
		         self.genecatalog_udb = checkdb
		      else: 
			 msg = "-makegdb option is False, could not find the genecatalog.udb file in databases/ folder, you need to specify a .udb file name"
			 self.invokeexit(self.s2,msg)
		   else:
		      if type(self.genecatalog) == list:
		          self.genecatalog_udb = self.config.get(self.s3, 'genecatudb').strip("\n").split(",")
		      else:    
			  self.genecatalog_udb = self.config.get(self.s3, 'genecatudb')
		   
		
		############## RUN FUNCTION ####################
	        #    if self.gcrun.upper() == "TRUE":
		  
		try:												        
														        
		   #self.printer("Running ublast gene catalog against itself")  				        
		   self.GC_vs_GC = {}										        
		   if self.makedatabase_g:									        
		      if type(self.genecatalog) == list:							        														        
			 for sample_j in self.genecatalog:							        
			     sbase = ".".join(basename(sample_j).split(".")[0:-1])				        
			     udb = savedb  + sbase + ".udb"							        
			     self.ublast_makedb(sample_j, udb)  						        
														        
														        
			 for sample_k in self.genecatalog:							        
			     sample_vs_all = savedir + basename(sample_k).split(".")[0]  + ".allvsall.ublast"        
			     if os.path.exists(sample_vs_all):  						        
				os.remove(sample_vs_all)							        
			     os.system("touch " +  sample_vs_all)						        														        
			     self.ublast_genecat_list(sample_k, savedir, savedb, sample_vs_all) 		        							       
			     self.GC_vs_GC[sample_k] = sample_vs_all						        
		      else:											        
			 genecatalog_base = ".".join(basename(self.genecatalog).split(".")[0:-1])		        
			 self.genecatalog_udb = savedb  + genecatalog_base + ".udb"				        
			 self.ublast_makedb(self.genecatalog, self.genecatalog_udb)				        
			 self.GC_vs_GC_ublast= savedir + self.filebase + ".allvsall.ublast"			        
			 self.ublast(self.genecatalog,self.genecatalog_udb,self.GC_vs_GC_ublast )		        
			 self.GC_vs_GC[self.genecatalog] =  self.GC_vs_GC_ublast				        
														        
														        
		   else:
		      if self.gcrun.upper() == "TRUE":
		         if type(self.genecatalog) == list:
			    for sample_k in self.genecatalog:	
			        filebase = basename(sample_k).split(".")[0]  						        
				sample_vs_all = savedir + basename(sample_k).split(".")[0]  + ".allvsall.ublast"        
				if os.path.exists(sample_vs_all):  						        
				   os.remove(sample_vs_all)							        
				os.system("touch " +  sample_vs_all)
		        	for udb in self.genecatalog_udb:
				    sk_vs_sj_ublast= savedir + filebase +  "_vs_" + sbase + ".ublast"
                                    self.ublast(sample_k, udb, sk_vs_sj_ublast)
				    os.system("cat " +  sk_vs_sj_ublast + " >> " + sample_vs_all)
				    os.remove(sk_vs_sj_ublast)		 		    
				self.ublast_genecat_list(sample_k, savedir, savedb, sample_vs_all) 		 
				self.GC_vs_GC[sample_k] = sample_vs_all						 
		         else: 
		             genecatalog_base = ".".join(basename(self.genecatalog).split(".")[0:-1])			      
			     self.GC_vs_GC_ublast= savedir + self.filebase + ".allvsall.ublast"			      
			     self.ublast(self.genecatalog,self.genecatalog_udb,self.GC_vs_GC_ublast )		      
			     self.GC_vs_GC[self.genecatalog] =  self.GC_vs_GC_ublast				      
									        
		except Exception,e:										        
		   print str(e) 										        
		   traceback.print_exc()									        
          #else: 
	   #  self.invokeexit(self.s3)

     
#####################################################################################################################################
####    															 ####
####    						  BLASTDECIDE FUNCTIONS   						 ####  DONE
####    															 ####
#####################################################################################################################################      
          
      def blastdecide(self, optionslist):                   # abbr. : BD
      
          from bldecide import BlastDecision

	  try:    
              BD_obj = BlastDecision()     	      
              BD_obj.opts=BD_obj.parser.parse_args(optionslist)
	      BD_obj.printer = self.printer
	      self.printer("\n### bldecide.py initialized at " + BD_obj.timestarted  )
	      self.printer("### OPTIONS : " + str(BD_obj.opts) )
	      BD_obj.mainthing()

	       #  obj.parser.print_help()
	  except Exception,e: 
              print str(e)
              traceback.print_exc()
	      sys.exit()


      def blastdecide_uniprot(self,samplefasta, savedir):
          self.printer("IN THE BLASTDECIDE UNIPROT")
          filebase = basename(samplefasta).split(".")[0]
	  lengthfile = self.startdir + filebase + ".length"
	  self.printer("OOOOOOOOOOOOOself.GC_vs_UN[samplefasta]" + self.GC_vs_UN[samplefasta])
	  self.printer("XXXXXXXXXXXXXLengthfile 	       " + lengthfile)
	  self.printer("YYYYYYYYYYYYYSamplefile 	       " + samplefasta)
          if self.config.has_section(self.s4):
	      self.P_obj.updateConfig(self.s4 , "-it", self.GC_vs_UN[samplefasta])
              self.P_obj.updateConfig(self.s4 , "-l", lengthfile)
	  GC_vs_UN_blastdecide =  savedir + filebase + ".genesvsdb.blastdecide"
          self.P_obj.updateConfig(self.s4 , "-o", GC_vs_UN_blastdecide)	     			         
	  self.GC_vs_UN_blastdecide[self.GC_vs_UN[samplefasta]] = [lengthfile, GC_vs_UN_blastdecide]			 
	      							 
          return 1



      def runblastdecide_uni(self):		
 		self.s4 = 'blastdecide_uniprot'
		self.inputo[self.s4] = {}
		self.printer(self.s4)		
		savedir = self.workdir + "blastdecide/"
		os.system("mkdir -p " + savedir)
		self.GC_vs_UN_blastdecide = {}
                imsg = "Input file required. Edit configuration: either give the filename or run ublast ([ublast_uniprot] -run = True)"
                if self.config.has_section(self.s4):
		   if not self.config.has_option(self.s4 ,'-o') or self.config.get(self.s4 ,'-o') == "DEFAULT":
		      outopts = "DEFAULT"
                self.inputblastdecide = []
		############## RUN FUNCTION ####################		

		if self.config.get(self.s4,'-run').upper() == "TRUE":
		   try:		
			if self.config.has_option(self.s4 ,'-it'):
                	   if self.config.get(self.s4 ,'-it') == "DEFAULT":  				   
		        	 if self.unirun.upper() == "TRUE":					       
		        		if type(self.genecatalog) == list:
					   self.outputall_bd =   savedir + self.filebase + ".ALL.genesvsdb.blastdecide.gz"
					   if os.path.exists(self.outputall_bd) : os.remove(self.outputall_bd)
		        		   if outopts == "DEFAULT":
				               for sample in self.genecatalog:
		        			   self.blastdecide_uniprot(sample, savedir) 
						   optslist = self.P_obj.converttolist(self.s4)
						   self.blastdecide(optslist)
						   zipcmd = "gzip -c " + self.GC_vs_UN_blastdecide[self.GC_vs_UN[sample]][1] + " >> " + self.outputall_bd
						   self.inputblastdecide.append(self.GC_vs_UN_blastdecide[self.GC_vs_UN[sample]][1])
						   os.system(zipcmd)
					           #os.remove(self.GC_vs_UN_blastdecide[self.GC_vs_UN[sample]][1])
		        		   else:
					       for sample in self.genecatalog:
					           self.outputall_bd = self.config.get(self.s4 ,'-o')
				        	   self.blastdecide_uniprot(sample, savedir)
				        	   self.inputblastdecide.append(self.GC_vs_UN_blastdecide[self.GC_vs_UN[sample]][1])
						   if "," in self.outputall_bd:
						      omsg = "Input is set to DEFAULT but multiple output names are given. We do not know if \
					        	      each output corresponds to the right input. If input is single, enter a single output \
							      Otherwise, set output to DEFAULT" 
				        	      self.invokeexit(self.s4,omsg)
						   else: 

				        	      self.P_obj.updateConfig(self.s4 , "-o", self.outputall_bd)
						      optslist = self.P_obj.converttolist(self.s4)
						      self.blastdecide(optslist)
						      zipcmd = "gzip -c " + self.GC_vs_UN_blastdecide[self.GC_vs_UN[sample]][1] + " >> " + self.outputall_bd
						      os.system(zipcmd)
					              #os.remove(self.GC_vs_UN_blastdecide[self.GC_vs_UN[sample]][1])
						      

					else: 
		        		       self.blastdecide_uniprot(self.genecatalog, savedir)
					       self.inputblastdecide.append(self.GC_vs_UN_blastdecide[self.GC_vs_UN[self.genecatalog]][1])
					       self.outputall_bd =   savedir + self.filebase + ".ALL.genesvsdb.blastdecide.gz"
					       # self.P_obj.updateConfig(self.s4 , "-o", self.outputall_bd)
					       optslist = self.P_obj.converttolist(self.s4)				       
					       output = self.config.get(self.s4 ,'-o')
					       self.blastdecide(optslist)
					       zipcmd = "gzip -c " + output + " > " + self.outputall_bd
					       os.system(zipcmd)					       
					       #os.remove(output)
					      
					       		        		  
		        	 else: 		               
		        	    self.invokeexit(self.s4,imsg)

			   else:
			         #print self.config.get(self.s4 ,'-it')
				          
		        	 if type(self.genecatalog) == list:
				    self.outputall_bd =   savedir + self.filebase + ".ALL.genesvsdb.blastdecide.gz"
				    if os.path.exists(self.outputall_bd) : os.remove(self.outputall_bd)
				    #print self.outputall_bd
				    for line in file(self.config.get(self.s4 ,'-it')):
				          inpt = line.split("\t")[0]
					  lengthfile = line.split("\t")[1].strip("\n")
				          self.P_obj.updateConfig(self.s4 , "-it", inpt)
                                          self.P_obj.updateConfig(self.s4 , "-l", lengthfile)
				          filebase = basename(inpt).split(".")[0]
					  #print filebase
                                          GC_vs_UN_blastdecide =  savedir + filebase + ".genesvsdb.blastdecide"
					  #print GC_vs_UN_blastdecide
                                          self.P_obj.updateConfig(self.s4 , "-o", GC_vs_UN_blastdecide)	
					  optslist = self.P_obj.converttolist(self.s4)
                                          self.blastdecide(optslist)
					  
					  zipcmd = "gzip -c " + GC_vs_UN_blastdecide + " >> " + self.outputall_bd
					  os.system(zipcmd)
					  self.outputall_bd	
					  #os.remove(GC_vs_UN_blastdecide)
			               			  
				 
				 else:
				 
				     self.GC_vs_UN_ublast = self.config.get(self.s4 ,'-it')	
				     			 
				     output = savedir + self.filebase + ".ALL.genesvsdb.blastdecide"
				     self.outputall_bd = output + ".gz"
				     self.P_obj.updateConfig(self.s4 , "-o", output)		 
				     lengthfile = self.startdir + self.filebase + ".length"
				     if os.path.exists(lengthfile):
                        		self.P_obj.updateConfig(self.s4 , "-l", lengthfile)
					optslist = self.P_obj.converttolist(self.s4)
					self.blastdecide(optslist)
					zipcmd = "gzip -c " + output + " > " + self.outputall_bd
					os.system(zipcmd)
					os.remove(output)

				     else: 
					lmsg = "Expected length file : " + lengthfile + " NOT FOUND.Run fasta2index on genecatalog and put on workdir/samples folder. OR start over with default"
					self.invokeexit(self.s4,lmsg)


			else:		      
			   self.invokeexit(self.s4,imsg)


		     
		   except Exception,e:
		      print str(e)
		      traceback.print_exc()
		      

      def blastdecide_genecatalog(self, samplefasta, savedir):
          self.printer("IN THE BLASTDECIDE GC")
          filebase = basename(samplefasta).split(".")[0]
	 
	  lengthfile = self.startdir + filebase + ".length"
          if self.config.has_section(self.s5):
	      self.P_obj.updateConfig(self.s5 , "-it", self.GC_vs_GC[samplefasta])
              self.P_obj.updateConfig(self.s5 , "-l", lengthfile)
	  obase = basename(self.GC_vs_GC[samplefasta]).split(".")[0]	  
	  GC_vs_GC_blastdecide =  savedir + obase + ".allvsall.blastdecide"
          self.P_obj.updateConfig(self.s5 , "-o", GC_vs_GC_blastdecide)	     			         
	  self.GC_vs_GC_blastdecide[self.GC_vs_GC[samplefasta]] = [lengthfile, GC_vs_GC_blastdecide]			 
	      							 
          return 1


      def runblastdecide_gc(self):
	  self.s5 = 'blastdecide_genecatalog'
	  self.inputo[self.s5] = {}
	  self.printer(self.s5)
	  savedir = self.workdir + "blastdecide/"
	  os.system("mkdir -p " + savedir)
          self.GC_vs_GC_blastdecide={}
          imsg = "Input file required. It is not set to Default. Edit configuration: either give the filename or run ublast ([ublast_genecatalog] -run = True)"
          if self.config.has_section(self.s5):
	     if not self.config.has_option(self.s5 ,'-o') or self.config.get(self.s5 ,'-o') == "DEFAULT":
	  	outopts = "DEFAULT"
    	     if self.config.get(self.s5,'-run').upper() == "TRUE":
		try:

	            if not self.config.has_option(self.s5 ,'-it') : self.invokeexit(self.s5,imsg)
	            if self.config.get(self.s5 ,'-it') == "DEFAULT":

		            if self.gcrun.upper() == "TRUE":
		         	   if type(self.genecatalog) == list: 

		         	       self.outputall_gc = savedir + self.filebase + ".ALL.genesvsgenes.blastdecide.gz"
				       if os.path.exists(self.outputall_gc) : os.remove(self.outputall_gc)
		         	       if outopts == "DEFAULT":
		         		   for sample in self.genecatalog:
		         		       self.blastdecide_genecatalog(sample, savedir) 
		         		       optslist = self.P_obj.converttolist(self.s5)
		         		       self.blastdecide(optslist)
					       zipcmd = "gzip -c " + self.GC_vs_GC_blastdecide[self.GC_vs_GC[sample]][1] + " >> " + self.outputall_gc
					       os.system(zipcmd)
					       os.remove(self.GC_vs_GC_blastdecide[self.GC_vs_GC[sample]][1])
		         	       else:
				           for sample in self.genecatalog:
					       self.outputall_gc = self.config.get(self.s5 ,'-o')
		         		       self.blastdecide_genecatalog(sample, savedir)

		         		       if "," in self.outputall_gc:
		         			  omsg = "Input is set to DEFAULT but multiple output names are given. We do not know if \
		         				  each output corresponds to the right input. If input is single, enter a single output \
		         				  Otherwise, set output to DEFAULT" 
		         			  self.invokeexit(self.s5,omsg)
		         		       else: 
		         			   self.P_obj.updateConfig(self.s5 , "-o", self.outputall_gc)
		         			   optslist = self.P_obj.converttolist(self.s5)
		         			   self.blastdecide(optslist)  			
		         	   else: 
		         		  self.blastdecide_genecatalog(self.genecatalog, savedir)
		         		  optslist = self.P_obj.converttolist(self.s5)
					  self.outputall_gc = self.config.get(self.s5 ,'-o')
		         		  self.blastdecide(optslist)				       
		            else:  			    
		               self.invokeexit(self.s5,imsg)


	            else:

			   if type(self.genecatalog) == list:
		     	       self.outputall_gc =   savedir + self.filebase + ".ALL.genesvsgenes.blastdecide.gz"
			       if os.path.exists(self.outputall_gc) : os.remove(self.outputall_gc)
		     	       for line in file(self.config.get(self.s5 ,'-it')):			       
		     		     inpt = line.split("\t")[0]
		     		     lengthfile = line.split("\t")[1].strip("\n")
		     		     self.P_obj.updateConfig(self.s5 , "-it", inpt)
                     		     self.P_obj.updateConfig(self.s5 , "-l", lengthfile)
		     		     filebase = basename(inpt).split(".")[0]
		     		     #print filebase
                     		     GC_vs_GC_blastdecide =  savedir + filebase + ".genesvsgenes.blastdecide"
                     		     self.P_obj.updateConfig(self.s5 , "-o", GC_vs_GC_blastdecide) 
		     		     optslist = self.P_obj.converttolist(self.s5)
                     		     self.blastdecide(optslist)		     	       
		     		     zipcmd = "gzip -c " + GC_vs_GC_blastdecide + " >> " + self.outputall_gc
		     		     os.system(zipcmd)  
		     		     #os.remove(GC_vs_GC_blastdecide)


			   else:

				     self.GC_vs_GC_ublast = self.config.get(self.s5 ,'-it')
				     self.outputall_gc = savedir + self.filebase + ".ALL.genesvsgenes.blastdecide"
				     self.P_obj.updateConfig(self.s5 , "-o", self.outputall_gc)
				     lengthfile = self.startdir + self.filebase + ".length"
				     if os.path.exists(lengthfile):
                			self.P_obj.updateConfig(self.s5 , "-l", lengthfile)
					optslist = self.P_obj.converttolist(self.s5)
					self.blastdecide(optslist)
				     else: 
					lmsg = "Expected length file : " + lengthfile + " NOT FOGCD.RGC fasta2index on genecatalog and put on workdir/samples folder. OR start over with default"
					self.invokeexit(self.s5,lmsg)
		except Exception,e:
		     print str(e)
		     traceback.print_exc()  
      
      
      
      
      
      
      
      		   
#####################################################################################################################################
####    															 ####
####    						  CLUSTERING FUNCTIONS   						 ####  DONE
####    															 ####
#####################################################################################################################################      
 
      def clustering(self,optionslist):
	  from clustergenes import ClusterGenes
	
	  try:
              CL_obj = ClusterGenes()
	      CL_obj.opts = CL_obj.parser.parse_args(optionslist)
	      CL_obj.printer = self.printer
	      self.printer("\n### clustergenes.py initialized at " + CL_obj.timestarted)  
	      self.printer("### OPTIONS : " + str(CL_obj.opts))
              CL_obj.mainthing()
	  #except IOError as i:
    	   #   print "I/O error({0}): {1}".format(i.errno, i.strerror)
	  except Exception,e: 
	      print str(e)
	      traceback.print_exc()
	      sys.exit()	  


      def clustergenes(self):              ############ What happens if you use smaller gene names for the graph
	  self.s6 = 'clustergenes'
	  self.printer(self.s6)
	  self.inputo[self.s6] = {}
	  savedir = self.workdir + self.s6 + "/"
	  os.system("mkdir -p " + savedir)
	  if self.config.has_section(self.s6):

	         if self.config.get(self.s6 ,'-p') == "DEFAULT": #	    INPUT OPTIONS
		 	if not self.config.has_section(self.s5): invokeexit(self.s5)
		        if self.config.get(self.s5,'-run').upper() == "TRUE" or self.config.get(self.s3,'-run').upper() == "TRUE":
			   self.P_obj.updateConfig(self.s6 , "-p", self.outputall_gc)			
		 else:
                 	self.GC_vs_GC_blastdecide = self.config.get(self.s6 ,'-p')



                 if self.config.has_option(self.s6 ,'-o'):       #            OUTPUT OPTIONS
                    outname = self.config.get(self.s6 ,'-o')

  		 if outname.upper() == "DEFAULT" or not outname:
		      if self.config.get(self.s6 ,'-a') == 'True':   
                	 self.clusteroutput = savedir + self.filebase 
		      else: 
			 b = self.config.get(self.s6 ,'-b')
			 self.clusteroutput =  savedir + self.filebase + ".b" + b              #MGfunc_results/clustergenes/genecat.b3.all.sh
		      self.P_obj.updateConfig(self.s6 , "-o", self.clusteroutput)		       
		 else:
		      self.clusteroutput = savedir + self.config.get(self.s6 ,'-o')


                 self.CLopts = self.P_obj.converttolist(self.s6)   
	         
	
                 ############## RUN FUNCTION ####################
		 if self.config.get(self.s6,'-run').upper() == "TRUE":	
		    try:
			 self.clustering(self.CLopts)
                    except Exception,e:
			 print str(e)
			 traceback.print_exc()

#####################################################################################################################################
####    															 ####
####    						 CLUSTERFILTER FUNCTIONS   						 ####  DONE
####    															 ####
#####################################################################################################################################      
 

      def filter(self, optslist):
      
          import clusterFilter
	  try:	
	      CF_obj = clusterFilter.main()
              CF_obj.args = CF_obj.parser.parse_args(optslist)
              CF_obj.printer = self.printer
	      self.printer("\n### clusterFilter.py initialized at " + CF_obj.timestarted  )
              self.printer("### OPTIONS : " + str(CF_obj.args))
              CF_obj.mainthing()
 
	  except IOError as i:
    	      print "I/O error({0}): {1}".format(i.errno, i.strerror)
	  except Exception,e: 
	      print str(e)
	      traceback.print_exc()
	      sys.exit()	  


      def clusterFilter(self):
      
	  self.s7 = 'clusterFilter'
	  self.inputo[self.s7] = {}
	  savedir = self.workdir + self.s7 + "/"
	  os.system("mkdir -p " + savedir)
	  self.printer(self.s7)
	  if self.config.has_section(self.s7):


		############## INPUT FILE ####################################################
	        
		if self.config.get(self.s7 ,'-i') == "DEFAULT":		    
	              			 
                           if self.config.get(self.s6 ,'-m') == 'reciprocal':	 
                        	 if self.config.get(self.s6 ,'-a') == 'True':   
                        	      self.clusterfile = self.clusteroutput + ".all.scc"
                        	 else:  					       
                        	      self.clusterfile = self.clusteroutput + ".bh.scc" 
				      
				      #print self.clusterfile 
                           else:							 
                        	 self.clusterfile = self.clusteroutput + ".trihits.scc"  
	              
		           self.P_obj.updateConfig(self.s7, "-i", self.clusterfile)
					 
	        else:								 
	             self.clusterfile = self.config.get(self.s7 ,'-i')
		           
	     
	     
		############## FILTERING OPTIONS ###############################################

		#if self.config.get(self.s7 ,'-i') == "DEFAULT":


		############## OUTPUT OPTIONS ##################################################

		if self.config.has_option(self.s7 ,'-o'):
		   self.s7out= self.config.get(self.s7 ,'-o')
		
		
	        if not self.s7out or self.s7out == "DEFAULT":
			self.filteredclusters = savedir + self.filebase
			self.P_obj.updateConfig(self.s7 , "-o", self.filteredclusters)
		else:
			self.filteredclusters = self.config.get(self.s7 ,'-o')

		
        	############## RUN FUNCTION ####################################################
		if self.config.get(self.s7,'-run').upper() == "TRUE":
		   try: 
		        #self.config.remove_option(self.s7,'run')
		        CLFopts = self.P_obj.converttolist(self.s7)
			#print CLFopts
	        	self.filter(CLFopts)
			self.filteredclusters = self.filteredclusters + ".filtered"
        	   except Exception,e:
	        	print str(e)
	        	traceback.print_exc()
			
		else:   ##CHECH THIS PART
		    self.filteredclusters = self.clusterfile  
		    self.printer("### ClusterFilter bypassed!")
	      	     

#####################################################################################################################################
####    															 ####
####    	                                 PARSE CLUSTER INFO FUNCTIONS   						 ####  DONE (Upgrade: ADD THREADING-done)
####    															 ####  add -g, update script
#####################################################################################################################################      
 
 
      def getuniprotinfo(self, optslist):
          from parseclinfo import GetUniprotInfo 
	  try:
              M_obj = GetUniprotInfo()
              M_obj.opts = M_obj.parser.parse_args(optslist)
	      M_obj.printer = self.printer
	      self.printer("\n### parseclinfo.py initialized at " + M_obj.timestarted )
	      self.printer("### OPTIONS : " + str(M_obj.opts))
              M_obj.mainthing()
	  #except IOError as i:
	   #   print "I/O error({0}): {1}".format(i.errno, i.strerror)
	  except Exception,e: 
	      print str(e)
	      traceback.print_exc()
	      sys.exit()

      def parseclusterinfo(self):	      
	  self.s8 = 'parseclusterinfo'
	  self.inputo[self.s8] = {}
	  self.printer(self.s8)
	  savedir = self.workdir + self.s8 + "/"
	  os.system("mkdir -p " + savedir)
	  if self.config.has_section(self.s8):
	         if self.config.get(self.s8 ,'-c') == "DEFAULT":
		    self.inputcluster = self.filteredclusters 
		    self.P_obj.updateConfig(self.s8 , "-c", self.inputcluster)
		 else:
		    self.inputcluster = self.config.get(self.s8 ,'-c')
	         
	   
	         if self.config.get(self.s8 ,'-b') == "DEFAULT":
		    if self.config.get(self.s4,'-run').upper() == "TRUE" or self.config.get(self.s2,'-run').upper() == "TRUE": 
		       #self.inputblastdecide = self.outputall_bd
		       self.P_obj.updateConfig(self.s8 , "-b", self.inputblastdecide)
		 else:
		    if self.config.has_option(self.s8 ,'-bf'): 
		       bhfolder = self.config.get(self.s8 ,'-bf')
		       self.config.remove_option(self.s8,'-bf')
		       if "," in  self.config.get(self.s8 ,'-b'):
		          self.inputblastdecide = []
		          for fil in self.config.get(self.s8 ,'-b').split(","):
			      realfil = bhfolder + "/" + fil
			      self.inputblastdecide.append(realfil)
		       self.P_obj.updateConfig(self.s8 , "-b", " ".join(self.inputblastdecide))
		       self.printer(self.config.get(self.s8 ,'-b'))
		    else:
		       if "," in self.config.get(self.s8 ,'-b'):
			 self.inputblastdecide =  self.config.get(self.s8 ,'-b').split(",")
			 self.P_obj.updateConfig(self.s8 , "-b", " ".join(self.inputblastdecide))
		       else:
			 self.inputblastdecide = [self.config.get(self.s8 ,'-b')]
		 
		 
		 ######################### ALL GENE IDS ######################################
		 if self.config.has_option(self.s8 ,'-g'):
		    if self.config.get(self.s8 ,'-g') == "DEFAULT":  
		       self.allgeneidsfile = self.startdir + "allgeneids" 
		       self.P_obj.updateConfig(self.s8 , "-g", self.allgeneidsfile)
		 
		    
		 ######################### UNIPROT INPUT - TAB and INDEXES ###################   
		 self.P_obj.updateConfig(self.s8 , "-ut", self.uniprottab)	         
		 
		 if self.config.has_option(self.s8 ,'-uif'):
		    if self.config.get(self.s8 ,'-uif') == "DEFAULT":
		       testloc = "/".join(os.path.abspath(self.uniprotdat).split("/")[:-1]) + "/"
		       if os.path.exists(testloc + "uniprotindexes"):
		          if os.path.exists(testloc + "uniprotindexes" + "A.tab.index"):
			     self.P_obj.updateConfig(self.s8 , "-uif", testloc + "uniprotindexes")
			  else: 
			     self.printer("Index folder found but index files names are wrong\n")
			     self.printer("Using large index file, no index folder found\n")
                             self.P_obj.updateConfig(self.s8 , "-ui", self.uniprottabindex)
		       else: 
		          self.printer("Using large index file, no index folder found\n")
			  self.P_obj.updateConfig(self.s8 , "-ui", self.uniprottabindex)
		       
		 else: 
		    self.printer("Using large index file, no index folder found\n")
		    self.P_obj.updateConfig(self.s8 , "-ui", self.uniprottabindex)	   
	         
		 enr = self.config.get(self.s8 ,'-e')
		 if enr.upper() == "DEFAULT" or enr.upper() == "TRUE":
		    if self.preclusterfile: 
		       self.precluster = self.preclusterfile
		    else:   
		       self.precluster = self.genecatpath + self.filebase + ".precluster"
		       if not os.path.exists(self.precluster):
		          self.invokeexit(self.s8,"Precluster file not found! Check your options!\n")
		    self.P_obj.updateConfig(self.s8 , "-e", self.precluster)
		 if self.config.get(self.s8 ,'-e').upper() == "FALSE":
		    self.config.remove_option(self.s8,'-e')
		 else: 
		    self.precluster = self.config.get(self.s8 ,'-e')
		    		 

		 if not self.config.has_option(self.s8 ,'-o') or self.config.get(self.s8 ,'-o') == "DEFAULT":
			self.Poutbase = savedir + self.filebase
			self.P_obj.updateConfig(self.s8 , "-o", self.Poutbase)
		 else:
			self.Poutbase = savedir + self.config.get(self.s8 ,'-o')
			self.P_obj.updateConfig(self.s8 , "-o", self.Poutbase)

		 self.PIopts = self.P_obj.converttolist(self.s8)
		 #print self.PIopts
		 
		 ############## RUN FUNCTION ####################       
		
		 if self.config.get(self.s8,'-run').upper() == "TRUE":
		    try:
	        	 self.getuniprotinfo(self.PIopts)

        	    except Exception,e:
	        	 print str(e)
	        	 traceback.print_exc()



#####################################################################################################################################
####    															 ####
####    						GENE ENRICHMENT		   						 ####  DONE
####    															 ####
#####################################################################################################################################    



      def gene_enrichment(self, optslist):
          
	  try:
	      import goawrapper
              G_obj = goawrapper.rungoatools()
              G_obj.args = G_obj.parser.parse_args(optslist)
	      G_obj.printer = self.printer
	      self.printer("\n### goawrapper.py initialized at " + G_obj.timestarted )
	      self.printer("### OPTIONS : " + str(G_obj.args))
              G_obj.mainthing()
	  except IOError as i:
	      print "I/O error({0}): {1}".format(i.errno, i.strerror)
	  except Exception,e: 
	      print str(e)
	      traceback.print_exc()
	      sys.exit()



      def enrichment(self):
          self.s9 = 'enrichment'
	  self.inputo[self.s9] = {}
	  self.printer(self.s9)
	  savedir = self.workdir + "enrichment/"
	  os.system("mkdir -p " + savedir)

	  if self.config.has_section(self.s9):
	         if self.config.get(self.s9 ,'-i') == "DEFAULT":
		    self.inputcluster = self.Poutbase + ".allwithhits.txt"
		    self.P_obj.updateConfig(self.s9 , "-i", self.inputcluster)
		 else:
		    self.inputcluster = self.config.get(self.s9 ,'-i')
	   
	   
	         if self.config.get(self.s9 ,'-a') == "DEFAULT" and self.config.get(self.s9 ,'-p') == "DEFAULT":
		    self.aso = self.Poutbase + ".association"
		    self.popl = self.Poutbase + ".population"
		    self.P_obj.updateConfig(self.s9 , "-a", self.aso)
		    self.P_obj.updateConfig(self.s9 , "-p", self.popl)
		 else:
		    self.aso = self.config.get(self.s9 ,'-a')
		    self.popl = self.config.get(self.s9 ,'-p')
	  

                 if self.config.get(self.s9 ,'-o') == "DEFAULT" or not self.config.get(self.s9 ,'-o'):
                    self.enrichout = savedir + self.filebase
		    self.P_obj.updateConfig(self.s9 , "-o", self.enrichout)
		 else: 
                    self.enrichout = savedir + self.config.get(self.s9 ,'-o')
		    self.P_obj.updateConfig(self.s9 , "-o", self.enrichout)
                 
		 
		 
		 self.ENopts = self.P_obj.converttolist(self.s9)
		 ############## RUN FUNCTION ####################       
		
		 if self.config.get(self.s9,'-run').upper() == "TRUE":
		    try:
	        	 self.gene_enrichment(self.ENopts)

        	    except Exception,e:
	        	 print str(e)
	        	 traceback.print_exc()



#####################################################################################################################################
####    															 ####
####    				EXTRACT CLUSTER SEQUENCES FROM FASTA FILE   						 ####  DONE
####    															 ####
#####################################################################################################################################      

      def extract(self, inputlist,THREAD_LIMIT):
	  self.jobs = Queue.Queue(0)
	  self.singlelock = threading.Lock()
	  myrange =  range(THREAD_LIMIT)
          self.ThreadNo_e = list(myrange[::-1])
          self.target = self.Cluster2Fasta	 
	  #self.StartQueueAndThreads(inputlist, THREAD_LIMIT)
	  try: 
	     self.StartQueueAndThreads(inputlist, THREAD_LIMIT)
	  except:
	     traceback.print_exc()
	     sys.exit()
          
      def Cluster2Fasta(self):                   ## Target function	  
	  import cluster2fasta
	  import traceback
	  therandom = self.ThreadNo_e.pop() + 1
	  self.singlelock.acquire()
          self.printer("starting thread number:" + str(therandom))
          self.singlelock.release()
	  #infile = self.jobs.get(True,1)
	  while 1: 
	   
	       try:
	       		   
        	   infile = self.jobs.get(True,1)
		   self.singlelock.acquire()
		   self.P_obj.updateConfig("cluster2fasta", "-c", infile)		   
 		   self.P_obj.updateConfig("cluster2fasta", "-o", self.extoutlist[infile])
		   self.config.remove_option("cluster2fasta","-t")
		   Eopts = self.P_obj.converttolist("cluster2fasta")
		   self.singlelock.release()
		   		   
        	   E_obj = cluster2fasta.main()
        	   E_obj.args = E_obj.parser.parse_args(Eopts)
		   E_obj.printer = self.printer
        	   self.printer("\n### cluster2fasta.py initialized at "+ E_obj.timestarted + "\n")
        	   self.printer("### OPTIONS: "+str(E_obj.args)+"\n")
        	   E_obj.mainthing()	
		   outfile = self.extoutlist[infile]
		   os.system("cat " + outfile + ".genecatalog.fasta " + outfile + ".uniprotids.fasta > " + outfile + ".fasta" )		   
		   self.jobs.task_done()

	       except: 
	           traceback.print_exc()
	           break      
		      


      def runcluster2fasta(self):
	  self.s10 = 'cluster2fasta'
	  self.inputo[self.s10] = {}
          self.printer(self.s10)
	  savedir = self.workdir + self.s10 + "/"
	  os.system("mkdir -p " + savedir)
	  self.extoutlist={}
          if self.config.has_section(self.s10):
                
		 
		if not self.config.get(self.s10 ,'-c') or self.config.get("cluster2fasta" ,'-c') == 'DEFAULT' :
		      if int(self.config.get(self.s10 ,'-t')) > 1:
			 self.getclusterids_fromCLfile(self.Poutbase + ".allinfo.txt")
			 CL1 = self.Poutbase + ".nohits.txt"           #N
			 CL2 = self.Poutbase + ".semihits.txt"         #S   
			 CL3 = self.Poutbase + ".allhits.txt"	       #A   
			 inputlist = [CL1, CL2, CL3]
			 
		      else: 
		         inputlist = [self.Poutbase + ".allinfo.txt"]
		      
		else: 
		      if "," in self.config.get("cluster2fasta" ,'-c'):
	                 inz = self.config.get("cluster2fasta" ,'-c')
			 inputlist = inz.rstrip().split(",")
		      else: 
		         inputlist = [self.config.get("cluster2fasta" ,'-c')]			      
                      
		THREAD_LIMIT = int(self.config.get(self.s10 ,'-t'))
		self.c2f_thread = THREAD_LIMIT
                outlist = []
		
		if not self.config.get(self.s10 ,'-o') or self.config.get("cluster2fasta" ,'-o') == 'DEFAULT' :
		      for i in inputlist: 
		          ibase= basename(i)
		          outlist.append(savedir + ibase[:-4])
		      
		else: 
                     if "," in self.config.get("cluster2fasta" ,'-o'):
	                 out = self.config.get("cluster2fasta" ,'-o')
			 outlist = out.rstrip().split(",")
		     else: 
		         outlist = [self.config.get("cluster2fasta" ,'-o')]	
		     
		
		for x in range(0, len(inputlist)): 
		    self.extoutlist[inputlist[x]] = outlist[x]
		    
		if type(self.genecatalog) == list:
		   if self.config.get(self.s10 ,'-ki').upper() != "DEFAULT":
		      if self.config.has_option(self.s10 ,'-sfi'): self.config.remove_option(self.s10,"-sfi")
		   else: 
		      if not self.config.has_option(self.s10 ,'-sfi') or self.config.get(self.s10 ,'-sfi').upper() == "DEFAULT":
			 self.sfi = savedir + self.filebase + ".fastaindex.list"
			 with open(self.sfi,"w") as sfi_fh:
			      for fasta in self.genecatalog:
		        	  fastabase = basename(fasta).split(".")[0]
				  indexfile= self.startdir + fastabase + ".indexed"
                        	  sfi_fh.write(indexfile + "\t" + fasta + "\n")
			 self.P_obj.updateConfig(self.s10 , "-sfi", self.sfi)
			 self.config.remove_option(self.s10,"-ki")
			 self.config.remove_option(self.s10,"-kf")
		      else: 
		         self.sfi = self.config.get(self.s10 ,'-sfi')
		      
		else:
		   if self.config.has_option(self.s10 ,'-sfi'): self.config.remove_option("cluster2fasta","-sfi")
		   ### default values for -ki and -kf			    
		   if self.config.has_option(self.s10 ,'-ki') and self.config.get("cluster2fasta" ,'-ki').upper() == 'DEFAULT' :
		      self.P_obj.updateConfig(self.s10 , "-ki", self.startdir +  self.filebase + ".indexed")
		      #print self.startdir +  self.filebase + ".indexed"

		   if self.config.has_option(self.s10 ,'-kf') and self.config.get("cluster2fasta" ,'-kf').upper() == 'DEFAULT' :
		      self.P_obj.updateConfig(self.s10 , "-kf", self.genecatalog)
		   
		   ### default values for -ui and -uf
                   if self.config.has_option(self.s10 ,'-ui') and self.config.get("cluster2fasta" ,'-ui').upper() == 'DEFAULT' :
		      self.P_obj.updateConfig(self.s10 , "-ui", self.uniprotfasta + ".indexed")

		   if self.config.has_option(self.s10 ,'-uf') and self.config.get("cluster2fasta" ,'-uf').upper() == 'DEFAULT' :
		      self.P_obj.updateConfig(self.s10 , "-uf", self.uniprotfasta)
		
		
		######################### RUN FUNCITON ############################
		if self.config.get(self.s10,'-run').upper() == "TRUE":
		   try:
			self.extract(inputlist,THREAD_LIMIT)

        	   except Exception,e:
	        	print str(e)
	        	traceback.print_exc()
			sys.exit()

#####################################################################################################################################
####    															 ####
####    				SPLIT CLUSTERS FASTA FILE   								 ####  DONE
####    															 ####
#####################################################################################################################################      


      def splitc(self, inputdict,THREAD_LIMIT):
          inputlist = inputdict.keys()	
	  self.jobs = Queue.Queue(0)
	  self.singlelock = threading.Lock()
	  myrange =  range(THREAD_LIMIT)
          self.ThreadNo_e = list(myrange[::-1])
          self.target = self.SplitClusters	
	  
	  try: 
	     self.StartQueueAndThreads(inputlist, THREAD_LIMIT)
	  except:
	     traceback.print_exc()
	     sys.exit()
	  
      def getclusterids_fromfasta(self, filename):
          os.system('grep ">" ' + filename + ' | cut -f1 -d":" | tr -d ">" > ' + filename[0:-4] + '.fasta.clid'  )

      def getclusterids_fromCLfile(self, filename):
          filebase= ".".join(filename.split(".")[:-2])
          os.system("grep -P '\tN\t' " + filename + " > " + filebase + ".nohits.txt")
	  os.system("grep -P '\tS\t' " + filename + " > " + filebase + ".semihits.txt")
	  os.system("grep -P '\tA\t' " + filename + " > " + filebase + ".allhits.txt")
	  os.system("cut -f1 " + filebase + ".nohits.txt > " + filebase + ".nohits.clid")
	  os.system("cut -f1 " + filebase + ".semihits.txt > " + filebase + ".semihits.clid")
	  os.system("cut -f1 " + filebase + ".allhits.txt > " + filebase + ".allhits.clid")
	  
	  
	  
      def SplitClusters(self):
           
           import splitclusters
	   therandom = self.ThreadNo_e.pop() + 1
	   self.singlelock.acquire()
           self.printer("starting thread number:" + str(therandom))
           self.singlelock.release()
	   while 1: 

		  try:    
		         infile = self.jobs.get(True,1)
			 self.singlelock.acquire()
		         self.P_obj.updateConfig("splitclusters", "-i", infile)			   	      
                         self.P_obj.updateConfig("splitclusters", "-cl", self.splitdict[infile])     
		         self.P_obj.updateConfig("splitclusters", "-o", self.splitpaths[infile])
			 self.config.remove_option("splitclusters","-t")
		         SPopts = self.P_obj.converttolist("splitclusters")
			 self.singlelock.release()
        	         obj = splitclusters.Split()
			 obj.printer = self.printer        
        	         obj.opts=obj.parser.parse_args(SPopts)   # mode = splitclusters or changeheaders
        	         #self.printer("\n### splitclusters initialized at " + obj.timestarted )
        	         #self.printer("### OPTIONS : " + str(obj.opts))
        	         obj.mainthing()
                         self.jobs.task_done()
			 
	          except: 
	              break 			 
#		  except Exception,e: 
#        	      print str(e)
#        	      traceback.print_exc()
#		      sys.exit()




      def splittclusters(self):
	  self.s11 = 'splitclusters'
	  self.inputo[self.s11] = {}
          self.printer(self.s11)
	  #print "I'm running splitclusters"
	  savedir = self.workdir + self.s11 + "/"
	  os.system("mkdir -p " + savedir)
	  CL1 = self.Poutbase + ".nohits.txt"  
	  CL2 = self.Poutbase + ".semihits.txt"
	  CL3 = self.Poutbase + ".allhits.txt"	 
	  #print CL1
	  #print CL2
	  #print CL3
	  
	  self.splitdict = {} 
	  self.splitpaths = {}
          if self.config.has_section(self.s11):
		## INPUT 1 FASTA FILES
		if self.config.get(self.s11 ,'-i') == 'DEFAULT' :
		   if self.c2f_thread > 1:      		          
			  F1 = self.extoutlist[CL1] + ".fasta"
			  F2 = self.extoutlist[CL2] + ".fasta"
			  F3 = self.extoutlist[CL3] + ".fasta"		      
			  inputlist = [F1, F2, F3]			  
			  
	           else:
		          inputlist = [self.extoutlist[self.Poutbase + ".allinfo.txt"] + ".fasta"]
		      
	        else:
		      if "," in self.config.get(self.s11 ,'-i'):
	                 ins = self.config.get(self.s11 ,'-i')
			 inputlist = ins.rstrip().split(",")
		      else: 
		         inputlist = [self.config.get(self.s11 ,'-i')]			      
 		
		
		## INPUT 2 CLID FILES     
	        if not self.config.get(self.s11 ,'-cl') == 'DEFAULT' :
		      if "," in self.config.get(self.s11 ,'-cl'):
	                 inc = self.config.get(self.s11 ,'-cl')
			 clfiles = inc.rstrip().split(",")
		      else: 
		         clfiles = [self.config.get(self.s11 ,'-cl')]
		else: 
		      
		      if self.c2f_thread > 1: 
			 FCL1 = CL1[0:-4] + ".clid"
			 FCL2 = CL2[0:-4] + ".clid"
			 FCL3 = CL3[0:-4] + ".clid"
			 clfiles = [FCL1, FCL2,FCL3]
			 
		      else: 
		         if self.config.get(self.s10,'-run').upper() == "TRUE": os.system("cut -f1 " + self.Poutbase + ".allinfo.txt > " + self.Poutbase + ".allinfo.clid")
		         clfiles = [self.Poutbase + ".allinfo.clid"]
		
		opti =  self.config.get(self.s11 ,'-i')
		optcl = self.config.get(self.s11 ,'-cl')
		
		for x in range(0, len(inputlist)): 
		    self.splitdict[inputlist[x]] = clfiles[x]
		    infile = inputlist[x]
		    m1 = re.search("nohits", infile)
		    m2 = re.search("semihits", infile)
		    m3 = re.search("allhits", infile)	
		    if self.c2f_thread > 1:
		       if m1: 
		    	   path = savedir + "nohits/"			   
		       if m2: 
		    	   path = savedir+ "semihits/"     
		       if m3: 
		    	   path = savedir + "allhits/"        
		    else: 
		    	if "allinfo" in infile:
		    	   path = mainpath	    
                    os.system("mkdir -p " + path )
		    self.splitpaths[infile] = path
		
		thlimit = int(self.config.get(self.s11 ,'-t'))
				
		######################### RUN FUNCITON ############################
		if self.config.get(self.s11,'-run').upper() == "TRUE":
		   try:		   
	        	self.splitc(self.splitdict,thlimit)		      
        	   except Exception,e:
	        	print str(e)
	        	traceback.print_exc()
	   
#####################################################################################################################################
####    															 ####
####    				CHANGE FASTA HEADERS	   								 ####  DONE
####    															 ####
#####################################################################################################################################      


      def changeh(self,inputlist):
         THREAD_LIMIT = self.THREADLIMIT_ch
         self.jobs = Queue.Queue(0)
	 self.singlelock = threading.Lock()
         self.target = self.changeheaders_inside	 
	 self.StartQueueAndThreads(inputlist, THREAD_LIMIT)
      
      def changeheaders_outside(self):
          while 1: 
	       try:    
		    
		   infile = self.jobs.get(True,1)		   
		   self.singlelock.acquire()
		   self.P_obj.updateConfig('changeheaders', "-i", infile)		   		
		   opts = self.P_obj.converttolist('changeheaders')
		   self.singlelock.release()
		   
		   optslist = " ".join(opts)
		   
	           cmdline = "python2.7 " + self.pkgdir + "changeheaders.py " + optslist 
	           os.system(cmdline)
		   
		   self.jobs.task_done()
			 
	       except Exception,e: 
        	    print str(e)
		    if str(e) == "Empty":
        	    #traceback.print_exc()
		       sys.exit()
		   
		   
      def changeheaders_inside(self):
          from changeheaders import Change 
	  while 1: 
	       try:    
		   infile = self.jobs.get(True,1)
		   base = ".".join(basename(infile).split(".")[:-1])		   
		   self.singlelock.acquire()
		   N= self.chfolder + base + ".new.fsa"
		   HT= self.chfolder + base + ".headertable"  
		   self.P_obj.updateConfig('changeheaders', "-i", infile)
		   self.P_obj.updateConfig('changeheaders', "-n", N)
		   self.P_obj.updateConfig('changeheaders', "-ht", HT)		   		
		   opts = self.P_obj.converttolist('changeheaders')
		   self.singlelock.release()	      
        	   obj = Change()        
        	   obj.opts=obj.parser.parse_args(opts)
		   self.singlelock.acquire()
		   obj.printer = self.printer2
        	   obj.printer("\n### " + sys.argv[0] + " initialized at " + obj.timestarted  )
        	   obj.printer("### OPTIONS : " + str(obj.opts))
        	   self.singlelock.release()
		   
		   obj.mainthing()
        	   self.jobs.task_done()
        	    #  obj.parser.print_help()
	       except Exception,e: 
        	   #print str(e)
                   #traceback.print_exc()
	           break

      def runchangeheaders(self):
	  self.s12 = 'changeheaders'
	  self.inputo[self.s12] = {}
          self.printer(self.s12)
	  self.chfolder = self.workdir + self.s12 + "/"
	  os.system("mkdir -p " + self.chfolder)

	  self.importlist = []
	  if self.config.has_section(self.s12):
	      self.THREADLIMIT_ch = int(self.config.get(self.s12 ,'-T'))
	      self.config.remove_option(self.s12 ,'-T')
	      if self.config.get(self.s12 ,'-i') == 'DEFAULT' :	         		 
		 for fasta in self.splitpaths:		 
		     folder = self.splitpaths[fasta] 			     	     
		     fastalist = glob.glob(folder + "/*.fsa")
		     self.importlist = self.importlist + fastalist
		     
	         
	      else: 
	         folder = self.config.get(self.s12 ,'-i')
	         self.importlist = glob.glob(folder + "/*.fsa")
	   
	      if self.config.get(self.s12 ,'-n') == 'DEFAULT' :
	         self.config.remove_option(self.s12 ,'-n')
              if self.config.get(self.s12 ,'-ht') == 'DEFAULT' :
	         self.config.remove_option(self.s12 ,'-ht')
		 
	      try:
	        if self.config.get(self.s12,'-run').upper() == "TRUE":
	           self.changeh(self.importlist)		    
              except Exception,e:
	           print str(e)
	           traceback.print_exc()



#####################################################################################################################################
####    															 ####
####    				MAKE MULTIPLE ALIGNMENTS using THREADING     s13					 ####  ???
####    															 ####
#####################################################################################################################################      


      def alignments(self, inputlist, THREAD_LIMIT):
         self.jobs = Queue.Queue(0)
	 self.singlelock = threading.Lock()
         self.target = self.AlignMuscle	 
	 self.StartQueueAndThreads(inputlist, THREAD_LIMIT)	                
	  
      def AlignMuscle(self):
           "method designed to do the actual alignment with Muscle"
	   while 1: 

		try: 
	           infile = self.jobs.get(True,1)
		   inf = basename(infile)		   
		   indir = os.path.dirname(infile)
		   base = ".".join(inf.split(".")[0:-1])
		   outfile =  self.savedir13 + base + ".aln"
	           com = "muscle -in %s -out %s -clwstrict %s -quiet" % (infile, outfile, self.alnopts)
                   os.system(com)
		   
		   self.jobs.task_done()
		except:
		   break
      def runalignments(self):
      
	  self.s13 = 'runalignment'
	  self.inputo[self.s13] = {}
          self.printer(self.s13)
	  self.savedir13 = self.workdir + self.s13 + "/"
	  os.system("mkdir -p " + self.savedir13)
          self.alignlist = []
	  self.namechange = {}
	  self.clustersize = {}
	  
	  if self.config.has_section(self.s13):
	      THREADLIMIT = int(self.config.get(self.s13 ,'-T'))
	      if self.config.has_option(self.s13 ,'-i'):
	         inputfolder = self.config.get(self.s13 ,'-i')		 
		 if inputfolder == 'DEFAULT' :
		    if self.config.has_section(self.s12): 
		       importlist =  glob.glob(self.chfolder + "*.new.fsa")
		       for newfasta in importlist: 
		           #newfasta = fasta[:-4] + ".new.fsa"
			   fastasize = open(newfasta, "r").read().count(">")
			   ht = newfasta[:-8] + ".headertable"
		           self.alignlist.append(newfasta)
			   clusterid = basename(newfasta).split(".")[0]
			   self.namechange[clusterid] = ht
			   self.clustersize[clusterid] = fastasize
		    else: 
		       msg = "NO input folder name or list found!" 
		       self.invokeexit(self.s13,msg)        		 		             
		 else: 
                    importlist = glob.glob(inputfolder + "/*.new.fsa")
		    for newfasta in importlist: 
		           #newfasta = fasta[:-4] + ".new.fsa"
			   fastasize = open(newfasta, "r").read().count(">")
			   ht = newfasta[:-8] + ".headertable"
		           self.alignlist.append(newfasta)
			   clusterid = basename(newfasta).split(".")[0]
			   self.namechange[clusterid] = ht
			   self.clustersize[clusterid] = fastasize
		    if len(importlist) == 0: 
		        msg = "NO input folder name or list found!" 
              		self.invokeexit(self.s13,msg)	  
			
	      if self.config.has_option(self.s13 ,'fast'):
                 fastopts = self.config.get(self.s13 , 'fast')
		 if fastopts.upper() == "TRUE": 
		    self.alnopts = " -maxiters 3"
		 else: 
		    self.alnopts = " -maxiters 8"
		    
	      else:	 
		 msg = "NO alignment options are found, change config file, add an option or keep DEFAULT!" 
                 self.invokeexit(self.s13,msg)	
		 		  

	      try:
	        if self.config.get(self.s13,'-run').upper() == "TRUE":
	           self.alignments(self.alignlist, THREADLIMIT )		    
              except Exception,e:
	           print str(e)
	           traceback.print_exc()     
            
#####################################################################################################################################
####    															 ####
####    				MAKE TREES in PAUP - THREAD for making MULTIPLE TREES s14				 ####  ???
####    															 ####
#####################################################################################################################################      
 
      
      
      def treemaker(self,inputlist,THREAD_LIMIT): #(includes paup & treenames)
         #import runtrees as rt	 
         self.jobs = Queue.Queue(0)
	 self.singlelock = threading.Lock()
         self.target = self.PaupNJ	 
	 self.StartQueueAndThreads(inputlist, THREAD_LIMIT)	                 	 
        

      def PaupNJ(self):           
	   while 1: 
		try: 
		   from os.path import basename
                   import commands, tempfile
		   print "Started PAUPNJ"
	           ALNFILE = self.jobs.get(True,1)		   	           	
        	   basename = basename(ALNFILE)[:-4]
        	   NEXUSFILE = self.savedir14 + basename + ".nexus"
        	   PIRFILE =  self.savedir14 + basename + ".pir"
        	   LOGFILE =  self.savedir14 + basename + ".lognj"
        	   NJFILE =  self.savedir14 + basename + ".njtree"
		   STDOUT =  self.savedir14 + basename + ".stdout"
		   print 'ALNFILE', ALNFILE
		   #self.singlelock.release()
        	   com = "clustalw  -convert -infile=" + ALNFILE + " -output=pir -outfile=" + PIRFILE
        	   os.system(com)

        	   mall = '''tonexus fromfile=PIRFILE tofile = NEXUSFILE format = pir datatype = protein interleave = yes replace = yes;

        	   BEGIN paup;
        	   LOG file=LOGFILE REPLACE;
        	   SET CRITERION=distance;
        	   SET TAXLABELS=full;
        	   bootstrap search=nj nreps=1000 keepall = yes;
        	   contree all /le50 = yes;
        	   savetrees from=1 to=1 FMT=phylip FILE=NJFILE REPLACE=YES brlens=Yes;
        	   [showdist;]
        	   describetrees all /brlens = yes;
        	   quit;
        	   endblock;
        	   '''

        	   for i in ['ALNFILE', 'PIRFILE', 'LOGFILE', 'NJFILE', 'NEXUSFILE']:
        	       print i, eval(i)
        	       mall = string.replace(mall, i, eval(i))

        	   tmpfile = tempfile.mktemp()
        	   fid = open(tmpfile, 'w+')
        	   lines = string.split(mall, '\n')
        	   fid.write("#NEXUS\n" + lines[0] + '\n')

        	   fid.close()
		   os.system("paup %s  -u -n > %s" % (tmpfile, STDOUT))
        	   os.remove(tmpfile)
        	   txt = open(NEXUSFILE).read()
        	   txt = string.replace(txt, basename, "'%s'" % basename)

        	   fid = open(NEXUSFILE, 'a')
        	   fid.write("\n" + mall[len(lines[0]):] + "\n")

        	   fid.close()

		   os.system("paup %s  -u -n >> %s" % (NEXUSFILE, STDOUT))
		   os.remove(PIRFILE)
		   self.jobs.task_done()
		except:
		   traceback.print_exc()
		   break

      #def runpaup(self,ALNFILE): 

	      
      def runtrees(self):
      
	  self.s14 = 'trees'
	  self.inputo[self.s14] = {}
          self.printer(self.s14)
	  self.savedir14 = self.workdir + self.s14 + "/"
	  os.system("mkdir -p " + self.savedir14)
          treelist = []
	  self.treefiles = {}
	  if self.config.has_section(self.s14):
	      THREADLIMIT = int(self.config.get(self.s14 ,'-T'))
	      inputlist =  glob.glob(self.savedir13 + "CL_*.new.aln")
	      #print self.clustersize
	      for i in inputlist: 
	         clusterid = basename(i).split(".")[0]
		 
                 if self.clustersize[clusterid] >  2: 
		     treelist.append(i)   
	      try:
	        if self.config.get(self.s14,'-run').upper() == "TRUE":
	           self.treemaker(treelist, THREADLIMIT)		    
              except Exception,e:
	           print str(e)
	           traceback.print_exc() 
		   
		   
		    	      
	      
#####################################################################################################################################
####    															 ####
####    				               CHANGE TREE NAMES	       						 ####  ???
####    															 ####
#####################################################################################################################################      
    
      def treenamechange(self,inputlist,THREAD_LIMIT): #(includes paup & treenames)               
         self.jobs = Queue.Queue(0)
	 self.singlelock = threading.Lock()
         self.target = self.ChangeTreeNames 
	 self.StartQueueAndThreads(inputlist, THREAD_LIMIT)	                 	 

      def ChangeTreeNames(self):
	 
	 while 1: 
		try: 
		   import treenames 
	           infile = self.jobs.get(True,1)
		   #self.singlelock.acquire()		   	   
		   clusterid = basename(infile).split(".")[0]
		   headerfile = self.namechange[clusterid]
		   output = self.savedir14_1 + clusterid
		   T_obj=treenames.main()
		   obj = self.P_obj
#		   print clusterid		   
#		   print infile
#		   print headerfile
#		   print output
		   obj.updateConfig('treenames', "-t", infile) 
		   obj.updateConfig('treenames', "-i", headerfile)	   
                   obj.updateConfig('treenames', "-o", output )
		   optlist = obj.converttolist('treenames')
		   #self.singlelock.release()		   		   		   		  
		   T_obj.args = T_obj.parser.parse_args(optlist)
		   self.singlelock.acquire()	
		   T_obj.printer = self.printer2
		   T_obj.printer("### OPTIONS : " + str(T_obj.args))
		   self.singlelock.release()		   
		   T_obj.mainthing()		    		   		   		     
		   self.jobs.task_done()
		except:
		   traceback.print_exc() 
		   break

      def runtreenames(self):
          self.s14_1 = 'treenames'
	  self.inputo[self.s14_1] = {}
          self.printer(self.s14_1)
	  self.savedir14_1 = self.workdir + self.s14_1+ "/"
	  os.system("mkdir -p " + self.savedir14_1)
          self.alignlist = []
	  if self.config.has_section(self.s14_1):
	     #print self.config.items(self.s14_1)
	     THREADLIMIT = int(self.config.get(self.s14 ,'-T'))
	     self.P_obj.updateConfig('treenames', "-p", "True" )
	     if self.config.get('treenames', "-u").upper() == "DEFAULT":	        
		self.P_obj.updateConfig('treenames', "-u", self.Poutbase + ".uniprottab" )
	     inputlist =  glob.glob( self.savedir14 + "CL_*.njtree")
             try:
	        if self.config.get(self.s14_1,'-run').upper() == "TRUE":
		   starttime= time.time()
		   self.printer("\n### treenames initialized at " + dt.today().strftime("%d-%m-%Y %H:%M:%S")+ " with THREAD LIMIT: " + str(THREADLIMIT) )	           	     
		   self.treenamechange(inputlist, THREADLIMIT )	
		   timeused = (time.time()-starttime)
	           self.printer("### Time used: "+str(round(timeused)) + " seconds ("+str(round(timeused/60))+ " min)\n") 
             except Exception,e:
	           print str(e)
	           traceback.print_exc() 

#####################################################################################################################################
####    															 ####
####    				              QUEUE AND THREAD STARTER	       						 ####  DONE
####    															 ####
#####################################################################################################################################      



      def StartQueueAndThreads(self,inputlist,THREAD_LIMIT): 
              # Spawn the threads
	      # print "Startqueueandthreads() started"			        
              myrange = xrange(THREAD_LIMIT)	        
              #self.ThreadNo = list(myrange)[::-1]	        
              for x in myrange: 			        
             	  threading.Thread(target=self.target).start()  
              # Put stuff in queue			        
              for i in inputlist: #Can be changed to anything, parsing the input args
                    fullQ = 1 
                    while fullQ != 0:
                           try:
                           	   self.jobs.put(i,block=True,timeout=1)
                           	   fullQ = 0
				   
                           except:
                           	   print "Sleeping"
                                   time.sleep(10)
				   
              self.jobs.join() # This command waits for all threads to finish.
       



#####################################################################################################################################
####    															 ####
####    				  MAINTHING :  WHERE THE MAGIC HAPPENS       						 ####  DONE
####    															 ####
#####################################################################################################################################      





      def mainthing(self):
      
##################################### SECTION 0: BEGINNING ################################## DONE-working	
          self.printer("########################### SECTION 0: BEGINNING ##################################")
          self.printer("Parsing configuration file")
	  self.parseconfig()
	  self.printer("Reading Beginning Section options")	  
	  self.beginning()
	  self.genecatcheck()
	  self.startmain()
	  	  
##################################### SECTION 1: FORMAT UNIPROT ############################# DONE-working
          self.printer("########################### SECTION 1: FORMAT UNIPROT #############################")
          self.printer("Checking Uniprot Input")	            
	  self.uniprotcheck()
	  self.printer("Running FormatSwiss if necesary")
	  self.formatswiss()
	  
##################################### SECTION 2: RUN UBLAST UNIPROT########################## DONE-working
          self.printer("########################### SECTION 2: RUN UBLAST UNIPROT #########################")
	  self.printer("Running UBLAST against Uniprot")
          self.ublast_uniprot()
##################################### SECTION 3: RUN UBLAST GENECATALOG###################### DONE-working
          self.printer("########################### SECTION 3: RUN UBLAST GENECATALOG #####################")
	  self.printer("Running UBLAST against Genecatalog")
          self.ublast_genecatalog()

##################################### SECTION 4: BLASTDECIDE for GC_vs_UN ################### DONE-working
          self.printer("########################### SECTION 4: BLASTDECIDE for GC_vs_UN ###################")
	  self.printer("Running Blastdecide for Genecatalog vs Uniprot")
          self.runblastdecide_uni()
		   
##################################### SECTION 5: BLASTDECIDE for GC_vs_GC ################### DONE-working
          self.printer("########################### SECTION 5: BLASTDECIDE for GC_vs_GC ###################")
	  self.printer("Running Blastdecide for Genecatalog vs Genecatalog")
          self.runblastdecide_gc()

##################################### SECTION 6: Clustergenes ############################### DONE-working
	  self.printer("########################### SECTION 6: Clustergenes ###############################")
	  self.printer("Clustering genes")
	  self.clustergenes()

##################################### SECTION 7: Filter Clusters ############################ DONE-working
          self.printer("########################### SECTION 7: Filter Clusters ############################")
	  self.printer("Filtering clusters")
	  self.clusterFilter()

##################################### SECTION 8: Get uniprot info for clusters ############## DONE-working
	  self.printer("########################### SECTION 8: Get uniprot info for clusters ##############")
	  self.printer("Fetching information for each clusters")
          self.parseclusterinfo()
 
##################################### SECTION 9:  GOATOOLS ################################## DONE-working
          self.printer("########################### SECTION 9:  GOATOOLS ##################################")
	  self.printer("Running enrichments for clusters")
	  self.enrichment()
	  
##################################### SECTION 10: Extract Sequences ######################### DONE-working
          self.printer("########################### SECTION 10:  CLUSTER2FASTA ##################################")
	  self.printer("Running clusters2fasta")   # 3 files max
	  self.runcluster2fasta()
		   
##################################### SECTION 11: Split clusters ############################ DONE-working
	  self.printer("########################### SECTION 11: SPLITCLUSTERS ##################################")
	  self.printer("Splitting clusters")  # 3 files max
	  self.splittclusters()
	  
##################################### SECTION 12: Change headers ############################ DONE-working
	  self.printer("########################### SECTION 12: CHANGEHEADERS ##################################")
          self.printer("Changing headers temporarily")   # 1 file per cluster
	  self.runchangeheaders()
	   

##################################### SECTION 13: Alignments ################################
	  self.printer("########################### SECTION 13: ALIGNMENTS ##################################")
          self.printer("Running alignments for the clusters") # 1 file per cluster
          self.runalignments()

##################################### SECTION 14: Trees #####################################
	  self.printer("########################### SECTION 14: NJ TREES ##################################") 
          self.printer("Making trees for clusters") # 1 file per cluster
	  self.runtrees()
	  self.printer("########################### SECTION 14: CONVERTING TREE NAMES ##################################") 
	  self.runtreenames()
	  self.makeoptslog()
          timeused = (time.time()-self.start)
	  self.printer("### Time used: "+str(round(timeused)) + " seconds ("+str(round(timeused/60))+ " min)\n")	     
               
if __name__ == "__main__":
    try:
        myobj = MGFunc()
        myobj.opts = myobj.parser.parse_args(sys.argv[1:])
        myobj.printer("\n### "+sys.argv[0]+" initialized at "+ myobj.timestarted + "\n")
        myobj.printer("### OPTIONS: "+str(myobj.opts)+"\n")
        myobj.mainthing()
    #except IOError as i:
     #   print "I/O error({0}): {1}".format(i.errno, i.strerror)
    except Exception,e:
        print str(e)
        traceback.print_exc()
	 
	 
#######################################################################
#  TO DO
#
# finish  last two sections
# finish config file with default values - Done
#
# Upgrade
# 
# Print all errors to .runerr, also the exceptions
# Add threading whenver you can
# Maybe change the clustering script with threading
# Add graphical results
# Blast cluster results to an additional database if given : Section X: Compare clusters to EGGNOG or PFAM



#######PARSER REMOVED ITEMS
# self.parser.add_argument("-p", metavar="path", default="./", type=str, help="path to the fasta sequence files")
#self.parser.add_argument("-w", metavar="workdir", default="./MGFunc_results", type=str, help="Directory where all the results will be saved")
#self.parser.add_argument("-u", metavar="usearchpath", default="", type=str, help="Path to the usearch program")	  
#self.parser.add_argument("-mu", "--makedatabase_u", action="store_true" , help="If chosen, usearch ublast make_db will run on uniprot database file")
#self.parser.add_argument("-mg", "--makedatabase_g", action="store_true" , help="If chosen, usearch ublast make_db will run gene catalog fasta file")





