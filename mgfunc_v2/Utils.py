#!/usr/bin/env python
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen (asli@cbs.dtu.dk)
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
import sys, gzip
#import re, string
#import argparse
import os
#from Bio.Blast import NCBIStandalone
#from Bio import SeqIO
from operator import itemgetter, attrgetter
import pickle
#sys.path.append('/home/projects5/pr_53035/people/asli/bin/lib/python2.7/site-packages')
#import networkx as nx
from os.path import basename
import time
from datetime import datetime as dt
import fasta2index
import traceback
import ConfigParser
import sys
import cPickle
from fasta2index import Fasta2Index

#df = {blastdecide.py: inputopts: "-bh"..., outputopts: "-o" , clustergenes: "-i", "-o"..., }

class ParseConfig:

     def __init__(self, configfile):
        self.config = configfile

     def parser(self):
         self.cparser = ConfigParser.ConfigParser()
         self.cparser.read(self.config)

        
	 
     def makedict(self):
         #print "Putting config-file arguments into directory of list\n"
      	 
         dpickle=open("configdict.p","wb")
         self.d = {}
	 for section in self.cparser.sections():
        	 self.d[section] = []
        	 for var,var2 in self.cparser.items(section):                	 
                	 #if var2 == 'False' : do nothing
			 if var2 == 'True':
			    self.d[section].append(var)
			    
			 if var2 != 'False' and var2 != 'True':
			    self.d[section].append(var)
			    self.d[section].append(var2)
			    			   
	 cPickle.dump(self.d,dpickle)
         
	 dpickle.close()
     
     def converttolist(self, section):
         l=[]

	 for var,var2 in self.cparser.items(section):
	     if var != "-run":
	     #if var2 == 'False' : do nothing
		if var2 == 'True':
		   l.append(var)		
		if var2 != 'False' and var2 != 'True':
	           l.append(var)
	           l.append(var2)

	 return l
	 


     def updateConfig(self, sectionname, newopt, newvalue):      
          if not self.cparser.has_option(sectionname, newopt):		 
             self.cparser.set(sectionname, newopt, newvalue)
	  else: 
	     if newvalue != self.cparser.get(sectionname, newopt):
	        self.cparser.set(sectionname, newopt, newvalue)
	     
#example dictionary values for blastdecide_ublast
#as it will look in command-line form
#print " ".join(d["blastdecide_ublast"])



########## have a current input/output 



     def parseconfig_all(self, sectionname, workdir, previoussec, nextsec):
         if self.cparser.has_section(sectionname):
	     savedir = workdir + sectionname
	     os.system("mkdir -p " + savedir)
	     
	     for var,var2 in self.cparser.items(sectionname):	         
		    if var2.upper() == "DEFAULT": 
		       newopt = var
		       newvalue= previoussec
		       self.updateConfig(self, sectionname, newopt, newvalue)

             optlist = self.converttolist(sectionname)

         return optlist
	 
	 
	 	        

class StartMain:
     def __init__(self, genecatalog, workdir):
        self.workdir = workdir
        self.genecatalog = genecatalog
	
     # TO DO : prepare for multiple input files, samples might be in separate files	

     def prepareinput(self):  
        '''A function that would read the input or prepare the sample names dictionary: 
	read the gene catalogue IDs, and get the Sample IDs. 
	'''
		
	if os.path.exists(self.workdir + self.genecatalog + ".sampledb.p"): 
		self.sampledb= pickle.load( open(self.workdir + self.genecatalog + ".sampledb.p", "rb" ) )
	else:
	        self.sampledb = {}
	        
		if os.path.exists(self.workdir + self.genecatalog + ".db_index.p"): 
		   self.db_index= pickle.load( open(self.workdir + self.genecatalog + ".db_index.p", "rb" ) )
		else:		
		   obj=Fasta2Index()	
	           obj.opts = obj.parser.parse_args(["-f", self.genecatalog, "-o", self.workdir])		
	           obj.createDBindex()
                
		
		indexfile= self.workdir + self.genecatalog + ".indexed"
		
	        samplefile =  self.workdir +  self.genecatalog + ".samples"
	        
		fid = open(indexfile, "r")
	        ind = open(samplefile, "w")  
	        line = fid.readline() # line = N_revised_scaffold10520_gene1_1.1.1 \t 456 \t 4345
		try:
	            header= line.split("\t")[0]    # header =  N_revised_scaffold10520_gene1_1.1.1		    
		    sampleid=header.split("_")[0]  # sampleid= N		    
		    geneid = "_".join(header.split("_")[1:]) #geneid = revised_scaffold10520_gene1_1.1.1
		    
		except Exception,e: 
	            print str(e) , "Your gene catalog headers are wrong!! It should look like this: >SampleID_geneid. \n Yours look like this: " , line.split("\t")[0] , "\n"
		    sys.exit()
			
	        while 1:
	            line = fid.readline()

	            if not line or not len(line):
		       if sampleid in self.sampledb:
		    	   self.sampledb[sampleid].append(geneid)
		       else:		   
	            	   self.sampledb[sampleid] = [geneid]
		       ind.write("%s\t%s\n" %(sampleid, geneid))
	               break

  		    if sampleid in self.sampledb:			
  		       self.sampledb[sampleid].append(geneid)		
  		    else:					
  		       self.sampledb[sampleid] = [geneid]

  		    ind.write("%s\t%s\n" %(sampleid, geneid))	 
	            header= line.split("\t")[0]    # header =  N_revised_scaffold10520_gene1_1.1.1		
		    sampleid=header.split("_")[0]  # sampleid= N		
		    geneid = "_".join(header.split("_")[1:]) #geneid = revised_scaffold10520_gene1_1.1.1	       

                ##########################################
		#
		#  N	revised_scaffold10520_gene1_1.1.1		
		#  N	revised_scaffold10520_gene1_2.1.1
		#
		###########################################
		ind.close()
	        fid.close()   
                pdbindex= open(self.workdir + self.genecatalog + ".sampledb.p", "wb" )
	        pickle.dump(self.sampledb, pdbindex)
	        pdbindex.close()

	self.sampleids=self.sampledb.keys()  
	
	
	
class Execute:
      def __init__(self, genecatalog,workdir):
              self.workdir = workdir
              self.genecatalog = genecatalog


    ##############################################################################################################################
      # 															   #
      #							    RUN UBLAST		  					     	   #  DONE
      # 															   #
      ##############################################################################################################################
      
     
      def ublast(self, fasta, db, out):
          
	  self.UBevalue =0.00001	             
	  ublastcmd=self.usearchpath + "usearch -ublast "+ fasta +  " -db " + db + " -evalue " + str(self.UBevalue) + " -blast6out " + out
	  os.system(ublastcmd)
          
	  return 1
	  
      def ublast_makedb(self, fasta, udbname):
      
          udbcmd= self.usearchpath + "usearch -makeudb_ublast " + fasta + " -output " + udbname
	  os.system(udbcmd)
          return 1

      def ublast_genes_vs_uniprot(self,mu):
          self.makedatabase_u = mu
          savedir = self.workdir + "blast/"
	  savedb = self.workdir + "databases/"
	  os.system("mkdir -p " + savedir)
	  os.system("mkdir -p " + savedb)
          self.printer("Running ublast genes catalog against uniprot")	
	  self.uniprotdb_basename = ".".join(basename(self.uniprotfasta).split(".")[0:-1])  
	  if self.makedatabase_u: 	     
	     self.uniprotdb_udb = savedb  + self.uniprotdb_basename + ".udb"  # "/udb_files
	     self.ublast_makedb(self.uniprotfasta, self.uniprotdb_udb)
	  
	  self.GC_vs_UN_ublast = savedir + self.filebase + ".genesvsdb.ublast"  # NAMES could be shortened / results saved to /blast instead
          self.ublast(self.genecatalog, self.uniprotdb_udb, self.GC_vs_UN_ublast)
          return 1
	  
      def ublast_genes_vs_genes(self,mg):
          self.makedatabase_g = mg
          savedir = self.workdir + "blast/"
	  savedb = self.workdir + "databases/"
	  os.system("mkdir -p " + savedir)
	  os.system("mkdir -p " + savedb)
          self.printer("Running ublast gene catalog against itself")
	  self.genecatalog_base = ".".join(basename(self.genecatalog).split(".")[0:-1])
	  if self.makedatabase_g: # or the opposite 
	     
	     self.genecatalog_udb = savedb  + self.genecatalog_base + ".udb"    # "/udb_files
	     self.ublast_makedb(self.genecatalog, self.genecatalog_udb)
          
          
	  self.GC_vs_GC_ublast= savedir + self.filebase + ".allvsall.ublast"
	  self.ublast(self.genecatalog,self.genecatalog_udb, self.GC_vs_GC_ublast)
          return 1
      
      
      
      
      
      
      ##############################################################################################################################
      # 															   #
      #						FORMAT UNIPROT DATABASE FILES							   #
      # 															   #
      ##############################################################################################################################
      
 
      def SwissFasta(self,optslist): 
	  import swiss2fasta
	  try:
              SF_obj = swiss2fasta.main()
              SF_obj.args = SF_obj.parser.parse_args(optslist)
              SF_obj.printer("\n### swiss2fasta initialized at "+ SF_obj.timestarted + "\n")
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
              SP_obj.printer("\n### swiss2tab initialized at "+ SP_obj.timestarted + "\n")
              self.printer("### OPTIONS: "+str(SP_obj.args)+"\n")
              SP_obj.mainthing()
	  except IOError as i:
              print "I/O error({0}): {1}".format(i.errno, i.strerror)
	  except Exception,e:
              print str(e)
              traceback.print_exc()
              sys.exit()
      
      ##############################################################################################################################
      # 															   #
      #						BLASTDECIDE FUNCTIONS								   #  DONE
      # 															   #
      ##############################################################################################################################

      def blastdecide(self, optionslist):                   # abbr. : BD
      
          from blastdecide_ublast import BlastDecision

	  try:    
              BD_obj = BlastDecision()     	      
              BD_obj.opts=BD_obj.parser.parse_args(optionslist)
	      if BD_obj.opts.verbose: 
		 print "\n### blastdecide_ublast.py initialized at " + BD_obj.timestarted  
		 print "### OPTIONS : " + str(BD_obj.opts)
	      BD_obj.mainthing()

	       #  obj.parser.print_help()
	  except Exception,e: 
              print str(e)
              traceback.print_exc()
	      sys.exit()
 
 
      ##############################################################################################################################
      # 															   #
      #						CLUSTERING FUNCTIONS								   #
      # 															   #
      ##############################################################################################################################

      def clustering(self,optionslist):
	  from clustergenes import ClusterGenes
	
	  try:
              CL_obj = ClusterGenes()
	      CL_obj.opts = CL_obj.parser.parse_args(optionslist)
	      if CL_obj.opts.verbose: 
		print "\n### clustergenes.py initialized at " + CL_obj.timestarted  
		print "### OPTIONS : " + str(CL_obj.opts)
              CL_obj.mainthing()
	  #except IOError as i:
    	   #   print "I/O error({0}): {1}".format(i.errno, i.strerror)
	  except Exception,e: 
	      print str(e)
	      traceback.print_exc()
	      sys.exit()	  

      ##############################################################################################################################
      # 															   #
      #						CLUSTERFILTER FUNCTIONS								   #
      # 															   #
      ##############################################################################################################################

      def filter(self, optslist):
      
          import clusterFilter
	  try:	
	      CF_obj = clusterFilter.main()
              CF_obj.args = CF_obj.parser.parse_args(optslist)
              
	      if CF_obj.args.v: 
		 print "\n### clusterFilter.py initialized at " + CF_obj.timestarted  
		 print "### OPTIONS : " + str(CF_obj.args)
              CF_obj.mainthing()
 
	  except IOError as i:
    	      print "I/O error({0}): {1}".format(i.errno, i.strerror)
	  except Exception,e: 
	      print str(e)
	      traceback.print_exc()
	      sys.exit()	  

      ##############################################################################################################################
      # 															   #
      #						GET UNIPROT INFORMATION FOR CLUSTERS						   #
      # 															   #
      ##############################################################################################################################

      def getuniprotinfo(self, optslist):
          from parseclusterinfo import GetUniprotInfo 
	  try:
              M_obj = GetUniprotInfo()
              M_obj.opts = M_obj.parser.parse_args(optslist)
	      if M_obj.opts.verbose: 
		 print "\n### parseclusterinfo.py initialized at " + M_obj.timestarted 
		 print "### OPTIONS : " + str(M_obj.opts)
              M_obj.mainthing()
	  except IOError as i:
	      print "I/O error({0}): {1}".format(i.errno, i.strerror)
	  except Exception,e: 
	      print str(e)
	      traceback.print_exc()
	      sys.exit()
	      
	      
      ##############################################################################################################################
      # 															   #
      #					EXTRACT CLUSTER SEQUENCES FROM FASTA FILE						   #
      # 															   #
      ##############################################################################################################################

      def extract(self, inputlist):
          THREAD_LIMIT = len(inputlist)
	  self.jobs = Queue.Queue(0)
	  self.singlelock = threading.Lock()
          self.target = self.Cluster2Fasta	 
	  self.StartQueueAndThreads(inputlist, THREAD_LIMIT)
          
      def Cluster2Fasta(self):	  
	  import cluster2fasta
	  
	  while 1: 
	   
	       try:
		   
        	   infile = self.jobs.get(True,1)
		     
		   self.singlelock.acquire()
		   self.P_obj.updateConfig("cluster2fasta", "-c", infile)		   
 		   self.P_obj.updateConfig("cluster2fasta", "-o", self.extoutlist[infile])
		   Eopts = self.P_obj.converttolist("cluster2fasta")
		   self.singlelock.release()
		   
		   
        	   E_obj = cluster2fasta.main()
        	   E_obj.args = E_obj.parser.parse_args(Eopts)
        	   E_obj.printer("\n### cluster2fasta.py initialized at "+ E_obj.timestarted + "\n")
        	   self.printer("### OPTIONS: "+str(E_obj.args)+"\n")
        	   E_obj.mainthing()
		   os.system("cat " + outfile + ".koalagenes.fasta " + outfile + ".uniprotids.fasta > " + outfile + ".fasta" )
		   self.jobs.task_done()
	       except IOError as i:
        	   print "I/O error({0}): {1}".format(i.errno, i.strerror)
	       #except Exception,e:
        	   #if str(e) == "Empty":
		    #  sys.exit()
	       except: 
	           break      
		      
        	 #  traceback.print_exc()
		  # sys.exit()


      ##############################################################################################################################
      # 															   #
      #						SPLIT CLUSTERS FASTA FILES			 				   #  !!!!!!!!!!!!!!!!
      # 															   #
      ##############################################################################################################################

      def splitc(self, inputdict):
          inputlist = inputdict.keys()
	  self.splitdict = inputdict
          THREAD_LIMIT = len(inputlist)
	  self.jobs = Queue.Queue(0)
	  self.singlelock = threading.Lock()
          self.target = self.SplitClusters	 
	  self.StartQueueAndThreads(inputlist, THREAD_LIMIT)
	  
	  
      def getclusterids_fromfasta(self, filename):
          os.system('grep ">" ' + filename + ' | cut -f1 -d":" | tr -d ">" > ' + filename[0:-4] + '.fasta.clid'  )

      def getclusterids_fromCLfile(self, filename):
          os.system("grep -P '\tN\t' " + filename + " > " + filename[0:-4] + ".nohits.txt")
	  os.system("grep -P '\tS\t' " + filename + " > " + filename[0:-4] + ".semihits.txt")
	  os.system("grep -P '\tA\t' " + filename + " > " + filename[0:-4] + ".allhits.txt")
	  os.system("cut -f1 " + filename[0:-4] + ".nohits.txt > " + filename[0:-4] + ".nohits.clid"
	  os.system("cut -f1 " + filename[0:-4] + ".semihits.txt > " + filename[0:-4] + ".semihits.clid"
	  os.system("cut -f1 " + filename[0:-4] + ".allhits.txt > " + filename[0:-4] + ".allhits.clid"
	  
	  
	  
      def SplitClusters(self):
           self.splitpaths = []
           from splitclusters import Split 
	   while 1: 

		  try:    
		         infile = self.jobs.get(True,1)
			 
			 self.singlelock.acquire()
		         self.P_obj.updateConfig("splitclusters", "-i", infile)			   	      
                         self.P_obj.updateConfig("splitclusters", "-cl", self.splitdict[infile])
			 
		         if "nohits" in infile: 
	                     path = self.workdir + "nohits"
		         elif "semihits" in infile: 
	                     path = self.workdir + "semihits"		             		 
		         elif "allhits" in infile: 
	                     path = self.workdir + "allhits"		             		 
		         else: 
	                     path = self.workdir + "splitcluster"	     

                         os.system("mkdir -p " + path )
		         self.P_obj.updateConfig("splitclusters", "-o", path)
			 self.splitpaths.append(path)
		         SPopts = self.P_obj.converttolist("splitclusters")
			 self.singlelock.release()
			 
        	         obj = Split()        
        	         obj.opts=obj.parser.parse_args(SPopts)   # mode = splitclusters or changeheaders
        	         if obj.opts.verbose: 
        	            print "\n### splitclusters initialized at " + obj.timestarted  
        	            print "### OPTIONS : " + str(obj.opts)
        	         obj.mainthing()
                         self.jobs.task_done()
			 
			 
		  except Exception,e: 
        	      print str(e)
        	      traceback.print_exc()
		      sys.exit()

      ##############################################################################################################################
      # 															   #
      #						CHANGE FASTA HEADERS								   # !!!!!!!!!!!!!!!!!!!
      # 															   #
      ##############################################################################################################################
      def changeh(self,inputlist):
         THREAD_LIMIT = self.THREADLIMIT_ch
         self.jobs = Queue.Queue(0)
	 self.singlelock = threading.Lock()
         self.target = self.changeheaders_outside	 
	 self.StartQueueAndThreads(inputlist, THREAD_LIMIT)
      
      def changeheaders_outside(self):
          while 1: 
	       try:    
		    
		   infile = self.jobs.get(True,1)
		   
		   self.singlelock.acquire()
		   self.P_obj.updateConfig("splitclusters", "-i", infile)		   		
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
		   self.singlelock.acquire()
		   self.P_obj.updateConfig("splitclusters", "-i", infile)		   		
		   opts = self.P_obj.converttolist('changeheaders')
		   self.singlelock.release()
	      
        	   obj = Change()        
        	   obj.opts=obj.parser.parse_args(opts)   # mode = splitclusters or changeheaders
        	   if obj.opts.verbose: 
        	      print "\n### changeheaders initialized at " + obj.timestarted  
        	      print "### OPTIONS : " + str(obj.opts)
        	   obj.mainthing(sys.argv[1:])
        	   self.jobs.task_done()
        	    #  obj.parser.print_help()
	       except Exception,e: 
        	   print str(e)
                   traceback.print_exc()
	           sys.exit()





      ##############################################################################################################################
      # 															   #
      #						MAKE MULTIPLE ALIGNMENTS using THREADING					   #
      # 															   #
      ##############################################################################################################################


      def alignments(self):
               
	 inputlist =  glob.glob(self.workdir + "/CL_*.new.fsa")	 
	 self.target = self.AlignMuscle	 
	 self.StartQueueAndThreads(inputlist)
	  
      def AlignMuscle(self, options):
           "method designed to do the actual alignment with Muscle"
           base = ".".join(infile.split(".")[0:-1])

	   while 1: 

		try: 
	           infile = self.jobs.get(True,1)
		   outfile = base + ".aln"
	           com = "muscle -in %s -out %s -clwstrict %s" % (infile, outfile, options)
                   os.system(com)
		   self.jobs.task_done()
		except:
		   break
      
            
      ##############################################################################################################################
      # 															   #
      #					MAKE TREES in PAUP - THREAD for making MULTIPLE TREES 		                           #
      # 															   #
      ##############################################################################################################################
      
      
      def treemaker(self): #(includes paup & treenames)
         import runtrees as rt
         
	 inputlist =  glob.glob(self.workdir + "/CL_*.new.aln")	 
	 self.target = ra.PaupNJ	 
	 self.StartQueueAndThreads(inputlist)
      
      
      def treenamechange(self): #(includes paup & treenames)
               
	 inputlist =  glob.glob(self.workdir + "/CL_*.tree") 	 
	 self.target = self.ChangeTreeNames	 
	 self.StartQueueAndThreads(inputlist)
      

      def PaupNJ(self, options): 
           "method designed to do the actual alignment with Muscle"

	   while 1: 
		try: 
	           infile = self.jobs.get(True,1)
		   
	           com = "python2.7 bb_phylo.py %s  " % (infile)
                   os.system(com)
		   self.jobs.task_done()
		except:
		   break



      ##############################################################################################################################
      # 															   #
      #						CHANGE TREE NAMES								   #
      # 															   #
      ##############################################################################################################################

      def ChangeTreeNames(self, options):
         import treenames2 
	
         optlist = self.P_obj.converttolist("treenames")
	 
	 while 1: 

		try: 
	           infile = self.jobs.get(True,1)
		   base = ".".join(infile.split(".")[0:-1])
		   headerfile = base + ".headertable"
		   T_obj=treenames2.main()
                   optlist = optlist + ["-t", infile] + ["-i", headerfile] 
		   T_obj.args = T_obj.parser.parse_args(optlist)
		   T_obj.mainthing()
		   
		   self.jobs.task_done()
		except:
		   break
		   
      ##############################################################################################################################
      # 															   #
      #						QUEUE AND THREAD STARTER							   #
      # 															   #
      ##############################################################################################################################


      def StartQueueAndThreads(self,inputlist,THREAD_LIMIT): 
               # Spawn the threads			        
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
                                   time.sleep(20)
				   
              self.jobs.join() # This command waits for all threads to finish.
       
