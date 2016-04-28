#!/usr/bin/env python2.7
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen and Kosai Al Nakeeb 
# Email: sosayweall.github@gmail.com
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
import sys, gzip
import re, string
import argparse
import os
#from Bio.Blast import NCBIStandalone
from operator import itemgetter, attrgetter
#sys.path.append('/usr/bin/lib/python2.7/site-packages')
from os.path import basename
import time
from datetime import datetime as dt
import pickle
import subprocess
import glob
import threading,Queue
import traceback
usage= sys.argv[0] + " [-h] [-c C] [-d D] [-b B] [-bh] [-i] [-v]\n"
helpstr= '''
description: 
This script identifies the uniprot hits for each cluster.
If the option -bh is given, only outputs best hit for a gene in a cluster.\n
First an index file for the genes_vs_database.blastdecide file is created in order 
to fasten the process.
'''

epi="Author: Asli I. Ozen and Kosai Al-Nakeeb (email: sosayweall.github@gmail.com)" 

class GetUniprotInfo:
    def __init__(self):
        self.timing=""
	self.start = time.time()
	self.d_ = dt.today()
	self.timestarted = self.d_.strftime("%d-%m-%Y %H:%M:%S")
        self.parseArgs()
	self.cluster = {}
	self.nocluster = {}
	self.clfile={}
	self.preclust = {}
	self.GO = {}
	self.EC = {}
	self.KEGG = {}

    def printer(self,string): 
        if self.opts.verbose:
            print string,	
    def parseArgs(self):
        self.parser = argparse.ArgumentParser(description=helpstr, epilog=epi)
        #parser.add_option("-p", "--path", dest="P", default="./", type=str, help="path to the blast output files")
        #parser.add_argument("-o", "--opath", dest="OP", type="str", help="path to the parse output files")    	
        
        #self.parser.add_argument("-d", "--uniprotdb", type=str, help="uniprot database tab separated format",required=True)
	self.parser.add_argument("-c", metavar="clusterfile", type=str, help="cluster file",nargs=1,required=True)
        self.parser.add_argument("-b", metavar="blastdecide", type=str, help="blastdecide file",nargs='*',required=True)
	self.parser.add_argument("-ut", metavar="uniprottab", type=str, help="uniprot tab formatted data file",required=True)
	self.parser.add_argument("-ui", metavar="uniprotindex", type=str, help="uniprot tab index data file")
	self.parser.add_argument("-g", metavar ="allgeneids", type=str, help="Allgeneids used in the study",nargs=1)
	self.parser.add_argument("-uif", metavar="uniprotindex", type=str, help="folder name uniprot tab index data files. Files must end with '.tab.index'")
	self.parser.add_argument("-bh", "--besthit", action = "store_true", help="If chosen only best uniprot hits will be output")
	self.parser.add_argument("-na", "--addna", action = "store_true", help="Include results for unknown genes")
	self.parser.add_argument("-e", metavar ="expand", help="Include the actual size of each gene from precluster file (homology reduction)",nargs=1)
      	self.parser.add_argument("-i", "--index", action = "store_true", help="If chosen, blastdecide file is indexed first.")
	self.parser.add_argument("-t",type= int, help="Thread limit",default=2)
	self.parser.add_argument("-o", metavar="outputbase", type=str, help="output FILE basename. results will be saved as \
	                                                                     'basename.extention'. Path name can be included")									             

	self.parser.add_argument("-v", "--verbose", action="store_true" , help="increase output verbosity")
        return 1


    def StartQueueAndThreads_1(self):
	    myrange =  range(self.THREAD_LIMIT)
	    self.ThreadNo = list(myrange[::-1])	    
            for x in myrange:
                    threading.Thread(target=self.mythread).start()           	    
	    self.jobs1.join()	    

    def mythread(self):
            import traceback
            therandom = self.ThreadNo.pop() + 1
	    self.singlelock.acquire()
	    self.printer("starting thread number:" + str(therandom))
	    self.singlelock.release()
	    while True:
		    try:
			    clid = self.jobs1.get(True,1)
			    self.singlelock.acquire()
			    #print "Working on " + clid
			    self.singlelock.release()
												        
			    a = self.getblastres4cluster(clid)  				        
			    if a==1 :								        
			       b = self.fetchinfo4cluster(clid) 				        
			    else:								        
			       self.printer("There was an error with getting uniprot id for cluster : " + clid )
			       b=0
			       #self.jobs1.task_done()						        

			    if b == 1: 
			       self.writeresults(clid)
			       self.jobs1.task_done()
			    else: 
			       self.printer("There was an error with fetching uniprot info : " + clid)
			       self.jobs1.task_done()

		    except Exception,e: 
			   print str(e)
			   import traceback
			   traceback.print_exc()
			   break
    def StartQueueAndThreads_2(self):
	    myrange =  range(self.THREAD_LIMIT)
	    self.ThreadNo = list(myrange[::-1])	    
            for x in myrange:
                    threading.Thread(target=self.mythread2).start()	    
	    self.jobs2.join()
    def mythread2(self):
            import traceback
            therandom = self.ThreadNo.pop() + 1
	    self.singlelock.acquire()
	    self.printer("starting thread number:" + str(therandom))
	    self.singlelock.release()
	    while True:
		    try:
			    
			    geneid = self.jobs2.get(True,1)
			    a = self.getblastres4gene(geneid)
			    if a==1 : 
			       b = self.fetchinfo4gene(geneid)
			    else: 
			       self.printer("There was an error with getting uniprot id for cluster" + geneid)
			       self.jobs2.task_done()
			       break			    
			    if b == 1: 
			       self.writegeneresults(geneid)
			       self.jobs2.task_done()
			    else: 
			       self.printer("There was an error with fetching uniprot info" + geneid)
			       self.jobs2.task_done()
		   
		    except Exception,e: 
			   print str(e)   
			   traceback.print_exc()
			   break
    def gzipopen(self,fileID):
            if fileID[-3:] == ".gz":
                    return gzip.open(fileID)
            else:
                    return open(fileID,"rU")    
    def grep(self,word, filename):
	command = "grep -m1 -w '" +  word + "' " + filename 
	#print command
	process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True)
        
	#fout = os.popen(command)
	if process:
	   output = process.communicate()
	   #print output	   
	   if output[0]: 
	      return output[0].strip("\n")
	   else:
	      return None
	      
    def readincluster(self):
            print "Reading clusters"
            self.clgeneids = []
	    cid = self.gzipopen(self.opts.c[0])
	    for line in cid.readlines():
	        if len(line.split("\t")) > 2:
		   linelist=line.strip("\n").split("\t")  					       
		else:
		   linelist=line.strip("\n").split("\t")[0:2]	    

		clid=linelist[0]
		geneids=linelist[1].split(",")
		self.clgeneids = self.clgeneids + geneids
		self.clfile[clid]=geneids
		self.cluster[clid] = {}	
		#print clid	
		self.jobs1.put(clid,block=True,timeout=1)

    def genesoutcluster(self):
        print "Getting gene ids outside clusters"
        self.singeneids=[]
        for i in self.allgeneids:
	    if not i in self.clgeneids:
    		self.singeneids.append(i)
		self.jobs2.put(i,block=True,timeout=1)
        print "Single gene ids length: " , len(self.singeneids), "\n" 
    def createbdecide_index(self):
        indexlist= {}
	for blastdecide in self.blastdecide:
	    posdict={}
	    sampleid= basename(blastdecide).split(".")[0]
	    nh= open(blastdecide + ".index", "w") 
            fh = self.gzipopen(blastdecide)
	    genepos = fh.tell()
	    line=fh.readline()
	    genename = line.split("\t")[0]  
	    posdict[genename] = [str(genepos)]	    
	    while 1:
	        genepos = fh.tell()        	
        	line=fh.readline()
		if not line or not len(line): 
		   break
        	genename = line.split("\t")[0]                  
        	if genename in posdict:
		   posdict[genename].append(str(genepos))
		else:
		   posdict[genename] = [str(genepos)]
	    for gene in posdict:
        	s=posdict[gene]
        	positions=",".join(s)
        	nh.write("%s\t%s\n" %(gene,positions))
	    fh.close()
	    nh.close()
	    indexlist[sampleid] = str(blastdecide) + ".index"
	return indexlist

#          def indextab(self):
#                  fid = open(self.args.o[0],"r")
#                  fout = open(self.args.o[0]+".indexed","w")
#                  line = fid.readline()
#                  while 1:
#                          start = fid.tell()
#                          line = fid.readline()
#                          if not line or not len(line):
#  #                               stop = fid.tell()
#  #                               header = line.split("\t")[0]
#  #                               fout.write(header + "\t" + str(start) + "," + str(stop)+"\n")
#                                  break
#                          stop = fid.tell()
#                          header = line.split("\t")[0]
#                          fout.write(header + "\t" + str(start) + "," + str(stop)+"\n")
#                  fout.close()
#                  fid.close()
#

    def uniprottab_indexes(self):
        uindexes = glob.glob(self.utif + "/*.tab.index")
	self.uniprottabindex = {}
	for i in uindexes: 
	    first= basename(i)[0]
	    self.uniprottabindex[first] = i
	
	
	
    def geneexpansion(self):    
	for i in file(self.prcl):
	    clid = i.split("\t")[0]
	    geneid = i.split("\t")[1]
	    centroid = i.split("\t")[2].strip("\n")
	    if centroid == "*":	       
	       self.preclust[geneid] = [centroid]
	    else: 
	       self.preclust[centroid].append(geneid)
	    
	    
    def getactualsize(self,genelist):
        sum=0
        
	for g in genelist:
	    if g in self.preclust:
	       size = len(self.preclust[g])
	       sum+=size
	    else:
	       size=1
	       sum+=size
	      
	return sum

    def getblastres4cluster(self,clusterid):
        #fh = gzipopen(self.blastdecide)
	try: 							    
	    cf = open(self.outputbase + ".hitcounts.txt", "a")
	    nh = open(self.outputbase + ".cl_nohitslist.txt", "a")				        		 
	    u=[]
	    nohits=[]      
	    idswithhits=[]				       
	    members_names=self.clfile[clusterid]
	    clustersize=len(members_names)		
            if self.opts.e : sizecluster = self.getactualsize(members_names)			
	    sizehits = 0
	    sizenohits = 0

	    for i in members_names:
		sampleid = i.split("_")[0]
        	indexline = self.grep(i,self.indexfile[sampleid])  ### GREP FROM BLASTDECIDE INDEX FILE
		uniprotids=[]    
		l=[]
		if indexline:
		   positions= indexline.split("\t")[1].split(",")							       
        	   idswithhits.append(i)
		   if self.opts.e: 
	              if i in self.preclust:			     
	        	 sizehits+= len(self.preclust[i])
	              else:
			 sizehits+= 1

		   for p in positions:	          
		      self.blasthfiles[sampleid].seek(int(p))                 ### SEEK in BLASTDECIDE FILE
                      blasthit = self.blasthfiles[sampleid].readline()  
		                 
	              try:
			 uniprotid = blasthit.split("\t")[1]  
	        	 score = float(blasthit.split("\t")[4])
	              except: 
		         traceback.print_exc()
			 self.singlelock.acquire()
			 self.printer("Problematic blasthit: " + blasthit )
			 self.printer(self.blastdfiles[sampleid] )		
			 self.singlelock.release()	       
    	              uniprotids.append([uniprotid,float(score)])

		   l = sorted(uniprotids, key=lambda item: float(item[1]),reverse=True) 
		   #sorted(self.matrix.keys(), key=lambda item: (int(item.partition('_')[1]))
		else:
		   nohits.append(i)
		   self.singlelock.acquire()		   
		   nh.write("%s\t%s\n" %(clusterid, i))
		   self.singlelock.release()
		   if self.opts.e: 
	              if i in self.preclust:			     
	        	 sizenohits+= len(self.preclust[i])
	           else:
		      if i in self.preclust:
			 sizenohits+= 1

		if self.opts.besthit:
		    if len(l) > 0:
		       self.cluster[clusterid][i] = [l[0]]
	            else: self.cluster[clusterid][i] = l
		else:		    
		    self.cluster[clusterid][i] = l 

            self.singlelock.acquire()
	    
	    if len(nohits) == clustersize:
	       flag = "N"
	       self.cluster[clusterid]["flag"] = flag
	       if self.opts.e:		   
		  cf.write("%s %s %s %s %s\n" %(clusterid,  flag, sizehits,sizenohits, sizecluster))
	       else:		   
        	  cf.write("%s %s %s %s %s\n" %(clusterid,  flag, len(idswithhits),len(nohits), clustersize))  


	    if len(idswithhits) >= clustersize:
	       flag = "A"
	       self.cluster[clusterid]["flag"] = flag
	       if self.opts.e:		   
		  cf.write("%s %s %s %s %s\n" %(clusterid,  flag, sizehits,sizenohits, sizecluster))
	       else:
		  cf.write("%s %s %s %s %s\n" %(clusterid,  flag, len(idswithhits),len(nohits), clustersize))
	    else:
               if len(idswithhits) > 0:
		  flag = "S"
		  self.cluster[clusterid]["flag"] = flag
        	  if self.opts.e:	      
		     cf.write("%s %s %s %s %s\n" %(clusterid,  flag, sizehits,sizenohits, sizecluster))
		  else:
        	     cf.write("%s %s %s %s %s\n" %(clusterid,  flag, len(idswithhits),len(nohits), clustersize))

            self.singlelock.release()	
	    cf.close()
	    #fh.close() 					        
            return 1
	except:   
	     import traceback
	     traceback.print_exc()
	
	
	
    def fetchinfo4cluster(self,keys):
	try:
            if self.opts.e:	   
	           self.expandclusters(keys)
            if self.cluster[keys]["flag"] != "N":
	       for g in self.cluster[keys]:
        	   if g != "flag":
#        	      if self.opts.e:	   
#	        	     if i in self.preclust:			     
#	        		sizehits+= len(self.preclust[g])
#	        	     else:
#				sizehits+= 1		   		   
		      for l in self.cluster[keys][g]:
                	  if len(l) > 0:
			     uid= l[0].split("|")[2]
			     if type(self.uniprottabindex) == dict:	
				positions = self.grep(uid, self.uniprottabindex[uid[0]])
		             else: 
				positions = self.grep(uid, self.uniprottabindex)
			    
			     if positions:
				start= positions.split("\t")[1].split(",")[0]				
				self.singlelock.acquire()
				self.uh.seek(int(start))
				tabinfo = self.uh.readline()
				with open(self.outputbase + ".uniprottab", "a") as fhw:
				   fhw.write(keys + "\t" + g + "\t" + tabinfo)
				self.GO[uid] = tabinfo.split("\t")[-2]
				self.KEGG[uid] = tabinfo.split("\t")[-1]
				self.singlelock.release()
				description = tabinfo.split("\t")[2]
				m = re.search(r"EC=((\d+|\-)\.)*(\d+|\-)*", description)
				if m: 
				   ecnumber = m.group()
				   self.EC[uid] = ecnumber
                        	else:
				   self.EC[uid] = "N/A"
           
	    return 1
	
	except:   
	     import traceback
	     traceback.print_exc()
	


    def writeresults(self, clusterid):
	try:

	    GN,UN,EC,GO,KE=[],[],[],[],[]		
	    for gene in self.cluster[clusterid].keys():
		   if gene != "flag":	       
  			 GN.append(gene)					########## ADD gene ids for the clusters ####
	    		 l =  self.cluster[clusterid][gene]
			 if len(l) > 0:

        		    ########## ADD Uniprot IDs for the clusters. if besthit: one ID/gene ####

	    		    for d in l:
	    			u = d[0].split("|")[2]
	    			UN.append(d[0])
	    			go = self.getitem(self.GO, u)
				self.singlelock.acquire()
	    			if len(go) > 0 : GO.append(go)				
	    			self.singlelock.release()
				ec = self.getitem(self.EC, u)
	    			if len(ec) > 0 : EC.append(ec)
	    			kegg = self.getitem(self.KEGG, u)
	    			if len(kegg) > 0 : KE.append(kegg.strip("\n"))

	    		 else:  
	    		    if self.opts.addna: UN.append("N/A")					 
            GNall = self.getjoinedlist(GN)
	    UNall = self.getjoinedlist(UN)
	    GOall = self.getjoinedlist(GO)
	    ECall = self.getjoinedlist(EC)
	    KEall = self.getjoinedlist(KE)
	    self.singlelock.acquire()
	    with open(self.outputbase + ".allinfo.txt","a") as ah:
	       ah.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (clusterid, self.cluster[clusterid]["flag"], ",".join(GN), ",".join(list(set(UN))),",".join(GO),",".join(KE),",".join(EC)))
            if self.cluster[clusterid]["flag"] != "N" : 
	       with open(self.outputbase + ".allwithhits.txt","a") as awh:       
	          awh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (clusterid, self.cluster[clusterid]["flag"], ",".join(GN), ",".join(UN),",".join(GO),",".join(KE),",".join(EC)))	      	       
	       with open(self.outputbase + ".allwithGO.txt","a") as awg: 
	          GOlist= ",".join(GO)
	          m = re.search(r"GO:", GOlist)
		  if m:     
	             awg.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (clusterid, self.cluster[clusterid]["flag"], ",".join(GN), ",".join(UN),",".join(GO),",".join(KE),",".join(EC)))
		  
            self.singlelock.release()
#	    awh.close()  
#	    ah.close()  
	except:   
	     import traceback
	     traceback.print_exc()
	

	
    def expandclusters(self,clusterid):
  	#for clusterid in self.cluster.keys():
	    for gene in self.cluster[clusterid].keys():
		   if gene != "flag":		   
	        	 if gene in self.preclust:
			     for members in self.preclust[gene]:
			         if members != "*":
				    newgeneid = members
				    self.cluster[clusterid][newgeneid] = self.cluster[clusterid][gene]  

########################################################## GET RESULTS FOR GENE #######################################

    def getblastres4gene(self,geneid):
	try:
              
              #fh = gzipopen(self.blastdecide)
	      nohits=[]      
	      idswithhits=[]
	      sampleid = geneid.split("_")[0]
              indexline = self.grep(geneid,self.indexfile[sampleid])  ### GREP FROM BLASTDECIDE INDEX FILE
	      uniprotids=[]
	      l=[]
	      sizehits = 0
	      sizenohits = 0
	      #print geneid
	      if indexline:
		 positions= indexline.split("\t")[1].split(",")			    
		 idswithhits.append(geneid)
		 for p in positions:
		    self.blasthfiles[sampleid].seek(int(p))		      ### SEEK in BLASTDECIDE FILE
        	    blasthit = self.blasthfiles[sampleid].readline()
		    #print blasthit
	            try:
		       uniprotid = blasthit.split("\t")[1]  
	               score = float(blasthit.split("\t")[4])
	            except: 
		       traceback.print_exc()
		       self.singlelock.acquire()
		       self.printer("Problematic blasthit: " + blasthit )
		       self.printer(self.blastdfiles[sampleid] )	      
		       self.singlelock.release()	     
                    uniprotids.append([uniprotid,int(score)])
		 l = sorted(uniprotids, key=lambda item: int(item[1]),reverse=True) 
		 #sorted(self.matrix.keys(), key=lambda item: (int(item.partition('_')[1]))
	      else:
		 nohits.append(geneid)
		 if self.opts.e: 
		    if geneid in self.preclust:				 
		       sizenohits+= len(self.preclust[geneid])
		 else:
		       sizenohits+= 1
		 
		 self.singlelock.acquire()
		 with open(self.outputbase + ".single_nohitslist.txt", "a") as nh:
		      nh.write("%s\t%s\n" % (geneid,str(sizenohits)))
		 self.singlelock.release()

	      if self.opts.besthit:
		  if len(l) > 0:
		     self.nocluster[geneid] = [l[0]]
		  else: self.nocluster[geneid] = l
	      else:  		 
		  self.nocluster[geneid] = l 

	      #fh.close()
	      #nh.close()
	      return 1	
	except:   
	     import traceback
	     traceback.print_exc()
	
    def fetchinfo4gene(self,geneid):
	try: 
	     sizehits=0	     
	     if self.opts.e: 
	     	if geneid in self.preclust:			     
	     	   sizehits+= len(self.preclust[geneid])
	     	else:
	     	   sizehits+= 1	   
	        self.expandgenes(geneid)

             for l in self.nocluster[geneid]:
        	 if len(l) > 0:
		    uid= l[0].split("|")[2]
		    
		    if type(self.uniprottabindex) == dict:   
	               positions = self.grep(uid, self.uniprottabindex[uid[0]])
        	    else: 
	               positions = self.grep(uid, self.uniprottabindex)
		    if positions:
	               start= positions.split("\t")[1].split(",")[0]	               
		       self.singlelock.acquire()
		       self.uh.seek(int(start))
	               tabinfo = self.uh.readline()	               
	               with open(self.outputbase + ".uniprottab.singles.txt", "a") as fhw:
	        	  fhw.write(geneid + "\t" + str(sizehits) + "\t" + tabinfo) 		       	               
		       self.GO[uid] = tabinfo.split("\t")[-2]
	               self.KEGG[uid] = tabinfo.split("\t")[-1].strip("\n")		       
	               description = tabinfo.split("\t")[2]
		       m = re.search(r"EC=((\d|\-)\.)*(\d|\-)*", description)
	               if m: 
	        	  ecnumber = m.group()
	        	  self.EC[uid] = ecnumber
                	  #fec.write("%s\t" % self.EC[uid])
        	       else:
	        	  self.EC[uid] = "N/A"
                       self.singlelock.release()
	     
	     return 1
    
	except:   
	     import traceback
	     traceback.print_exc()
	
    def expandgenes(self,gene):    
	if gene in self.preclust:
	    for members in self.preclust[gene]:
		if members != "*":
		   newgeneid = members
		   self.nocluster[newgeneid] = self.nocluster[gene]
 
    def writegeneresults(self, gene):
	try: 
	    ah = open(self.outputbase + ".single_allinfo.txt", "a")
	    GN,UN,EC,GO,KE=[],[],[],[],[]		


	    l =  self.nocluster[gene]
	    if len(l) > 0:
	       for d in l:
		   u = d[0].split("|")[2]
		   UN.append(d[0])
		   go = self.getitem(self.GO, u)
		   self.singlelock.acquire()
		   if len(go) > 0 : 
		      GO.append(go)
		   self.singlelock.release()
		   ec = self.getitem(self.EC, u)
		   if len(ec) > 0 : EC.append(ec)

		   kegg = self.getitem(self.KEGG, u)
		   if len(kegg) > 0 : KE.append(kegg.strip("\n"))

	    else:  
	       if self.opts.addna: UN.append("N/A") 				    

	    UNall = self.getjoinedlist(UN)
	    GOall = self.getjoinedlist(GO)
	    ECall = self.getjoinedlist(EC)
	    KEall = self.getjoinedlist(KE)
	    self.singlelock.acquire()
	    ah.write("SINGLETON\t%s\t%s\t%s\t%s\t%s\n" % (gene,",".join(UN),",".join(GO),",".join(KE),",".join(EC)))
	    self.singlelock.release()
	    ah.close()  
	except:   
	     import traceback
	     traceback.print_exc()

 
####################################################### UTILS #################################    	

    def getitem(self,dic,item):       
        if item in dic:
	   res = dic[item]
	   if res == "N/A":
	     if self.opts.addna: 
	         res = dic[item]
             else: 
	         res = ''     	   
	else: 
	   if self.opts.addna: 
	      res="N/A"
	   else: 
	      res=''
	   
        return res
	  

    def print_files(self):
        aso = open(self.outputbase + ".association", "w")
	popu = open(self.outputbase + ".population", "w")
        for clusterid in sorted(self.cluster.keys()):
	     #self.printer("Printing for " + clusterid + "\n")
	     for gene in self.cluster[clusterid]:
	         if gene != "flag":	    
		    if len(self.cluster[clusterid][gene]) > 0:      
	  	        d = self.cluster[clusterid][gene][0]		    		 
			uid = d[0].split("|")[2]			
			go = self.getitem(self.GO, uid)	
			#print gene, go		
	        	if len(go) > 0: 
			   aso.write("%s\t%s\n" % (gene, go))
			   popu.write("%s\n" % gene)	             
		    else: 
		        if self.opts.addna:
	  	           aso.write("%s\tN/A\n" % (gene))
			   popu.write("%s\n" % gene)
	aso.close()
	popu.close()


    def getjoinedlist(self, idlist):
        alllist=[]
	for i in idlist:
	    m = re.search(r",", i)
	    if m:
	       a=i.split(",")
	       for f in a:
	           alllist.append(f)
            else:
	       alllist.append(i)	       	
	return ",".join(alllist)

    def printer(self,mystring):
	if self.opts.verbose: print mystring

	
    def cleanup(self,filename):
        if os.path.getsize(filename) == 0:
	   os.system("rm " + filename)
	   
	   
###########################################################################################################
#			MAIN FUNCTION
#
###########################################################################################################
		     
    def mainthing(self):
    
        ################### EXPANSION : PRECLUSTERFILE #######################
        if self.opts.e : 
	   self.prcl = self.opts.e[0]	   
	   self.geneexpansion()

	########################## CLUSTERFILE ###############################
	
	self.clusterfile=self.opts.c[0]
	
	##################### UNIPROT FILES ##################################
	self.uniprottab = self.opts.ut
	self.uh=open(self.uniprottab,"r")
	if self.opts.ui:
	   self.uniprottabindex = self.opts.ui	   
	elif self.opts.uif:
	   self.utif=self.opts.uif  # If there are multiple idnex folders, specify, uniprot tab index folder
	   self.uniprottab_indexes()
	else: 
	   sys.exit("Uniprot tab index file(s) could not be found")
	
	########################## BLASTDECIDE ###############################
	
	self.blastdecide = [] 
	self.blastdfiles = {}
	self.blasthfiles = {}
	
	
	if len(self.opts.b) == 1:	    
	    self.opts.b = self.opts.b[0].split(" ")
	
	for i in range(0,len(self.opts.b)): 
	    base = basename(self.opts.b[i])
	    self.blastdecide.append(self.opts.b[i])
	    self.blastdfiles[base.split(".")[0]] = self.opts.b[i]
	    self.blasthfiles[base.split(".")[0]] = open(self.opts.b[i], "r")
	    
	##################### BLASTDECIDE INDEX ##############################
	if self.opts.index:
	   self.indexfile = self.createbdecide_index()
	else:
	   self.indexfile = {}
	   for blastdecide in self.blastdecide: 
	       base = basename(blastdecide)
	       indexfile = str(basename(blastdecide)) + ".index"
	       indexfile2 = str(blastdecide) + ".index"
	       if os.path.exists(indexfile): # or os.path.exists(indexfile_):
		  self.indexfile[base.split(".")[0]] = indexfile
	       else:
		  if os.path.exists(indexfile2):
		     self.indexfile[base.split(".")[0]] = indexfile2
		  else:
		     import warnings
        	     warnings.warn("Index file for blastdecide couldn't be found. Index will be created now.\nThis may take a few minutes depending on the blast output size....\n")
                     self.indexfile = self.createbdecide_index()

	########################## CLUSTERFILE ###############################
	
	if self.opts.g:
	   self.allgeneids=[]
	   for line in file(self.opts.g[0]):
        	self.allgeneids.append(line.strip("\n"))

       ##################### OUTPUT OPTIONS #################################

        uniprottabhead= "AC\tID\tDE\tGN\tTaxonomy\tAccession\tOrganism\tncbi_taxID\tGO_term\tKEGG_id"

        if self.opts.o:
	    self.outputbase = self.opts.o
	else:
	    self.outputbase = self.clusterfile
	
	ph=open(self.outputbase + ".clusterdict.p", "wb")      
	ph2=open(self.outputbase + ".preclusterdict.p", "wb")
	pnh=open(self.outputbase + ".noclusterdict.p", "wb") 
        	
	with open(self.outputbase + ".hitcounts.txt", "w") as cf:
	    cf.write("CL_ID FLAG HITS NOHITS CLUSTERSIZE\n")	    
        with open(self.outputbase + ".cl_nohitslist.txt", "w") as cnh:
	    cnh.write("GENES THAT ARE IN A CLUSTER BUT HAVE NO UNIPROT HITS\n")
        with open(self.outputbase + ".single_nohitslist.txt", "w") as nh:
	    nh.write("GENES THAT ARE NOT IN CLUSTER AND HAVE NO UNIPROT HITS\n")		
	with open(self.outputbase + ".uniprottab", "w") as nh: 
	     nh.write("CLUSTERID\tGENEID\t"+uniprottabhead+"\n")	
	with open(self.outputbase + ".uniprottab.singles.txt", "w") as nh: 
	     nh.write("GENEID\tACTUAL_SIZE\tUNIPROTTAB\n")
	     
	with open(self.outputbase + ".allwithhits.txt", "w") as nh: nh.write("")
	with open(self.outputbase + ".allinfo.txt", "w") as nh: nh.write("")
        with open(self.outputbase + ".single_allinfo.txt", "w") as nh: nh.write("")
	    
	##################### START THREADS ##################################
        self.THREAD_LIMIT = self.opts.t                # This is how many threads we want
        
 
	if self.opts.c:
		self.jobs1 = Queue.Queue()         # This sets up the queue object to use 5 slots
        	self.readincluster()
		self.singlelock = threading.Lock()      # This is a lock so threads don't print trough each other (and other reasons)        
                self.StartQueueAndThreads_1()#self.mythread,self.readincluster())
                stop = time.time()
                timeused = (stop-self.start)
                self.printer("### Threading took: {0} seconds ({1} threads)\n".format(timeused,self.THREAD_LIMIT))	
        else:
                sys.exit("**ERROR!** No input cluster-file registered. Provide one with -c or see -h for the help file\n")

	self.print_files()
	if self.opts.g:
		self.jobs2 = Queue.Queue() 
		self.genesoutcluster()
		self.singlelock = threading.Lock()
		self.StartQueueAndThreads_2()#self.mythread2,self.genesoutcluster())
                stop2 = time.time()
                timeused = (stop2-stop)
                self.printer("### Threading took: {0} seconds ({1} threads)\n".format(timeused,self.THREAD_LIMIT))

	
	pickle.dump(self.cluster, ph)
	pickle.dump(self.preclust, ph2)
	pickle.dump(self.nocluster, pnh)
	ph2.close()
	ph.close()
	pnh.close()
	self.uh.close()
	for handle in self.blasthfiles:
	    self.blasthfiles[handle].close()
        self.timing= (time.time() - self.start) /60
	self.printer("### [parseclinfo.py] Time used: "+str(round(self.timing*60)) + " seconds ("+str(round(self.timing)) + " min)\n")

if __name__ == '__main__':
  
     try:
         obj = GetUniprotInfo()
         obj.opts = obj.parser.parse_args(sys.argv[1:])
	 obj.printer("\n### " + sys.argv[0] + " initialized at " + obj.timestarted )
	 obj.printer("### OPTIONS : " + str(obj.opts))
         obj.mainthing()
     except IOError as i:
	 print "I/O error({0}): {1}".format(i.errno, i.strerror)
	 import traceback
	 traceback.print_exc()

     except Exception,e: 
	 print str(e)
	 import traceback
	 traceback.print_exc()

