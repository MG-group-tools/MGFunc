import threading,Queue
import time
from datetime import datetime as dt
import argparse
import traceback
import random
import subprocess
import shutil
import os
import re
import sys
import glob
from operator import itemgetter
import shutil
#from fetch_goid import GOFetch



helpstr= '''
description: 
This script runs the Goatools gene enrichment analysis on the given cluster file using the Assosiation and Population files. 
'''

epi="Author: Asli I. Ozen (asli@cbs.dtu.dk) and Kosai Iskold(kosai@cbs.dtu.dk)" 


class rungoatools():

        def __init__(self):
            self.start = time.time()
            d_ = dt.today()
	    self.timestarted = d_.strftime("%d-%m-%Y %H:%M:%S")
	    self.parseArgs()
	    self.EnrichedGOs = {}
        def parseArgs(self):
                parser = argparse.ArgumentParser(description=helpstr, epilog=epi)
                parser.add_argument("-i", help="Your cluster-file with GO-terms. (Obligatory).",nargs=1)
                parser.add_argument("-a",help="Association file. A file with gene names and their GO terms in a tab separated format",nargs=1)
		parser.add_argument("-p",help="Population file. If you have an association file, population is generally the first column",nargs=1)
              # parser.add_argument("-u",help="Uniprot file. If -a is not given, you need to provide a Uniprot tab-file.",nargs=1)
	        parser.add_argument("-s", help="Relative location or abs. path for the find_enrichment.py", default="/home/projects8/pr_53035/people/asli/Enzymes/Pipeline/scripts")
                parser.add_argument("-o",help="Name of output file base",nargs=1)
		parser.add_argument("-t",type= int, help="Thread limit",default=2)
		parser.add_argument("-c", "--cleanup", action="store_true" , help="cleanup individual enrichment results")
		parser.add_argument("-v", "--verbose", action="store_true" , help="increase output verbosity")
#               print "\nWorkerBee initialized at " + self.d_.strftime("%H:%M %d/%m-%Y") + "\n"
#               self.args = parser.parse_args()
                self.parser = parser



        def organize(self, sourcelist,destination):
	    for source in sourcelist:
	        os.system ("mv"+ " " + source + " " + destination)

	def StartQueueAndThreads(self,inputlist):
		#self.population = "population.dat"
		#fpop = open(self.population,"w")
		self.cl2file = {}
		myrange =  range(self.THREAD_LIMIT)
		self.ThreadNo = list(myrange[::-1])  ### fx if self.THREAD_LIMIT= 3, self.ThreadNo= [3,2,1]
                for x in myrange:
                        # This is the thread class that we instantiate.
#                       workerbeast().start()
                        threading.Thread(target=self.mythread).start()
		genelist = []
		outstring = ""
		self.output = open(self.out, "w")
		self.outputnice = open(self.outnice, "w")
		with open(inputlist,"r") as fid:
			c = 0
			line = fid.readline().rstrip()
			while True:
				c+=1
				if not line or line== "":
					break
				lineS = line.split()

				stud = str(lineS[0]) +".geneids"
				#print stud
				#if typ != "N" :
				with open(stud,"w") as fout:
					clid = lineS[0]
					self.cl2file[stud] = clid
				     	genelist = lineS[2].split(",")
				     	outstring  = "\n".join(genelist)+"\n"
				     	fout.write(outstring)
				     	#fpop.write(outstring)
			        self.jobs.put(stud,block=True,timeout=1)
				line = fid.readline().rstrip()
				
		self.jobs.join()
		self.output.close()
		self.outputnice.close()
		destination="/".join(os.path.abspath(self.args.o[0]).split("/")[:-1])
		self.organize(glob.glob("*.enrichment"), destination)
		self.organize(glob.glob("*.goinfo"), destination)
		self.organize(glob.glob("*.goinfo.ranked"), destination)
		self.organize(glob.glob("*.goids"), destination)
		self.organize(glob.glob("*.geneids"), destination)

	def mythread(self):
                therandom = self.ThreadNo.pop() + 1
		self.singlelock.acquire()
		self.printer("starting thread number:" + str(therandom))
		self.singlelock.release()
		
		while True:
			try:
				study = self.jobs.get(True,1)
				a= self.run_goatools(study)
				#print a
				if a == 1: 
				   self.jobs.task_done()
				else: 
				   print "There was an error with " + study 
				   self.jobs.task_done()
			except:
				break



        def run_goatools(self,study):

	    try:
	    
	       ############## RUN find enrichment ############################
	       tmp_sh = "tmp.cmd." + study + ".sh"
	       cmd = "python2.7 " + self.args.s + "/goatools_scripts/find_enrichment.py --alpha=0.05 --fdr " + study + " " + self.population + " " + self.association + " > " + study.split(".")[0] + ".enrichment"   
	       ftout = open(tmp_sh,"w")
	       ftout.write(cmd)
	       ftout.close()
	       err=open(study + ".stderr.txt", "w")	
	       dummy = subprocess.call(["sh",tmp_sh], stderr = err)
	       os.remove(tmp_sh)
	       #print "Dummy: ", dummy
	       worked=1
	       if dummy != 0:
	          worked=0	          
		  print "**ERROR!** Unable to run find_enrichment with file: ", study
		  return 0
	       
	       
	       ############## RUN FetchGoid ############################
	       if worked != 0:
	          #print "YAYY WORKED" + study
		  enrgolist_file = self.goaParser(study.split(".")[0] + ".enrichment")	       

		  tmp_sh2 = "tmp.cmd2." + enrgolist_file + ".sh"
		  cmd = "python2.7 " + self.args.s + "/fetch_goid.py -i " + enrgolist_file + " -d > " + enrgolist_file + ".goinfo"	     
		  with open(tmp_sh2,"w") as ftout:
	               ftout.write(cmd)
		  err2=open(enrgolist_file + ".stderr2.txt", "w")
		  dummy = subprocess.call(["sh",tmp_sh2], stderr = err2)
		  os.remove(tmp_sh2)
		  #print "Second Dummy: ", dummy
		  if dummy != 0:	          
		     print "**ERROR!** Unable to run fetch_goid with file: ", enrgolist_file
		     return 0

		  studygolist = []

		  for line in file(enrgolist_file + ".goinfo"):	       
	              l = line.strip("\n").split("\t")
	              goterm = l[2]
		      gores = self.EnrichedGOs[study.split(".")[0]][goterm]
	              studygolist.append([goterm, l[0], l[1],l[3], gores[1], gores[2],gores[3]])

		  studygolist.sort(key=itemgetter(4))
		  feout= open( enrgolist_file +  ".goinfo.ranked", "w")
		  for el in studygolist: 
	             for g in el: 
	        	 feout.write("%s\t" % (g))
		     feout.write("\n")
        	  feout.write("\n")
		  feout.close()


		  self.singlelock.acquire()	    	    
		  out = open(study.split(".")[0] + ".enrichment", "r")		 
		  self.output.write("## CLUSTER: %s\n%s\n" % (study.split(".")[0], out.read()))
		  self.singlelock.release()

		  self.singlelock.acquire()  
		  #print "WRITTEN ENRICHALL", study  	    
		  out2 = open(enrgolist_file + ".goinfo.ranked", "r")
		  self.outputnice.write("## CLUSTER: %s\n%s\n" % (study.split(".")[0], out2.read()))
		  self.singlelock.release()
 

 	       err.close()
	       err2.close() 
	       
	       
	       os.remove(study + ".stderr.txt")
	       os.remove(enrgolist_file + ".stderr2.txt")

	       return 1
	       
	    except IOError as i:					    
	    	    err.write("I/O error({0}): {1}".format(i.errno, i.strerror))
		    self.singlelock.release()
		    traceback.print_exc()
		    return 0
	    except Exception,e:
	    	    err.write(str(e))
		    self.singlelock.release()
		    traceback.print_exc()
		    return 0
	    	    



        def goaParser(self,enrfile):
                fid = self.gzipopen(enrfile)
		
                fout = open(enrfile.split(".")[0] + ".goids","w")
		#ClusterPat = re.compile("CL\_(\d+)")
		EntryEnd = 0
		
		whitecount = 0
		clid = enrfile.split(".")[0]
		while 1:
			fidline = fid.readline()
			fidline = fidline.rstrip()
			if not fidline:
                                whitecount += 1
                                if whitecount == 5:
                                        break
			elif fidline[0:2] == "id":
				whitecount = 0

				EntryEnd = 0
				self.EnrichedGOs[clid] = {}
				header = fid.readline()
				fidline = fid.readline().rstrip() #goterm1 blabla
				whitelinecount = 0
				while EntryEnd == 0:
					if fidline:    
						if fidline[0:2] == "GO":  
							Entries = fidline.split("\t")
							GOterm = Entries[0]    
							enr = Entries[1]
							name=Entries[2]
							p_holm=float(Entries[7])
							p_fdr=Entries[9]
							
							if GOterm in self.clfile[clid]:
							   self.EnrichedGOs[clid][GOterm] = [name, enr, p_holm,p_fdr]

							fidline = fid.readline().rstrip()
							continue
                                	        else:
							break
					if not fidline:
						break
				
		for term in self.EnrichedGOs[clid]: 
		    fout.write(term+"\n")
		   
		fout.close()
		return enrfile.split(".")[0] + ".goids"

        def gzipopen(self,fileID):
                if fileID[-3:] == ".gz":
                        return gzip.open(fileID)
                else:
                        return open(fileID,"rU")
	
        def readincluster(self):
	    self.clfile={}
	    cid = self.gzipopen(self.args.i[0])
	    for line in cid.readlines():
	        linelist= line.split("\t")
		clid=linelist[0]
		goterms=linelist[4].split(",")
		self.clfile[clid]=[]
		for i in goterms:
		    new_i= i.split(";")
		    #print new_i
		    if len(new_i) > 0:
		       for j in new_i:
		           self.clfile[clid].append(j)	    
                #print self.clfile[clid]
        def mainthing(self):
	        if not self.args.o:
		   self.out = self.args.i[0] + ".enrichmentall"
		   self.outnice = self.args.i[0] + ".enrichmentall.ranked"
		else:
		   self.out = self.args.o[0] + ".enrichmentall"
		   self.outnice = self.args.o[0] + ".enrichmentall.ranked"
                if not self.args.a:
                        sys.exit("**ERROR!** You need to provide an association file\n")
		else:	
			self.association = self.args.a[0]
		if not self.args.p:
		       tmp_file_name = "tmp.pop.file.sh"
		       self.population = "population.dat"
		       with open(tmp_file_name,"w") as ftout:
			       ftout.write("cut -f1 " +self.args.a[0]+"> "+self.population)
		       dummy = subprocess.call(["sh",tmp_file_name])
		       os.remove(tmp_file_name)
		       if dummy != 0:
			       sys.exit("**ERROR!** Unable to create population-file. Make sure association file is TAB separated.")
		else:
		   self.population= self.args.p[0]
		
		
		self.readincluster()
		self.THREAD_LIMIT = self.args.t                # This is how many threads we want
                self.jobs = Queue.Queue()         # This sets up the queue object to use 5 slots
                self.singlelock = threading.Lock()      # This is a lock so threads don't print trough each other (and other reasons)
                if self.args.i:
                        self.StartQueueAndThreads(self.args.i[0])
                        stop = time.time()
                        timeused = (stop-self.start)
                        self.printer("### Threading took: {0} seconds ({1} threads)\n".format(timeused,self.THREAD_LIMIT))
                else:
                        sys.exit("**ERROR!** No input cluster-file registered. Provide one with -i or see -h for the help file\n")
		

	       

                
	def printer(self,mystring):
		if self.args.verbose: print mystring


if __name__ == "__main__":
    try:
        myclass = rungoatools()
        myclass.args = myclass.parser.parse_args(sys.argv[1:])
	if myclass.args.verbose:
           myclass.printer("\n### "+sys.argv[0]+" initialized at "+ myclass.timestarted + "\n")
           myclass.printer("### OPTIONS: "+str(myclass.args)+"\n")
        myclass.mainthing()
    except IOError as i:
        print "I/O error({0}): {1}".format(i.errno, i.strerror)
    except Exception,e:
        print str(e)
        import traceback
        traceback.print_exc()


### TO DO

## Print the errors in a file

### HAve to devide association files for molecular function, biol process and cell components
