import threading,Queue
import time
from datetime import datetime as dt
import argparse
import random
import subprocess
import shutil
import os
import re
import sys
import glob



class mygoawrapper():
        start = time.time()
        d_ = dt.today()
	timestarted = d_.strftime("%d-%m-%Y %H:%M:%S")
        cleanup = 0
        def __init__(self):
                parser = argparse.ArgumentParser()
                parser.add_argument("-i", help="Your cluster-file with GO-terms. (Obligatory).",nargs=1)
                parser.add_argument("-a",help="Association file. If you have an association file, you do not need to provide the uniprot-file -U",nargs=1)
		parser.add_argument("-p",help="Population file. If you have an association file, population is generally the first column",nargs=1)
               # parser.add_argument("-u",help="Uniprot file. If -a is not given, you need to provide a Uniprot tab-file.",nargs=1)
	        parser.add_argument("-s", help="Relative location or abs. path for the find_enrichment.py", default="/home/projects8/pr_53035/people/asli/Enzymes/Pipeline/scripts/goatools_scripts")
                parser.add_argument("-o",help="Name of output file",nargs=1)
		parser.add_argument("-t",help="Thread limit",default=3)
#               print "\nWorkerBee initialized at " + self.d_.strftime("%H:%M %d/%m-%Y") + "\n"
#               self.args = parser.parse_args()
                self.parser = parser




	def StartQueueAndThreads(self,inputlist):
		#self.population = "population.dat"
		#fpop = open(self.population,"w")
		myrange =  range(self.THREAD_LIMIT)
		self.ThreadNo = list(myrange[::-1])  ### fx if self.THREAD_LIMIT= 3, self.ThreadNo= [3,2,1]
                for x in myrange:
                        # This is the thread class that we instantiate.
#                       workerbeast().start()
                        threading.Thread(target=self.mythread).start()
		genelist = []
		outstring = ""
		self.output = open(self.out, "w")
		with open(inputlist,"r") as fid:
			c = 0
			line = fid.readline().rstrip()
			while True:
				c+=1
				if not line or line== "":
					break
				lineS = line.split()
				tmpclust = "myTmpCluster."+str(lineS[0])
				
				with open(tmpclust+".tmpCL","w") as fout:
					genelist = lineS[2].split(",")
					outstring  = "\n".join(genelist)+"\n"
					fout.write(outstring)
					#fpop.write(outstring)
			        self.jobs.put(tmpclust+".tmpCL",block=True,timeout=1)
				line = fid.readline().rstrip()
				
		#fpop.close()
		
		self.jobs.join()
		self.output.close()

	def mythread(self):
                therandom = self.ThreadNo.pop() + 1
		self.singlelock.acquire()
		print "starting thread number:",therandom
		self.singlelock.release()
		
		while True:
			try:
				study = self.jobs.get(True,1)
				#print study,self.population,self.association
				#######				
				a= self.run_goatools(study)
				if a == 1: 
				   self.jobs.task_done()
				else: 
				   break
				
				
				
			except:
				break



        def run_goatools(self,study):
	    tmp_sh = "tmp.cmd." + study + ".sh"
	    cmd = "python2.7 " + self.args.s + "/find_enrichment.py --alpha=0.05 --fdr " + study + " " + self.population + " " + self.association	    
	    with open(tmp_sh,"w") as ftout:
	         ftout.write(cmd + "> " + study + ".out"  )
	    err=open(study + ".stderr.txt", "w")
	    dummy = subprocess.call(["sh",tmp_sh], stderr = err)
	    os.remove(tmp_sh)
	    if dummy != 0:
	       sys.exit("**ERROR!** Unable to run find_enrichment")
	    self.singlelock.acquire()
	    out = open(study + ".out", "r")	    
	    self.output.write("%s\n##\n" % out.read())	    
	    self.singlelock.release()
	    err.close()
	    os.remove(study + ".out")
	    os.remove(study + ".stderr.txt")
            return 1
	    
	    
        def check_args_and_run(self):
	        if not self.args.o:
		   self.out = self.args.i[0] + ".enrichment"
		else:
		   self.out = self.args.o[0] + ".enrichment"
                if not self.args.a:
                        sys.exit("**ERROR!** You need to provide either a uniprot tab-file or an association file\n")
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
		self.THREAD_LIMIT = self.args.t                # This is how many threads we want
                self.jobs = Queue.Queue()         # This sets up the queue object to use 5 slots
                self.singlelock = threading.Lock()      # This is a lock so threads don't print trough each other (and other reasons)
                if self.args.i:
                        self.StartQueueAndThreads(self.args.i[0])
                        stop = time.time()
                        timeused = (stop-self.start)
                        print "### Threading took: {0} seconds ({1} threads)\n".format(timeused,3),
                else:
                        sys.exit("**ERROR!** No input cluster-file registered. Provide one with -i or see -h for the help file\n")


	def printer(self,mystring):
		print mystring


if __name__ == "__main__":
    try:
        myclass = mygoawrapper()
        myclass.args = myclass.parser.parse_args(sys.argv[1:])
        myclass.printer("\n### "+sys.argv[0]+" initialized at "+ myclass.timestarted + "\n")
        myclass.printer("### OPTIONS: "+str(myclass.args)+"\n")
        myclass.check_args_and_run()
    except IOError as i:
        print "I/O error({0}): {1}".format(i.errno, i.strerror)
    except Exception,e:
        print str(e)
        import traceback
        traceback.print_exc()


