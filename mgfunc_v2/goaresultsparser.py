from __future__ import division
import argparse
from datetime import datetime as dt
import time
import os
import sys
import gzip
import re

class main:
        def __init__(self):
                self.start = time.time()
                self.d_ = dt.today()
                self.timestarted = self.d_.strftime("%d-%m-%Y %H:%M:%S")
                self.parseArgs()

        def parseArgs(self):###GETTING ARGUMENTS FROM COMMANDLINE###
                parser = argparse.ArgumentParser(prog="goaresult.py",usage="goaresultparser.py -i <goa result table> -o <ouput filename>",epilog="Example: python2.7 goaresultparser.py -i goa.enrichment -o enrichment_analyzed.txt\n\nWritten by Kosai+Asli, JUNE 2014. Last modified JUN 2014.",description="Desctription: Parses a goatools output table and outputs")
                parser.add_argument("-i",metavar="enrichment", help="Goatools enrichment table",nargs=1,required=True)
		parser.add_argument("-c",metavar="cluster", help="Cluster file with GOterms",nargs=1,required=True)
                parser.add_argument("-o",metavar="OUTPUT NAME",help="output-name",nargs=1,required=True)
                parser.add_argument("-v",help="Verbose. Prints out progress and details to stdout output. Write \"-v\" with no arguments in commandline. Default is off.",action="store_true")
#               self.self.args = parser.parse_self.args()
                self.parser = parser

        def goaParser(self):
                fid = self.gzipopen(self.args.i[0])
		
                fout = open(self.args.o[0] + ".results","w")
		ClusterPat = re.compile("CL\_(\d+)")
		EntryEnd = 0
		EnrichedGOs = {}
		i = 0
		whitecount = 0
		while 1:
			fidline = fid.readline()
			fidline = fidline.rstrip()
			if not fidline:
                                whitecount += 1
                                if whitecount == 5:
                                        break
			elif fidline[0] == "#":
				whitecount = 0
				#GoaList.append([])
				EntryEnd = 0
				CL_ID = re.findall(r"CL\_(\d+)",fidline)
				print "We are at !!!:",CL_ID[0]
				EnrichedGOs[CL_ID] = []
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
							p_holm=Entries[7]
							p_fdr=Entries[9]
							if GOterm in self.clfile[CL_ID] and enr == "e":
							   EnrichedGOs[CL_ID].append([GOterm, name, enr, p_holm,p_fdr])

							fidline = fid.readline().rstrip()
							continue
                                	        else:
							break
					if not fidline:
						break
				
				print "##########################"
				print GoaList[i]
				print "###########################"
				i+=1
		fid.close()
		fout.close()


## Clutesr = CL_1
#headers
#goterm1 blabla
#goterm2 blabla
#goterm3 blabla
## Cluster = CL_2
       
#  -h, --help     show this help message and exit
#  -i inputlist   list of GO ids to be fetched)
#  -o outputname  output FILE basename. results will be saved as
#                 'basename.extention'. Path name can be included
#  -g             get info for all the GO ids in the list
#  -p             get a list of parents for all the GO ids in the list
#  -cp common     get the common parents for all the GO ids in the list
#  -d             get all possible distancesto the root
#  -md maxdist    calculate max distance of each GO id in the list to the root
#  -v, --verbose  Increase print out size
	def outputnice(self):
	    from fetch_goid import GoFetch
	    obj = GOFetch()
	
	
	def StartQueueandThread
	
        def readincluster(self):
	    self.clfile={}
	    cid = self.gzipopen(self.args.c[0])
	    for line in cid.readlines():
	        linelist= line.split("\t")
		clid=linelist[0]
		goterms=linelist[4].split(",")
		self.clfile[clid]=[]
		for i in goterms:
		    new_i= i.split(";")
		    for j in new_i:
		        self.clfile[clid] = j
		
	    
<        def gzipopen(self,fileID):
                if fileID[-3:] == ".gz":
                        return gzip.open(fileID)
                else:
                        return open(fileID,"rU")


        def printer(self,string): #surpressing output print if -q (quiet) is on (NEW: shows output if -v is on)
                if self.args.v:
                        print string,

        def mainthing(self):
	        self.readincluster()
                self.goaParser()
                timeused = (time.time() - self.start) / 60
                self.printer("### Time used: "+str(round(timeused/60,1))+" hours)\n")



if __name__ == "__main__":
    try:
        myclass = main()
        myclass.args = myclass.parser.parse_args(sys.argv[1:])
        myclass.printer("\n### "+sys.argv[0]+" initialized at "+ myclass.timestarted + "\n")
        myclass.printer("### OPTIONS: "+str(myclass.args)+"\n")
        myclass.mainthing()
    except IOError as i:
        print "I/O error({0}): {1}".format(i.errno, i.strerror)
    except Exception,e:
        print str(e)
        import traceback
        traceback.print_exc()

