#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:57:08 2013

CLASS VERSION

Version with -N

@author: Kosai
"""

from __future__ import division
import argparse
from datetime import datetime as dt
import time
import sys
import os
import shutil

class main:
    def __init__(self):
        self.start = time.time()
        self.d_ = dt.today()
        self.timestarted = self.d_.strftime("%d-%m-%Y %H:%M:%S")
	self.parseArgs()

    def parseArgs(self):
        parser = argparse.ArgumentParser(prog="clusterFilter.py",description="Takes a cluster, and selects the clusters that have a specific user-specified number of koala/samples in them.",usage="clusterFilter.py -i cluster.txt -n 4 (-k k3mhn/K3MHN)",epilog="Written by Kosai+Asli, oct 2013. Last modified mar 2014.")
        parser.add_argument("-i",metavar="Cluster-name(s)",help="Cluster-file(s)",nargs="*",required=True)
        parser.add_argument("-n",metavar="sample number",help="Number of different samples's in each cluster",nargs="*")
        parser.add_argument("-k",metavar="sample-names",help="Sample ID's to be selected. If your samples are called S1_x.gene and S2_x.gene, example: '-k HKM'",nargs="*")
	parser.add_argument("-m",metavar="gene-number",help="Number of minimum genes per cluster",nargs="*")
        parser.add_argument("-o",metavar="Output",help="Output name",nargs="*")
	parser.add_argument("-p",metavar="operator",help="Logical operator, AND (&,AND,A) or OR (|,OR,O).")
        parser.add_argument("-v",help="Verbose. Prints out progress and details to stdout output. Write \"-v\" with no arguments in commandline. Default is off.",action="store_true")
        #return parser.parse_args("-i testclusters.txt testclusters2.txt -k hnk3 -o testcluster.argsO_tester".split()), parser #testing on windows
#        self.args = parser.parse_args()
        self.parser = parser
        
    def fileIn(self,infile):
        infile = open(infile,"r")
        return infile
        
    def fileOut(self,outfile):
        outfile = open(outfile, "w")
        return outfile
    
    def fileClose(self,cfile):
        cfile.close()


    def geneParser(self,tmp=1):
        tab = "\t"
        c=0
        outputname = -1
        for cluster in self.args.i:
            outputname += 1
            clusterI = self.fileIn(cluster)
            if tmp == 0:
                clusterO = self.fileOut(self.ori+".filtered")
            elif not self.args.o:
                clusterO = self.fileOut(cluster+".filtered")
		self.tmpFO = cluster+".filtered"
            elif len(self.args.o) == len(self.args.i):
                clusterO = self.fileOut(self.args.o[outputname])
		self.tmpFO = self.args.o[outputname]
            else:
                clusterO = self.fileOut(self.args.o[0] + str(outputname))
		self.tmpFO = self.args.o[0] + str(outputname)
            for line in clusterI:
                if len(line.rstrip()) == 0:
                    continue
 #               genes = line.rstrip().split(tab)[1].split(",") #gene list ie [x1_something,x2_something,x3_something]
		lineS = line.rstrip().split(tab)
		genes = lineS[1].split(",")
		numgenes = len(genes)
                if numgenes >= int(self.args.m[0]):
                    clusterO.write(lineS[0]+tab+lineS[1]+tab+str(numgenes)+"\n")
                    c+=1
            self.fileClose(clusterI)
        self.fileClose(clusterO)
        return c

    
    def nParser(self,tmp=1):
        tab = "\t"
        c=0
        outputname = -1
        for cluster in self.args.i:
            outputname += 1
            clusterI = self.fileIn(cluster)
            if tmp == 0:
                clusterO = self.fileOut(self.ori+".filtered")
            elif not self.args.o:
                clusterO = self.fileOut(cluster+".filtered")
                self.tmpFO = cluster+".filtered"
            elif len(self.args.o) == len(self.args.i):
                clusterO = self.fileOut(self.args.o[outputname])
                self.tmpFO = self.args.o[outputname]
            else:
                clusterO = self.fileOut(self.args.o[0] + str(outputname))
                self.tmpFO = self.args.o[0] + str(outputname)
            for line in clusterI:
                if len(line.rstrip()) == 0:
                    continue
		lineS = line.rstrip().split(tab)
		genes = lineS[1].split(",")
#                genes = line.rstrip().split(",")
                L=[]
                for el in genes:
                    L.append(el.split("_")[0])
                Lsize = len(set(L))
                if Lsize >= int(self.args.n[0]):
                    c += 1
                    clusterO.write(lineS[0]+tab+lineS[1]+tab+",".join(L)+tab+str(Lsize)+"\n")
            self.fileClose(clusterI)
        self.fileClose(clusterO)
        return c
    
    def findingK(self):
        L = []
        goodchars = ["H","M","N","K","3"]
        for el in self.args.k[0]:
            if el.upper() in goodchars:
                L.append(el.upper())
        if "K" in L and "3" in L:
            L.remove("K")
            L.remove("3")
            L.append("K3")
        elif "K" in L and "3" not in L:
            L.remove("K")
            L.append("K3")
        elif "K" not in L and "3" in L:
            L.remove("3")
            L.append("K3") #meaning it's ok to just write -k 3 (= K3)
        return set(L)

       

    def kParser2(self,tmp=1):
        tab = "\t"
        c=0
        outputname = -1
        for cluster in self.args.i:
            outputname += 1
            clusterI = self.fileIn(cluster)
	    if tmp == 0:
		clusterO = self.fileOut(self.ori+".filtered")
            elif not self.args.o:
                clusterO = self.fileOut(cluster+".filtered")
            elif len(self.args.o) == len(self.args.i):
                clusterO = self.fileOut(self.args.o[outputname])
            else:
                clusterO = self.fileOut(self.args.o[0] + str(outputname))
            for line in clusterI:
		if len(line.rstrip()) == 0:
		    continue
		lineS = line.rstrip().split(tab)
		genes = lineS[1].split(",") #gene list ie [x1_something,x2_something,x3_something]
                Sset=set([])
                for el in genes:
		    Sset.update([el.split("_")[0]])
		Lsize = len(Sset)
		if self.operator == "|":
		    if Sset.intersection(self.kset):
			c+=1
			clusterO.write(lineS[0]+tab+lineS[1]+tab+",".join(list(Sset))+tab+str(Lsize)+"\n")
		else:
		    if self.kset.issubset(Sset):
			c+=1
			clusterO.write(lineS[0]+tab+lineS[1]+tab+",".join(list(Sset))+tab+str(Lsize)+"\n")
            self.fileClose(clusterI)
        self.fileClose(clusterO)
        return c

    def determineOperator(self):
	if self.args.p:
            if self.args.p in ["|","OR","or","O"]:
                self.operator = "|"
            elif self.args.p in ["&","AND","and","A"]:
                self.operator = "&"
            else:
                self.printer("Invalid operator, using AND (default)\n")
        else:
            self.operator = "&"
	
    def makeKset(self):
	k = self.args.k
	self.kset = set([])
	for el in k:
	    k2 = el.split(",")
	    for el2 in k2:
		self.kset.update([el2])

    def printer(self,string): #surpressing output print if -q (quiet) is on
#        if not self.args.quiet:
        if self.args.v:
            print string,

    def mainthing(self):
#        self.printer("\nclusterParser initialized at " + self.d_.strftime("%H:%M %d/%m-%Y") + "\n")
	if not self.args.o:
	    self.ori = self.args.i[0]
	else:
	    self.ori = self.args.o[0]
        if self.args.m and not self.args.k and not self.args.n:
	    self.hits = self.geneParser(0)
	elif self.args.m and (self.args.k or self.args.n):
	    buffer1 = self.geneParser()
	    self.args.i = [self.tmpFO]
	    deleteme = self.tmpFO
	if self.args.n and not self.args.k:
	    self.hits = self.nParser(0)
	elif self.args.k and not self.args.n:
	    self.determineOperator()
       	    self.makeKset()
            self.hits = self.kParser2(0)
	elif self.args.k and self.args.n:
	    shutil.copy(self.args.i[0],self.args.i[0]+".tmp")
	    self.args.i = [self.args.i[0]+".tmp"]
	    deleteme2 = self.args.i[0]
	    buffer1 = self.nParser()
	    deleteme = self.tmpFO
	    self.args.i = [self.tmpFO]
	    self.determineOperator()
	    self.makeKset()
	    self.hits = self.kParser2(0)
	elif not self.args.m:
	    sys.exit("***ERROR!***: You didn't specify any filters, please type \"python2.7 clusterParser.py -h\" for help.\n")
	try:
	    if os.path.exists(deleteme):
		os.remove(deleteme)
	    if os.path.exists(deleteme2):
		os.remove(deleteme2)
	except Exception:
	    sys.exc_clear()
        self.timeused = (time.time() - self.start) / 60
        self.printer("Total "+str(self.hits)+" cluster(s) found!\n")
        self.printer("### Time used: "+str(round(self.timeused*60)) + " seconds ("+str(round(self.timeused)) + " min)\n")

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

#    myclass = main()
#    myclass.args = myclass.parser.parse_args(sys.argv[1:])
#    myclass.mainthing()


######################
'''
INPUT:
The User inputs a cluster file and some filters, for example the minimum number
of samples/Koalas in each cluster. The script will then extract all clusters
with this or a higher number of samples and put them in a new cluster-file.
You can also specify specific sample names, in the case of Koala's, H,K,M or B.

OUTPUT:
The output is a subset of the original cluster containing only the specified
samples (koala id's)


OPTIONS LIST:
-i: Cluster-file(s) (works with multiple input clusters)
-n: Number of different Koala's in each cluster
-k: Specific Koala ID's (sample) to be selected. Use [H,K,M,N], example: '-k HKM'
-o: Output name
'''


