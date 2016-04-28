#!/usr/bin/env python2.7

import sys
import os

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 13:13:45 2013

CLASS-VERSION

@author: Kosai
"""

import cPickle as pickle
from datetime import datetime as dt
import time
import argparse
import gzip

class main:
    '''
    Class version of the cluster2fasta program
    '''
    def __init__(self):
        self.start = time.time()
        self.d_ = dt.today()
	self.timestarted = self.d_.strftime("%d-%m-%Y %H:%M:%S")
	self._D = {}
        self.parseArgs()
    def parseArgs(self):
        parser = argparse.ArgumentParser(prog="cluster2fasta.py", usage="cluster2fasta.py -c mycluster.txt -o mycluster.output -num [-ui uniprot.index\/uniprot.index.p -uf uniprot.fasta] [-ki SAMPLE.index\/SAMPLE.index.p -kf SAMPLE.fasta]", epilog="Written by Kosai+Asli, oct 2013. Last modified apr 2014.")
        parser.add_argument("-ui",metavar="uniprot_index_file",help="Uniprot index file",nargs="*")
        parser.add_argument("-uf",metavar="uniprot_fasta",help="Fasta-file for all uniprot (from swiss2fasta)",nargs="*")
        parser.add_argument("-ki",metavar="sample_index_file",help="Genecatalog index file",nargs=1)
        parser.add_argument("-kf",metavar="sample_fasta",help="Fasta-file for all genecatalog sequences",nargs=1)
        parser.add_argument("-sfi",metavar="sample_list",help="A list of genecatalog index files and fasta files",nargs=1)
        #parser.add_argument("-sfl",metavar="sample_fasta_list",help="Fasta-files list for all genecatalog sequences",nargs=1)
        parser.add_argument("-c",metavar="Cluster-name",help="Cluster-file",nargs=1,required=True)
        parser.add_argument("-o",metavar="Output",help="Output name",nargs=1)
        parser.add_argument("-num",help="Adds 2 coloumns to a new file, with cluster_id\'s, number of sample-genes and number of uniprot ID\'s",action="store_true")
        parser.add_argument("-v",help="Verbose. Prints out progress and details to stdout output. Write \"-v\" with no arguments in commandline. Default is off.",action="store_true")
        #return parser.parse_args("-o testcluster.argsO_tester".split()), parser #testing on windows
        #return parser.parse_args("".split()), parser #testing on window
#        self.args = parser.parse_args()
        self.parser = parser
    
    def fileIn(self,infile):
        if infile[-3:] == ".gz":
            return gzip.open(infile,"r")
        else:
            return open(infile,"r")
    
    def fileOut(self,outfile):
        return open(outfile, "w")
    
    def fileClose(self,cfile):
        cfile.close()
    
    '''
    def dictMaker(i,D_ID): #Create dictionary from index-text file
        D = {}
        if i[0].split(".")[-1] == "index":
            indexline = ""
            for el in D_ID:
                indexline = el.rstrip().split("\t")
                D[indexline[0]] = [indexline[1],indexline[2]]
            self.printer("\n\nDICTIONARY DONE!!!!\n\n")
            return D
        else:
            return pickle.load(D_ID)
    '''
    
    def dictMaker(self,i,D_ID, j): #Create dictionary from index-text file
        	
        if i.split(".")[-1] == "indexed":
            indexline = ""
            for el in D_ID:
                indexline = el.rstrip().split("\t")
                self._D[indexline[0]] = [indexline[1],indexline[2], j]		
            self.printer("\nDictionary done, time used (so far): "+str(round((time.time() - self.start) / 60,3))+" min\n")
            return 1
#        else:
#	    print "Check index file names. :" + i
#            self._D = pickle.load(D_ID)
#            self.printer("\nDictionary done, time used (so far): "+str(round((time.time() - self.start) / 60,3))+" min\n")    
#            return 1
    
    def missingGeneLog(self,genecat,uniprot):
        log = self.fileOut(self.args.o[0] + ".missingGenes.log")
        for el in genecat:
            log.write(el[0]+"\t"+el[1]+"\n")
        for el in uniprot:
            log.write(el[0]+"\t"+el[1]+"\n")
        self.fileClose(log)
    
    def seqExtracter3(self,ID,myD,uni): #Dictionary look-up, one big dictionary
        
	if ID in myD:
            start = int(myD[ID][0])
            stop = int(myD[ID][1])
	    if uni == 1: 
               self.uniprotFasta.seek(start)
               seq = self.uniprotFasta.read(stop-start)
               seq = "".join(seq.split("\n"))
	       return seq,1
	    else:	    
	       fasta = self.fileIn(self._F[int(myD[ID][2])][1])
               fasta.seek(start)
               seq = fasta.read(stop-start)
               seq = "".join(seq.split("\n"))
	       self.fileClose(fasta)
               return seq,1
        else:
            return "",0
    def seqExtracter(self,ID,myD,fasta,uni): #Dictionary look-up, one big dictionary
        if ID in myD:
            start = int(myD[ID][0])
            stop = int(myD[ID][1])
            fasta.seek(start)
            seq = fasta.read(stop-start)
            seq = "".join(seq.split("\n"))
            return seq,1
        else:
            return "",0
    
    def seqExtracter2(self,ID,myD,fasta): #Dictionary look-up, each key is first gene letter
        start = int(myD[ID[0]][ID][0])
        stop = int(myD[ID[0]][ID][1])
        fasta.seek(start)
        seq = fasta.read(stop-start)
        seq = "".join(seq.split("\n"))
        return seq
    
    
    def genecat_list(self):
        clusterID =self.fileIn(self.args.c[0])
	output = self.fileOut(self.args.o[0]+".genecatalog.fasta")
	self._F = {}
	infiles=0
	for line in file(self.args.sfi[0]):
	     index = line.split("\t")[0]
	     fasta = line.split("\t")[1].strip("\n")
	     self._F[infiles] = [index,fasta]
             genecatID = self.fileIn(index) 
             a = self.dictMaker(index,genecatID,infiles) #takes time
             if a ==1 : self.printer("DictMaker worked for " + index)
	     else: self.printer("DictMaker did not work, check index files " + index)	     
	     self.fileClose(genecatID)
             infiles+=1
        suc = 0
        missing = []
        seq = ""
        for line in clusterID:
            L = line.rstrip().split("\t")
            C = str(L[0]) #clusterID
            L2 = L[2].split(",")
            for el in L2:		       
               seq,suc = self.seqExtracter3(el,self._D,0)
               
	       if suc == 1:
            	   output.write(">"+C+":"+el+"\n"+seq+"\n")
               else:
            	   missing.append([el,C])
        #print self._D
	self._D = {}
        self.fileClose(output)             
        self.fileClose(clusterID)
        return missing

    def genecat(self,args,parser):
        clusterID =self.fileIn(args.c[0])
        genecatID = self.fileIn(args.ki[0])
        genecatFasta = self.fileIn(args.kf[0])
        output = self.fileOut(args.o[0]+".genecatalog.fasta")
        a = self.dictMaker(args.ki[0],genecatID,0) #takes time
        if a ==1 : self.printer("DictMaker worked for " + args.ki[0])
	else: self.printer("DictMaker did not work, check index files " + args.ki[0])
        self.fileClose(genecatID)
        GenecatalogD = {}
        cGenecatalog = 1
        suc = 0
        missing = []
        seq = ""
        for line in clusterID:
            L = line.rstrip().split("\t")
            C = str(L[0]) #clusterID
            L2 = L[2].split(",")
            for el in L2:
                seq,suc = self.seqExtracter(el,self._D,genecatFasta,0)
                if suc == 1:
                    if el not in GenecatalogD:
                        GenecatalogD[el] = el[0]+str(cGenecatalog)
                        cGenecatalog += 1
                    #output.write(">"+C+"_"+GenecatalogD[el]+"\n"+seq+"\n")
                    output.write(">"+C+":"+el+"\n"+seq+"\n")
                else:
                    missing.append([el,C])
        #print self._D
	self._D = {}
    #    GenecatalogIDconversion(GenecatalogD)
        self.fileClose(output)
        self.fileClose(genecatFasta)
        self.fileClose(clusterID)
        return missing
 
 
    
    def uniprot(self,args,parser):
        clusterID = self.fileIn(args.c[0])
        uniprotID = self.fileIn(args.ui[0])
        self.uniprotFasta = self.fileIn(args.uf[0])
        ctotfile = os.popen("wc -l "+args.c[0])
        ctot = ctotfile.read()
        ctotfile.close()
        ctot = int(ctot.split(" ")[0])
        rangelist = range(0,ctot,1)
        output = self.fileOut(args.o[0]+".uniprotids.fasta")
        D = self.dictMaker(args.ui[0],uniprotID,0) #takes time
        if D ==1 : self.printer("DictMaker worked for " + args.ui[0])
	else: self.printer("DictMaker did not work, check index files " + args.ui[0])

        self.fileClose(uniprotID)
        seq = ""
        missing = []
        suc = 1
        c = 0
        for line in clusterID:
            c+=1
            L = line.rstrip().split("\t")
            C = str(L[0]) #clusterID
            if L[1] == "N":
	    	continue
            L2 = L[3].split(",")
            for el in L2:
                el = el.split("|")[2]
                seq,suc = self.seqExtracter3(el,self._D,1)
                if suc == 1:
                    output.write(">"+C+":"+el+"\n"+seq+"\n")
                else:
                    missing.append([el,C])
            #if c in rangelist:
                #self.printer("FINISHED "+str(c)+" ENTRIES out of "+str(ctot))
        del D
        self.fileClose(output)
        self.fileClose(self.uniprotFasta)
        self.fileClose(clusterID)
        return missing
    
    def GenecatalogIDconversion(self,D):
        self.printer("\nPrinting GeneConversionTable....")
        fout = self.fileOut("GeneConversionTable.txt")
        for key in D:
            fout.write(key+"\t"+D[key]+"\n")
        fout.close()
        self.printer("DONE!\n")
    
    
    def numberCounter(self,args,parser):
        clusterID = self.fileIn(args.c[0])
        if self.args.o:
            output = self.fileOut(args.o[0]+".genenumbers")
        else:
            output = self.fileOut(args.c[0]+".genenumbers")
        t = "\t"
        n = "\n"
        for line in clusterID:
            L = line.split("\t")
            output.write(L[0]+t+str(len(L[1].split(",")))+t+str(len(set(L[2].split(","))))+n)
        self.fileClose(clusterID)
        self.fileClose(output)

    def printer(self,string): #surpressing output print if -q (quiet) is on
#        if not self.args.quiet:
        if self.args.v:
            print string,
    def read_columns(self,  i, csv_file):
        item=""
        with open(csv_file, 'r') as csvfile:
             for line in csvfile.readlines():
                 array = line.strip("\n").split('\t')
                 item = item + "\n" + array[i]
        return item
    
    def mainthing(self):
#        self.printer("\n***cluster2fasta.py initialized at "\
#        + self.d_.strftime("%H:%M %d/%m-%Y") + "***\n")
#        self.printer("Arguments:\n")
#        self.parseArgs()
        no = 1
        missing1 = []
        missing2 = []
        if bool(self.args.ki)^bool(self.args.kf):
            self.printer("***ERROR!*** Only one of -ki and -kf was provided!\n")
        elif bool(self.args.ui)^bool(self.args.uf):
            self.printer("***ERROR!*** Only one of -ui and -uf was provided!\n")
        elif not self.args.c:
            self.printer("***ERROR!*** No cluster-files(s) provided!\n")
        elif (self.args.ki or self.args.ui) and not self.args.o:
            self.printer("***ERROR!*** No output-name provided!\n")
        else:
            if self.args.ki and self.args.kf and self.args.c and self.args.o:
                self.printer("\tCluster-file: "+self.args.c[0] +"\n\tGenecatalog-index file: "+self.args.ki[0]+"\n\tGenecatalog fasta-file: "+self.args.kf[0]+"\n\tOutput file-name: "+self.args.o[0]+".genecatgenes.fasta\n")
                no = 0
                missing1 = self.genecat(self.args,self.parser)
                self.printer("\nGenecatalog Genes Done! Time (so far): "+str(round((time.time() - self.start) / 60,3))+" min\n")
            if self.args.sfi and self.args.c and self.args.o:
                self.printer("\tCluster-file: \n\t\t"+self.args.c[0] +"\n\tGenecatalog-index files: \n\t\t"+self.read_columns(0, self.args.sfi[0])+"\n\tGenecatalog fasta-files: \n\t\t"+self.read_columns(1, self.args.sfi[0])+"\n\tOutput file-name: \n\t\t"+ self.args.o[0]+".genecatgenes.fasta.gz\n")
                no = 0
                missing1 = self.genecat_list()
                self.printer("\nGenecatalog Genes Done! Time (so far): "+str(round((time.time() - self.start) / 60,3))+" min\n")

            if self.args.ui and self.args.uf and self.args.c and self.args.o:
                self.printer("\tCluster-file: "+self.args.c[0] +"\n\tUniprot-index file: "+self.args.ui[0]+"\n\tUniprot fasta-file: "+self.args.uf[0]+"\n\tOutput file-name: "+self.args.o[0]+".uniprotids.fasta\n")
                no = 0
                missing2 = self.uniprot(self.args,self.parser)
                self.printer("\nUniprot ID\'s Done! Time (so far): "+str(round((time.time() - self.start) / 60,3))+" min\n")
            if self.args.num and self.args.c:
                if not self.args.o:
                    self.printer("\tCluster-file: "+self.args.c[0] +"\n\tOutput file-name: "+self.args.c[0][:-4]+".genenumbers\n")
                else:
                    self.printer("\tCluster-file: "+self.args.c[0] +"\n\tOutput file-name: "+self.args.o[0]+".genenumbers\n")
                no = 0
                self.numberCounter(self.args,self.parser)
                self.printer("\nNumber Calculations Done! Time (so far): "+str(round((time.time() - self.start) / 60,3))+" min\n")
            if no == 1:
                self.printer("none!\n")
            self.missingGeneLog(missing1,missing2)
            timeused = (time.time() - self.start) / 60
            self.printer("Time used: "+str(round(timeused*60))\
            + " seconds ("+str(round(timeused)) + " min)\n")

    def test(self,num):
        self.printer("test")

'''
if __name__ == "__main__":
    myclass = main
    myclass.mainthing
    myclass.test(2)
    self.printer("yoyoooyooyoo")
'''    
    
if __name__ == "__main__":
    try:
	myclass = main()
	myclass.args = myclass.parser.parse_args(sys.argv[1:])
	myclass.printer("\n### "+sys.argv[0]+" initialized at "+ myclass.timestarted + "\n")
        myclass.printer("### OPTIONS: "+str(myclass.args)+"\n")
	myclass.mainthing()
    #except IOError as i:
    #    print "I/O error({0}): {1}".format(i.errno, i.strerror)
    except Exception,e:
        print str(e)
        import traceback
        traceback.print_exc()



##############################
'''
INPUT:
The User inputs an index-file and a fasta-file.
The index file indexes each entry in the fasta file. In the case of -ui and -uf,
-ui would a pickle-file which contains the start and end for the sequences in each
entry of the uniprot file (-uf).

if -num is toggled, the script will not create a fasta-output, but instead
show the number of genecat-genes (sample-genes) and uniprot ID's in each cluster.

OUTPUT:
The output is a fasta file containing the sequences of each uniprot/genecat-gene in the input
from the clusters.


OPTIONS LIST:
"-ui"	"uniprot_index_file": Uniprot index file containing 
"-uf"	"uniprot_fasta": Fasta-file for all uniprot (from swiss2fasta)
"-ki"	"sample_index_file": Sample index file
"-kf"	"sample_fasta": Fasta-file for all sample sequences
"-c"	"Cluster-name": Cluster-file
"-o"	"Output fasta file": Output name
"-num": Adds 2 coloumns to a new file, with cluster_id's, number of sample-genes and number of uniprot ID's
'''
