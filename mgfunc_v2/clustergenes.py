#!/usr/bin/env python
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen (asli@cbs.dtu.dk)
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
import sys, gzip
import re, string
import argparse
import os
from Bio.Blast import NCBIStandalone
from operator import itemgetter, attrgetter
import pickle
import shutil
sys.path.append('/home/projects5/pr_53035/people/asli/bin/lib/python2.7/site-packages')
import networkx as nx
from os.path import basename
import time
from datetime import datetime as dt
import matplotlib.pyplot as plt

usage = sys.argv[0] + " [-h] -p P [-b B] [-a] [-o] [-m M] [-v]\n" 
helpstr = '''
description: This script clusters genes from a given blastdecide file using reciprocal best hits or triple-hits graphs.
'''

epi="Author: Asli I. Ozen (asli@cbs.dtu.dk)" 


class ClusterGenes:
    def __init__(self):
	self.start = time.time()
	d_ = dt.today()
        self.timestarted = d_.strftime("%d-%m-%Y %H:%M:%S")
	self.parseArgs()										
    def printer(self,string): 
        if self.opts.verbose:
           print string
    def parseArgs(self):	
	self.parser = argparse.ArgumentParser(description=helpstr, epilog= epi)
	#parser.add_option('-p', '--path', dest="P", default="./", type=str, help="path to the blast output files")
	#parser.add_argument('-o', '--opath', dest="OP", default="../parsed/", type="str", help="path to the parse output files")    
	#parser.add_argument('-d', '--decide', dest="D", type="str", help="blast tabular output file")

	self.parser.add_argument("-p", metavar ="preblast", type=str, help="output from blastdecide as input FILE (required)",nargs=1,required=True)
	self.parser.add_argument("-o", metavar ="outputname", type=str, help="output FILE basename. results will be saved as 'basename.extention'. Path name can be included",nargs=1)
	self.parser.add_argument("-b", metavar ="INT",type=int,default=3, help="number of besthits to use in reciprocal graph(default=3)",nargs=1)
	self.parser.add_argument("-m", metavar = "[reciprocal|tri]", default="reciprocal", type=str, help="All Reciprocal best-hits or min triplet reciprocal best hits method (default=reciprocal)")
	self.parser.add_argument("-a", "--allhits", action = "store_true", help="allhits will be used in reciprocal graph, not only besthits")
	self.parser.add_argument("-ob", action = "store_true", help="if option given, outputs only besthit for each gene, clusters are not generated")
	self.parser.add_argument("-v", "--verbose", action="store_true" , help="increase output verbosity")
	return 1
    
    def get_path(self, name):      
	return os.path.dirname(os.path.abspath(name)) 
															                  													  
    def read_hits(self): 												  
  															  
	if os.path.exists(self.o + ".bdict.p"):	  # saves time in the second round				   	  
	   														  
	      phl=open(self.o + ".bdict.p","rb")				   	  
	      self.blastDict = pickle.load(phl) 									  
   	      phl.close()												  
	else:														  
           try:														  
              import bldecide as BU										
              obj= BU.BlastDecision()											
              obj.opts = obj.parser.parse_args(['-id', self.blastdecide, '-n', 'nosave', '-v'])
	      self.printer(obj.opts)
	      obj.printer=self.printer
              obj.mainthing()												  
              self.blastDict=obj.blastDict

           except:													
              sys.stderr.write("Can not create blastdict\n")   
	      																								       
	      import traceback
	      traceback.print_exc()	   														
	      sys.exit()	      
	       														
	   try:														
	       ph=open(self.o + ".bdict.p","wb")				
	       pickle.dump(self.blastDict, ph)																									
	       ph.close() 												
	   except:
	       if os.path.exists(self.o + ".bdict.p"):	
	          os.system("rm " + self.o + ".bdict.p")												
	       raise IOError("Pickle doesn't work")
	       sys.exit()

 
    def MakebesthitsDict(self):                # besthits above certain threshold : self.nB

       self.besthits={}
       self.besthit={}
       output=open(self.o + ".besthitlist", "w")
       for i in self.blastDict:           
	   if len(self.blastDict[i]) >= self.nB: 
	      l = sorted(self.blastDict[i], key=itemgetter(3),reverse=True) 
              self.besthit[i] = l[0][0]
              for y in range(0,self.nB):	      	          		  
		  if i in self.besthits:
		      self.besthits[i].append(l[y][0])
		  else:
		      self.besthits[i] = [l[y][0]]
	      
	      output.write("%s\t%s\n" % (i, self.besthit[i]))
	   if len(self.blastDict[i]) < self.nB:
	      j = sorted(self.blastDict[i], key=itemgetter(3),reverse=True)
	      self.besthit[i] =j[0][0]
	      for x in range(0,len(j)):	      	          		  
		  if i in self.besthits:
		      self.besthits[i].append(j[x][0])
		  else:
		      self.besthits[i] = [j[x][0]]

	      output.write("%s\t%s\n" % (i, self.besthit[i]))
	   
       output.close()
       return None
    
    
    def MakeallhitsDict(self):
        self.allhits={}
	for i in self.blastDict:                   
	    l = sorted(self.blastDict[i], key=itemgetter(3),reverse=True)             
	    for y in l:	      	          		  
		  if i in self.allhits:
		      self.allhits[i].append(l[y][0])
		  else:
		      self.allhits[i] = [l[y][0]]
		      
		      
		      
    def MakebesthitDict(self):
       '''
       If Option --onlybesthit is given, this function is called to output
       only the best hits list for each gene in the data.
       '''     
       self.besthit={}
       output=open(self.o + ".besthitlist", "w")
       for i in self.blastDict:                   
	   l = sorted(self.blastDict[i], key=itemgetter(3),reverse=True) 
           self.besthit[i] = l[0][0]
	   output.write("%s\t%s\n" % (i, self.besthit[i]))   
       output.close()
       return None

    
    def Makebesthitsgraph(self, outputname):
        '''
	Prints out Reciprocal Best hits in the data. Data comes from self.besthits, which 
	is previously filtered using number of besthits option.
	''' 									   
        self.besthitsout= open(self.o + ".reciprocalbesthits", "w")
	self.reciprocal= open(outputname ,  "w") 
        #self.Gs=nx.Graph()
	self.Gs= {}
        keylist = self.besthits.keys()
        keylist.sort()
	keepset=set()
	self.edgeconvert = {}
	self.revconvert = {}
	e=1
        for i in keylist:                
	    a=keylist.index(i)
	    #self.Gs.add_node(a,gene=i)
	    for j in self.besthits[i]:    # all the besthits of i (5)
		if j in keylist:  
		   #print j
		   for k in self.besthits[j]:  # all the besthits of j
		       #print k
		       if k == i:    # if one of i's best hits(j) has also a best hit as i
   			  b=keylist.index(j)
			  #self.Gs.add_node(b,gene=j)
	        	  #self.Gs.add_edge(a,b)
			  self.Gs
			  if a in self.edgeconvert and b in self.edgeconvert:
			     e = e-2
			  elif a in self.edgeconvert and b not in self.edgeconvert: 
			     self.edgeconvert[b] = e
			     self.revconvert[e] = b
			     e = e-1
			  elif b in self.edgeconvert and a not in self.edgeconvert:  
			     self.edgeconvert[a] = e
			     self.revconvert[e] = a
			     e = e-1		
			  else: 
			      self.edgeconvert[a] = e
			      self.revconvert[e] = a
			      self.edgeconvert[b] = e + 1
			      self.revconvert[e+1] = b
			      
			  self.reciprocal.write("%d %d\n" % (self.edgeconvert[a], self.edgeconvert[b]))
	                  self.reciprocal.write("%d %d\n" % (self.edgeconvert[b], self.edgeconvert[a]))
			  e+=2
			  keepset.add(a)
	                  keepset.add(b)
			  
	#for i in self.Gs.edges():
           # self.besthitsout.write("%s\t%s\n" % (self.Gs.node[i[0]]["gene"],self.Gs.node[i[1]]["gene"]))	 
	self.besthitsout.close()
	#self.get_cc(self.Gs, self.o + ".besthits.cc")
        #self.get_scc(self.Gs, self.o + ".besthits.scc")
	self.reciprocal.close()
	self.R2sizefile =  outputname + ".size"     
	with open( self.R2sizefile ,  "w") as rsize:
	     rsize.write(str(len(keepset)+1))
	return None
	    
    def Makeallhitsgraph(self): 
        '''
	Prints out Reciprocal Best hits in all the data. Data comes from self.allhits,
	which is not filetered.
	''' 
	
        self.besthitout= open(self.o +".reciprocalallhits", "w")
	
        self.A=nx.Graph()
		      
        keylist = self.allhits.keys()
        keylist.sort()
        for i in keylist:                
	    a=keylist.index(i)
	    self.A.add_node(a,gene=i)
	    for j in self.allhits[i]:    # all the besthits of i (5)
	        if j in keylist:  
		   k=self.allhits[j]  # all the besthits of j
		   if i in k:    # if one of i's best hits(j) has also a best hit as i
   		      b=keylist.index(j)
		      self.A.add_node(b,gene=j)
	              self.A.add_edge(a,b)
	for i in self.A.edges():
            self.besthitout.write("%s\t%s\n" % (self.A.node[i[0]]["gene"],self.A.node[i[1]]["gene"]))	 
	self.besthitout.close()
	
	return None

    def Reciprocaldigraph(self, hitsdict, outputname, singlehit=False):
        '''
	This function produced the directed graphs from the reciprocal best hit information
	First it creates a graph between all genes if they are in each others best hits list. 
	Then it removes the single directed components and only keeps the double directions between 
	two nodes, resulting in the second digraph.
	
	'''
	
	self.reciprocal= open(outputname ,  "w") 	
        k=set()
        keeplist=[]
        self.R=nx.DiGraph()				   # Two digraphs are produced
	#self.R2=nx.Graph()
	keylist = hitsdict.keys()
        keylist.sort()
	keepset=set()
        
	if not singlehit:    
           for i in keylist: 
               a=keylist.index(i)
               self.R.add_node(a,gene=i)
               for j in hitsdict[i]: 
        	   if j in keylist:
        	      b=keylist.index(j)
		      #print a,b
        	      self.R.add_node(b,gene=j)
        	      self.R.add_edge(a,b)
        else:
	   for i in keylist: 
               a=keylist.index(i)
               self.R.add_node(a,gene=i)
               j=hitsdict[i] 
               if j in keylist:
        	  b=keylist.index(j)
        	  self.R.add_node(b,gene=j)
        	  self.R.add_edge(a,b)
       
        Rsize=self.R.size()
	self.printer("### R graph is done, graph size: " + str(Rsize) + "\n")
        ################################### PRINT RECIPROCAL EDGES ##############################
	x=0
	e=1
	edgelist = self.R.edges()
	newedge = [0,0]
	self.edgeconvert = {}
	self.revconvert = {}	
	while len(edgelist) > 0:
	   edge = edgelist.pop()	   
	   revedge = tuple([edge[1],edge[0]])
	   if revedge in edgelist: 
	      x+=1
	      
	      edgelist.remove(revedge)
	      if edge[0] in self.edgeconvert and edge[1] in self.edgeconvert:
	         e = e-2
	      elif edge[0] in self.edgeconvert and edge[1] not in self.edgeconvert: 
	         self.edgeconvert[edge[1]] = e
		 self.revconvert[e] = edge[1]
		 e = e-1
	      elif edge[1] in self.edgeconvert and edge[0] not in self.edgeconvert:  
	         self.edgeconvert[edge[0]] = e
		 self.revconvert[e] = edge[0]
		 e = e-1	      
	      else: 
  	          self.edgeconvert[edge[0]] = e
		  self.revconvert[e] = edge[0]
		  self.edgeconvert[edge[1]] = e + 1
	          self.revconvert[e+1] = edge[1]
	      
	      self.reciprocal.write("%d %d\n" % (self.edgeconvert[edge[0]], self.edgeconvert[edge[1]]))
	      self.reciprocal.write("%d %d\n" % (self.edgeconvert[edge[1]], self.edgeconvert[edge[0]]))
	      
	      
	      #t=frozenset([edge[0],edge[1]])
	      keepset.add(edge[0])
	      keepset.add(edge[1])
              e+=2
	      
	      
	      
	self.R2sizefile =  outputname + ".size"     
	with open( self.R2sizefile ,  "w") as rsize:
	     rsize.write(str(len(keepset)+1))
	self.printer("### Reciprocal connections are finished")
	self.reciprocal.close()       	
	self.printer("### Number of connections:{0} filename: {1}".format(str(x), outputname))
	return 1   
	
	
#	#del R: might save memory

    def SCC(self,connections, output):
        scriptpath =  self.get_path(sys.argv[0])
	path = self.get_path(self.o)
	conpath = self.get_path(connections)
	basecon = os.path.basename(connections)
	conn = conpath + "/" + basecon	
	os.system("ln -s " + scriptpath + "/SCC" + " ./SCC")	
	os.system("ln -s " + conn + " ./SCC.txt")	
	keylist = self.besthits.keys()
        keylist.sort()
	sccout = output + ".sccedges"	
	cmd = "./SCC < " + self.R2sizefile + " > " + sccout	
	os.system(cmd)
	clusters=[]
	for line in file(sccout):
	    if not line[0] == "F":
		   scc = line.split(" ")
		   clusters.append([line,len(scc)-1])
	l = sorted(clusters, key=lambda item: int(item[1]),reverse=True)	       
	with open(output ,  "w") as OUT:
	     a=0	
	     for line in l:
	         #if not line[0] == "F": 
	              a+=1         
		      scc = line[0].split(" ")
		      if int(scc[0]) != 0 :
			 OUT.write("CL_%s\t" % (a))
			 for n in xrange(0, len(scc)-1):
		               realnode = self.revconvert[int(scc[n])]
                	       if n == len(scc)-2 : 
                        	       OUT.write("%s\n" % (keylist[realnode]))
                	       else:
                        	       OUT.write("%s," % (keylist[realnode])) 
        os.system("rm SCC.txt")
	os.system("rm SCC")
        return None
	   
    def MakeReciprocaldigraphs(self):     	
        output=self.o +".bh.Redges"
	outputall=self.o + ".ah.Redges"
	sccoutput= self.o + ".bh.scc"
	sccoutputall= self.o + ".all.scc"
	
	if self.opts.allhits:
	   self.Reciprocaldigraph(self.allhits,outputall)
	   self.SCC(outputall, sccoutputall)
	   #self.get_scc(self.R2, sccoutputall)
	   self.drawgraph(self.R2)
	else:
	   #for i in self.besthits:
	   #    print "{0}\t{1}\n".format(i,self.besthits[i])
	   #self.Reciprocaldigraph(self.besthits,output)
	   self.Makebesthitsgraph(output)   # same thing, using dictionaries instead of graphs. 
	   self.SCC(output, 	 sccoutput)
	   #self.get_scc(self.R2, sccoutput)
	   #self.drawgraph(self.R2)
	   #self.MakeReciprocaldigraph(self.besthits,output2,singlehit=True) 				
        
	self.printer("### Di-graphs and Connected components generated.")							

	
    def get_cc(self,graph,name):
        out=open(name, "w")	
	cc_besthit=nx.connected_components(graph)	
	a=0
        for cc in cc_besthit:
	    a+=1	    
	    out.write("CL_%s\t" % (a))
	    for n in xrange(0, len(cc)):
	        if n == len(cc)-1 : 
			out.write("%s\n" % (graph.node[cc[n]]["gene"]))
	        else:
			out.write("%s," % (graph.node[cc[n]]["gene"]))
	out.close()

    def get_scc(self,graph,name):
        out2=open(name, "w")
	scc_besthit=nx.strongly_connected_components(graph)	
	a=0
        for scc in scc_besthit:
	    a+=1	    
	    out2.write("CL_%s\t" % (a))
	    for n in xrange(0, len(scc)):
	        if n == len(scc)-1 : 
			out2.write("%s\n" % (graph.node[scc[n]]["gene"]))
	        else:
			out2.write("%s," % (graph.node[scc[n]]["gene"]))
	out2.close()


    def TriHits(self):
        k=set()
        self.trigraph= open(self.o + ".trigraph", "w")  
	self.T=nx.Graph()
        self.Tkeep=set()
	self.Trikeep=[]
	if self.opts.allhits:
	   self.Makeallhitsgraph()
	   graph=self.A
	else:
	   self.Makebesthitsgraph()
	   graph=self.Gs
	
	for n in graph.nodes():
	   for i in graph.neighbors(n):
	     for j in graph.neighbors(i):
        	 if j != n:
		    #taxn=graph.node[n]["gene"].split("_")[-1]
		    #taxi=graph.node[i]["gene"].split("_")[-1]
		    #taxj=graph.node[j]["gene"].split("_")[-1]      
        	    if n in graph.neighbors(j):
		       #if taxn != taxi and taxi != taxj and taxn!= taxj:
		          self.T.add_node(n,gene=graph.node[n]["gene"])
			  self.T.add_node(i,gene=graph.node[i]["gene"])
			  self.T.add_node(j,gene=graph.node[j]["gene"])
			  self.T.add_edge(n,i)
			  self.T.add_edge(i,j)
			  self.T.add_edge(j,n)
			  
                	  t=set([n,i,j])
			  self.Tkeep=self.Tkeep|t         	  
			  k.add(frozenset(t))
                	 
	
	for t in self.Tkeep:
	    self.Trikeep.append(graph.node[t]["gene"])               
	#self.Trikeep=list(self.Tkeep)
	self.n_of_tri=len(k)
	self.printer("### Number of Triangles : ",len(k))
        for  tri in k:
	     for z in tri:
	         self.trigraph.write("%s\t" % (self.T.node[z]["gene"])) 
	     self.trigraph.write("\n")
	
	sccoutput= self.o + ".trihits.scc"
	self.get_scc(self.T, sccoutput)	
	self.drawgraph(self.T) 
	self.trigraph.close()   
	                


    def printkeep(self,file,keeplist,ostr): #self, filename, list of keep , output string
        output=open(str(file)+ostr, "w")
        for keys in keeplist:
	   output.write("%s\n" % (keys))
	output.close()
   
    def drawgraph(self,G):
	nx.draw(G)
	plt.savefig(self.opts.o[0] + ".cluster_graph.png")
	nx.draw_random(G)
	plt.savefig(self.opts.o[0] + ".cluster_graph_random.png")
	nx.draw_circular(G)
	plt.savefig(self.opts.o[0] + ".cluster_graph_circular.png")
	nx.draw_spectral(G)
        plt.savefig(self.opts.o[0] + ".cluster_graph_spectral.png")
    
    def mainthing(self):
        self.blastdecide= self.opts.p[0]										     	
        self.nB = self.opts.b[0]												 
	self.v= self.opts.verbose 
	
	if self.opts.o:
 	   self.o = self.opts.o[0]
	   self.opath = self.get_path(self.opts.o[0])
	else:
	   self.o = self.blastdecide + ".b"+ str(self.nB)
	   self.opath = self.get_path(self.blastdecide)
	
	
        if self.opts.ob:  
           self.read_hits()
	   self.MakebesthitDict()
	else: 
	   self.read_hits()
	   	   	     
	   if self.opts.m == "reciprocal":
	      
	      if self.opts.allhits:
	         self.MakeallhitsDict()	         
              else:
	         self.MakebesthitsDict()	         
	      
	      self.MakeReciprocaldigraphs()
	   else: 
	      self.TriHits()
        
	
	
	timeused = (time.time() - self.start) / 60
        self.printer("### Time used: "+str(round(timeused*60)) + " seconds ("+str(round(timeused)) + " min)\n")


if __name__ == '__main__':

     try:
         obj = ClusterGenes()
	 obj.opts = obj.parser.parse_args(sys.argv[1:])
	 if obj.opts.verbose: 
	   obj.printer("\n### " + sys.argv[0] + " initialized at " + obj.timestarted)
	   obj.printer("### OPTIONS : " + str(obj.opts))
         obj.mainthing()
     except IOError as i:
    	 print "I/O error({0}): {1}".format(i.errno, i.strerror)
     except Exception,e: 
	 print str(e)
	 import traceback
	 traceback.print_exc()

	

###############
# INPUT LIST
# pre-made blastdecide file : genecatalogue_vs_uniprot.blastdecide
# 
# OUTPUT LIST
# BlastDict: self.blastdecide.bdict.b5.p
# Reciprocalbesthit: basename(self.blastdecide) + ".besthitlist.txt"
# ReciprocalDigraph: nodes that are connected to each other: basename(self.blastdecide) + ".b"+ str(self.nB)+".Reciprocaldigraph"
# Recprocal-hit based Clusters: basename(self.blastdecide) + ".b"+ str(self.nB)+".Reciprocaldigraph.scc (strongly connected components) i -> j -> i
# Triple-hit based Clusters:basename(self.blastdecide) + ".b"+ str(self.nB)+ ".trihits.cc" i -> j -> k -> i
# 
# OLD functions' output, does the same job
# basename(self.blastdecide)+".b"+ str(self.nB)+ ".reciprocalbesthits", "w")
# basename(self.blastdecide)+ ".b"+ str(self.nB)+".reciprocalbesthit", "w")
#
#
#
# OPTIONS LIST
# ('-p', '--preblast', dest="P", default="", type=str, help="blastdecide output file",nargs=1,required=True)
# ('-b', '--besthits', dest="B", default="3", type=int, help="number of besthits to use in reciprocal graph",nargs=1)
# ('-o', '--onlybesthit', dest="O", action = "store_true", help="outputs only besthit for each gene")
# ('-m', '--method', dest="M", default="reciprocal", type=str, help="reciprocal or triple-reciprocal best-hits method")
# ('-v'...etc.)
#
#
#
#
# FUNCTIONS LIST
# class CLusterGenes:	
# def __init__(self):
# def parseArgs(self):
# def read_hits(self):  													      #
# def MakebesthitsDict(self):
# def MakeallhitsDict(self):
# def MakebesthitDict(self):
#
# def Makebesthitsgraph(self): - undirected graph
# def Makeallhitsgraph(self): - undirected graph
#
# def MakeReciprocaldigraph(self, hitsdict, outputname, singlehit=False): -directed graphs
# def MakeReciprocaldigraphs(self):	  -directed graphs
#    
# def get_cc(self,graph,name):   make cluster from connected components
# def get_scc(self,graph,name):  make clusters from strongly connected components
# def TriHits(self):
# def printkeep(self,file,keeplist,ostr): #self, filename, list of keep , output string
# def mainthing():
