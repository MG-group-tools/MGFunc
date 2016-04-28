#!/usr/bin/env python
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen (asli@cbs.dtu.dk)
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)

import sys
import MySQLdb
import urllib2
from operator import itemgetter 
import time
from datetime import datetime as dt
import argparse
import cPickle as pickle


helpstr= '''
description: 
This script retrieves information about a given GO term list from the mysql database, local or ebi.
'''

epi="Author: Asli I. Ozen (asli@cbs.dtu.dk) and Kosai Iskold(kosai@cbs.dtu.dk)" 
class KEGGFetch:
     def __init__(self):
     	 self.start = time.time()
	 d_ = dt.today()
         self.timestarted = d_.strftime("%d-%m-%Y %H:%M:%S")
	 self.parseArgs()										
         self.db = MySQLdb.Connect(db='KEGG_PATHWAY', user='asli', passwd='Pass9fy', host='mysql.cbs.dtu.dk')
	 self.dbko = MySQLdb.Connect(db='KEGG_PATHWAY_ko', user='asli', passwd='Pass9fy', host='mysql.cbs.dtu.dk')
         self.cursor = self.db.cursor()
	 self.cursorko=self.dbko.cursor()
	 self.pathways={}
	 self.keggs={}
	 self.cluster={}
	 	 
     def parseArgs(self):
         self.parser = argparse.ArgumentParser(description=helpstr, epilog=epi)
         self.parser.add_argument("-c", metavar ="list", type=str, help="ClusterFile",nargs=1,required=True)
	 self.parser.add_argument("-o", metavar ="outputname", type=str, help="output FILE basename. results will be saved as 'basename.extention'. Path name can be included",nargs=1)
	 self.parser.add_argument("-g",  action="store_true", help="get info for all the GO ids in the list")
	 self.parser.add_argument("-p", action="store_true",help="get a list of parents for all the GO ids in the list")
	 self.parser.add_argument("-cp", metavar ="common",help="get the common parents for all the GO ids in the list")
	 self.parser.add_argument("-d",  action="store_true", help="get all possible distancesto the root" )
	 self.parser.add_argument("-md", metavar ="maxdist", type=str, help="calculate max distance of each GO id in the list to the root")
	 self.parser.add_argument("-v", "--verbose", action="store_true" , help="Increase print out size")
	 return 1

     def readclusters(self):        
	 if os.path.exists(self.cldictfile):
	    self.cldict= pickle.load(open(self.cldictfile))
	 else: 	 
	     continue
	     #for each blabla: 
	 
     def getkeggids(self):       
	 for clid in self.cldict: 
	     for keggid in self.cldict[clid][3]:
	         if keggid in self.keggs:
		    print("keggid found")
		    self.cluster[clid][keggid] = self.keggs[keggid]
		 else: 
		    self.keggs[keggid] = self.getkegginfo(keggid)
	     	    self.cluster[clid][keggid] = self.keggs[keggid]


     def getkegginfo(self, keggid): 
     	 getinfo = "SELECT pathname,reaction FROM entry_name where name LIKE '%"+ keggid + "%';"
     	 return self.execute(getinfo)cursor.fetchall()
	 
	 
    sho 
     def getpathinfo(self, pathid): 
         getinfo = "SELECT * FROM pathway where name='"+ keggid + "';"
	 
	 
	     where name='"+ keggid + "';"
	 return self.execute(getinfo)
     
     
     
     
     
     
     
     
     
     
     
     def writeoutput(self):
     
	 
	 
	 
