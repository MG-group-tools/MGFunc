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


helpstr= '''
description: 
This script retrieves information about a given GO term list from the mysql database, local or ebi.
'''

epi="Author: Asli I. Ozen (asli@cbs.dtu.dk) and Kosai Iskold(kosai@cbs.dtu.dk)" 
class GOFetch:
     def __init__(self):
     	 self.start = time.time()
	 d_ = dt.today()
         self.timestarted = d_.strftime("%d-%m-%Y %H:%M:%S")
	 self.parseArgs()										

         self.db = MySQLdb.connect("mysql.ebi.ac.uk", "go_select", "amigo", "go_latest", port=4085)   
	# self.db = MySQLdb.Connect(db='go_termdb', user='asli', passwd='Pass9fy109', host='mysql.cbs.dtu.dk')
	 self.cursor = self.db.cursor()
	 self.parents={}
         self.parents["count"] = {}
         self.parents["info"] = {}
	 

     def parseArgs(self):
         self.parser = argparse.ArgumentParser(description=helpstr, epilog=epi)
         self.parser.add_argument("-i", metavar ="inputlist", type=str, help="list of GO ids to be fetched)",nargs=1,required=True)
	 self.parser.add_argument("-o", metavar ="outputname", type=str, help="output FILE basename. results will be saved as 'basename.extention'. Path name can be included",nargs=1)
	 self.parser.add_argument("-g",  action="store_true", help="get info for all the GO ids in the list")
	 self.parser.add_argument("-p", action="store_true",help="get a list of parents for all the GO ids in the list")
 # self.parser.add_argument("-cp", metavar ="common",help="get the common parents for all the GO ids in the list")
	 self.parser.add_argument("-d",  action="store_true", help="get all possible distancesto the root" )
	 self.parser.add_argument("-m", action="store_true", help="calculate max distance of each GO id in the list to the root")
	 self.parser.add_argument("-v", "--verbose", action="store_true" , help="Increase print out size")
	 return 1

     
     
     def getNODEinfo(self, goid): 
         getinfo = "SELECT * FROM term where acc='"+ goid + "';"
         return self.execute(getinfo)
    
     def fetchparents(self, goid):
         #PARENTS
	 getparents = "SELECT DISTINCT ancestor.*, graph_path.distance, graph_path.term1_id AS ancestor_id \
		       FROM term \
		       INNER JOIN graph_path ON (term.id=graph_path.term2_id) \
		       INNER JOIN term AS ancestor ON (ancestor.id=graph_path.term1_id)\
		       WHERE term.acc='"+ goid + "';"		      
	 return self.execute(getparents) 								        

     def getdist2root(self, goid):
     	 dist2root = "SELECT DISTINCT term.*,p.distance FROM\
		       term\
		       INNER JOIN graph_path AS p ON (p.term2_id=term.id)\
		       INNER JOIN term AS root ON (p.term1_id=root.id)\
		       WHERE root.is_root=1 AND term.acc='"+ goid + "';"
         return self.execute(dist2root) 
     def getmaxdist(self,goid):
         maxdist  = "SELECT distance as max from graph_path, term WHERE graph_path.term2_id=term.id \
            	 and term.acc='"+ goid + "' ORDER BY distance desc limit 1;"	
         return self.execute(maxdist)[0][0]
	      
     def execute(self, query): 
         try:       
            self.cursor.execute(query)
            results= self.cursor.fetchall()
	 except: 
	    sys.stderr.write("Query " + query + " failed")
	    import traceback
	    traceback.print_exc()
         return results



     def GO(self):
         if self.opts.p:
	    print "PARENT ID", "\tCHILD TERMS IN GOLIST\t" ,"DATABASE INFO FOR PARENT\n"
            for gonumber in self.golist:         
	        results = self.fetchparents(gonumber)				
		for i in results: 
        	    parent=i[3]
		    if parent in self.parents["count"]: 
		       self.parents["count"][parent].append(gonumber)
		    else: 
		       self.parents["count"][parent] = [gonumber]
        	    if parent not in self.parents["info"]:		     	   
		       self.parents["info"][parent] = i
	     
	    for keys,values in sorted(self.parents["count"].items(), key=itemgetter(1), reverse=True): 
		print keys, self.parents["count"][keys], self.parents["info"][keys]	      		       		       
         if self.opts.d:
	    for gonumber in self.golist:
	         results = self.getdist2root(gonumber)
		 if len(results)> 1:	           
		   print results[1][1] + "\t" + results[1][2] + "\t" + results[1][3] + "\t" + str(results[0][7]) + "," + str(results[1][7])
		 elif len(results) < 2 and len(results) > 0:
		   print results[0][1] + "\t" + results[0][2] + "\t" + results[0][3] + "\t" + str(results[0][7])
		 else: 
		   print "NA\tNA\t" + gonumber + "\tNA"
	 if self.opts.m:
	    for gonumber in self.golist:
	       results = self.getmaxdist(gonumber)	   
               print gonumber, results
	 
	 if self.opts.g: 
	    for gonumber in self.golist: 
		 results = self.getNODEinfo(gonumber)
		 if len(results)> 0:
		   print results[0][1] + "\t" + results[0][2] + "\t" + results[0][3] 
		 else: 
		   print "NA\tNA\t" + gonumber
	 

	 
     def mainthing(self):
         self.golist=[]
         for line in file(self.opts.i[0], "r"):
	     gonumber=line.strip("\n")
	     self.golist.append(gonumber)
	 self.GO()  
	     
if __name__ == '__main__':

     try:
         obj = GOFetch()
	 obj.opts = obj.parser.parse_args(sys.argv[1:])
	 if obj.opts.verbose: 
	   print "\n### " + sys.argv[0] + " initialized at " + obj.timestarted  
	   print "### OPTIONS : " + str(obj.opts)
         obj.mainthing()
     except IOError as i:
    	 print "I/O error({0}): {1}".format(i.errno, i.strerror)
     except Exception,e: 
	 print str(e)
	 import traceback
	 traceback.print_exc()



#import urllib
#from xml.etree import cElementTree as ElementTree
#
#def get_go_name(go_id):
#    #get the GO entry as XML
#    xml = urllib.urlopen("http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:"+go_id+"&format=oboxml")
#    #open in cElementTree, for fast XML parsing
#    for event, element in ElementTree.iterparse(xml):
#        #need to make sure we are getting the name contained within the 'term' entry
#        if element.tag == 'term':
#            for child in element.getchildren():
#                #this is the name of the GO ID in the URL above.
#                if child.tag == 'name':
#                    return child.text
#
##print get_go_name('0006915')
