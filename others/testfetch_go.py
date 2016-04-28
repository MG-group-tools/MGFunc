#!/usr/bin/env python
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen (asli@cbs.dtu.dk)
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)

import sys
import MySQLdb
import urllib2
from operator import itemgetter 

infile = file(sys.argv[1], 'rU').readlines()
taxa = {}
tax=''
db = MySQLdb.Connect(db='go_termdb', user='asli', passwd='Pass9fy109', host='mysql.cbs.dtu.dk')
cursor=db.cursor()
header=''
taxid=''
for line in infile:
     gonumber = line.strip("\n")
     getinfo = "SELECT * FROM term where acc='"+ gonumber + "';"
     getparents = "SELECT DISTINCT ancestor.*, graph_path.distance, graph_path.term1_id AS ancestor_id \
		       FROM term \
		       INNER JOIN graph_path ON (term.id=graph_path.term2_id) \
		       INNER JOIN term AS ancestor ON (ancestor.id=graph_path.term1_id)\
		       WHERE term.acc='"+ gonumber + "';"	
     cursor.execute(getparents)
     results= cursor.fetchall()
     #print results
     for i in results: 
        print i
     results = ""
     print "\n"
