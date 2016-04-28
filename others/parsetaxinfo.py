#!/usr/bin/env python

import sys, gzip
import re, string
import optparse
import os
from Bio.Blast import NCBIStandalone
from operator import itemgetter, attrgetter
sys.path.append('/home/projects5/pr_53035/people/asli/bin/lib/python2.7/site-packages')
import networkx as nx
import kyotocabinet as kc
from os.path import basename

def getuniprotresults(clusterfile,uniprotfile):
    #print "starting function"
    
    fh = open(uniprotfile, "r")
     
    
    if os.path.exists( "../blastdecideresults_normalblast/" + basename(uniprotfile) + '.blastdict.kch'): 
       db = kc.DB()
       db.open("../blastdecideresults_normalblast/" + basename(uniprotfile) + '.blastdict.kch', kc.DB.OREADER)		       
    else:                                        
       db=kc.DB()
       db.open(basename(uniprotfile) + '.blastdict.kch',  kc.DB.OWRITER | kc.DB.OCREATE)
       
       for line in fh:			     			     

    	   query=line.split("\t")[0]		     
    	   unipid=str(line.split("\t")[1]) + str(",")
	   #print query,unipid
    	   if int(db.check(query))>0:		     
    	      db.append(query,unipid)		     
    	   else:					     
    	      db.set(query,unipid) 		     



       fh.close()

    fc=open(clusterfile, "r")							        
    cl=open("Koalagene_taxid","w")
    for l in fc:
       #print l	
       #print "XXXXXXXXXXXXXXXX"							        
        							        
       members_text=l.strip("\n")						        
       members=members_text.split("\t") 					        
       			        
       					        
       for i in members:
           uniprotids=set()							        
    	   if int(db.check(i))>0:						        
    	      #if len(blastdict[i]) > 0:					        
    		 result = db.get(i).strip(",").split(",")			        
    		 uniprotids=set(result)  						        
    		 					        
           taxids={}
           for ids in uniprotids:
               t=fetchtaxinfo_uni(ids)
               taxinfo=list(eval(t))
               if ids in taxids:
	          taxids[ids]+=1
	       else:
	          taxids[ids]=1
	       
       	       			        
           membertax=sorted(taxids, key=taxids.get, reverse=True)[0]              
           cl.write("%s_%s\n" % (str(i), membertax))
	  
	  
    fc.close()
    cl.close()
    db.close()   					        
   
def fetchtaxinfo_uni(searchid): 
    
    command="grep -w " + searchid + "uniprothits_taxinfo"
    fout2 = os.popen(command)						        
    out2 = fout2.readline().strip("\n").split["\t"][1]  				        
    return out2      

if __name__ == '__main__':


   #try:
       #print sys.argv[1]
       if sys.argv[1]=="yes":
	  clusterfile=sys.argv[2]
	  uniprotfile=sys.argv[3]
	  #print clusterfile, uniprotfile
          getuniprotresults(clusterfile,uniprotfile)
  # except:
   #    print "Usage:", sys.argv[0],"<clusterfile>","<uniprot_blasdecide>"
