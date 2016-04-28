#!/usr/bin/env python2.7

from Bio import Phylo
#from Bio.Phylo import PhyloXML
import argparse
from datetime import datetime as dt
import time
import sys
import commands
import os


class main:
	def __init__(self):
		self.start = time.time()
		d_ = dt.today()
		self.timestarted = d_.strftime("%d-%m-%Y %H:%M:%S")
		self.parseArgs()
	def parseArgs(self):
		parser = argparse.ArgumentParser(prog="treenames2.py",usage="treenames2.py -i header-dict -t mytree.dnd -format phyloxml -v",epilog="Written by Asli and Kosai, nov 2013. Last modified mar 2014.")
		parser.add_argument("-i",help="Header-change file(s)",nargs="*",required=True)
		parser.add_argument("-t",help="tree file(s)",nargs="*",required=True)
		parser.add_argument("-u",help="tab-formatted Uniprotfile",nargs="?")
		parser.add_argument("-n",help="Node-names, can be \'D\' (Description/function), \'G\', (Gene name, DEFAULT), \'T\', (taxonomy), or \'A\', (All, tax,descript,genename).",nargs="?",default="G")
		parser.add_argument("-o",help="Output-name. If not specified, working-dir will be used",nargs=1)
		parser.add_argument("-format",help="format, can be \'phyloxml\',\'newick\',etc.",nargs=1)
		parser.add_argument("-p",help="prints trees out in STDOUT",action="store_true")
		parser.add_argument("-c",help="Deletes all temporary files",action="store_true")
		parser.add_argument("-l",help="Log file",nargs=1)
		parser.add_argument("-v",help="prints date and time out",action="store_true")
#		return parser.parse_args(),parser
		self.parser = parser
        def printert(self,string): 
            if self.args.v:
               if self.args.l: 
	          logfile= self.args.l[0]
	       else: 
	          logfile = self.args.t[0] + ".trenames.log"
	       
	       with open(logfile, "a") as log: 
	          log.write(string)
        def printer(self,string): 
            if self.args.v:
 	       print string
	def dictMaker(self,argsi,argso):
		#fh=open("../headerchange.CL_325", "r")
		fh = open(argsi,"r")
		headerids = argso +".geneids"
		fo = open( headerids,"w") #temp uniprot ID's file
		headerdict={}
		uniprot = ""
		for line in fh:
			line_s = line.split("\t")
			uniprot = line_s[1].strip("\n").split(":")[1]
			headerdict[line_s[0]]=uniprot
			if self.isUniprot(uniprot) == 1:
				fo.write(uniprot+"\n")
		fh.close()
		fo.close()
		return headerdict,headerids
	
	def makeTab(self,fidname,argsU,args):
		fout = fidname+".uniprotTab"
		self.printert(fidname)
		if os.path.exists(fout):
			if args.c:
				os.remove(fidname)
			return fout
		else:
			my_cmd = "fgrep -w -f "+fidname+" "+argsU+" > "+fout
			no_output = commands.getstatusoutput(my_cmd)
			if no_output[0] != 0:
				self.printer("Could not parse uniprot-file for %s. Probably has no uniprot hits\n" % fidname)
			if args.c:
				os.remove(fidname)
			return fout
	
	def makeUniprot(self,argsU,kd): #argsU has changed to a subset of the whole uniprot tab file
		d = 0		
		fh =open(argsU,"rb")
		myD = {}
		fh.seek(0)
		if fh.readline()[0:8] == "AC(name)": #hack..
			sup = fh.readline() #sup stands for suppression.
		else:
			fh.seek(0)
		temp = []
		descript = ""
		organism = ""
		myasciirange = range(65,91) + range(97,123) + [32] + range(48,58)
		d1 = 1
		d2 = 1 #d1 and d2 are identifiers to make nodes unique
		for line in fh:
			line = line.split("\t")
			#Something that changes line[2] to remove "Full="
			##NOTE: the Uniprottab file has changed, so has the characters here
			descript = str(d1) +" "+ self.removeWeirdChars(line[4].strip(";").split("=")[1],myasciirange)
			#print descript
			d1 += 1
			'''descript = re.sub("\:","_",line[2].strip(";").split("=")[1])
			descript = re.sub("\@","_",descript)
			descript = re.sub("\;","_",descript)'''
			temp = line[6][2:-2].split("', '")
			organism = " ".join(line[8].split(" ")[0:2])
			myD[line[2]] = [descript,str(d2)+ " " + temp[0][0].capitalize()+"_"+self.smallTaxonomy(temp)+"_"+organism] #Use shortenKingdom(temp[0]) instead of temp[0]		
			d2 += 1
		fh.close()
		if len(myD) != 0:
			d = 1
		if self.args.c:
			os.remove(argsU)
		#print myD
		return d,myD
	
	def openTree(self,argst):
		tree = Phylo.read(argst,"newick")
	#	tree = Phylo.read('example.dnd', 'newick')
		n2tree = tree.as_phyloxml()
		return tree, n2tree
	
	def updateTreeNames(self,n2tree,headerdict,uniDict,argsn):
		for clade in n2tree.get_terminals():
		        try:
			   genename = headerdict[clade.name]
			
			   if self.isUniprot(genename) == 0:
				   clade.name=genename
			   elif argsn == "D":
				   clade.name = genename + " | " + uniDict[genename][0]
			   elif argsn == "G":
				   clade.name = genename
			   elif argsn == "T":
				   clade.name = uniDict[genename][1]
			   elif argsn == "A":
				   clade.name = uniDict[genename][1]+" | "+uniDict[genename][0]+" | "+genename
		        except:
			   print self.args.i,clade.name + "\n"
			   #sys.exit()
		
		return n2tree
	
	def writeTree(self,outname,n2tree,form):
		fout = open(outname,"w")
		status = Phylo.write(n2tree, fout,form)
		fout.close()
		return status
	
	def isUniprot(self,myID):
		if len(myID.split("_")) == 2:
			return 1
		else:
			return 0
	
	def removeWeirdChars(self,s,myrange):
		snew = ""
		for ch in s:
			if ord(ch) in myrange:
				snew+=ch
			else:
				snew+=" "
		return snew
	
	def smallTaxonomy(self,mylist): #Tax = taxonomy, G-Unit!!
		if len(mylist) <2:
			return ""
		else:
			return mylist[1]
			
	
	def kingdomDic(self):
		myD = {"Bacteria":"B","Viruses":"V","Archaea":"A","Eukaryota":"E","unclassified sequences":"U","other sequences":"O"}
		return myD
	
	#from Bio.Phylo import PhyloXML
	#n2tree=tree.as_phyloxml()
	
	
	#for clade in n2tree.get_terminals():
	#	clade.name = headerdict[clade.name]
	
	
	def validate(self,n2tree):
		for clade in n2tree.get_terminals():
			key = clade.name
			self.printert(key)
	
	def printTree(self,n2tree,argso):		
	    with open(argso+".treefig","w") as drawascii:
	         Phylo.draw_ascii(n2tree, file=drawascii)

	def mainthing(self):
#		if self.args.v: #prints start time
#			self.printer("\n**treenames2.py initialized at " + self.timestarted + "**\n")
		uniDict = {} #uniprot ID's saved here
		self.out = []
		if self.args.i and self.args.t and len(self.args.i) != len(self.args.t):
			sys.exit("ERROR: Mismatched number of header-files ("+str(len(self.args.i))+") and tree-files ("+str(len(self.args.t))+")!\n")
		elif self.args.i and not self.args.t:
			sys.exit("ERROR: Mismatched number of header-files ("+str(len(self.args.i))+") and tree-files ("+"NONE"+")!\n")
		elif self.args.t and not self.args.i:
			sys.exit("ERROR: Mismatched number of header-files ("+"NONE"+") and tree-files ("+str(len(self.args.t))+")!\n")
		else:
			kd = self.kingdomDic()
			for c in range(0,len(self.args.i)): #looking through multiple files, pointed by the variable c
				#print self.args.i
				if self.args.o:
					self.out.append(self.args.o[0] + "." + self.args.n + "." + self.args.format[0])
				else:
					if self.args.t[c][0] != ".":
					    self.out.append(self.args.t[c].split(".")[0] + "." + self.args.n + "." +self.args.format[0])
			                else: 
					    self.out.append(".".join((self.args.t[c].split(".")[0:2])) + "." + self.args.n + "." +self.args.format[0])
				headerdict,fidname = self.dictMaker(self.args.i[c], self.out[c])
				
				if self.args.u:
					dic_status = 0
					fout_tab = self.makeTab(fidname,self.args.u,self.args)
					dic_status,uniDict = self.makeUniprot(fout_tab,kd)
				#if dic_status == 0:
					#sys.exit("There was a problem in creating a dictionary from the Uniprot tab-file"+self.args.i[c].split("/")[-1]+"\n")
				
				tree,n2tree = self.openTree(self.args.t[c])
				
				
				if not self.args.n: #in case user just wrote "-N" withouth arguments. Leave it for safety
					self.args.n = "G"
				
				
				n2tree = self.updateTreeNames(n2tree,headerdict,uniDict,self.args.n)	
			
				#if args.test: #prints tree notes out in STDOUT
				#	validate(n2tree)
				
								
				if not self.args.format:
					status = self.writeTree(self.out[c],n2tree,"phyloxml")
				else:
					status = self.writeTree(self.out[c],n2tree,self.args.format[0])
				if status != 1:
					self.printer("ERROR! Could not write tree for cluster :"+str(self.args.t[c])+"\n")
				if self.args.p:
					self.printTree(n2tree,self.out[c]) #prints tree to STDOUT
				if len(self.args.i) == 1:
					continue
				elif self.args.v:
					self.printer("*Finished tree "+self.args.t[c].split("/")[-1]+" [with headers: "+self.args.i[c].split("/")[-1]+"]\n")
		timeused = (time.time() - self.start) / 60
		if self.args.v: #prints end-time
			self.printert("### Time used: "+str(round(timeused*60)) + " seconds ("+str(round(timeused)) + " min)\n")
		



if __name__ == "__main__":
    try:
        myclass = main()
        myclass.args = myclass.parser.parse_args(sys.argv[1:])
        myclass.printer("\n### "+sys.argv[0]+" initialized at "+ myclass.timestarted + "\n")
        myclass.printer("### OPTIONS: "+str(myclass.args)+"\n")
        myclass.mainthing()
    except IOError as i:
        print "I/O error({0}): {1}".format(i.errno, i.strerror)
	sys.exit()
    except Exception,e:
        print str(e)
        import traceback
        traceback.print_exc()
	sys.exit()




#	myclass = main()
#	myclass.args = myclass.parser.parse_args(sys.argv[1:])
#	myclass.mainthing()



#'member to do a version of os.exists(uniprot_small_tab_file)

#"""
#Written by Asli and Kosai
#"""
#
#
############################
#'''
#INPUT
#Takes a tree in newick-format with headerchange-files (seqn->gene-ID) and converts
#the tree to any given format (default phyloxml). If you want other than gene-names, fx
#taxonomy, you need to give a tab-formatted uniprot-file.
#
#OUTPUT 	
#A new tree in a user-specified format. You can also print the tree in STDOUT if you want.
#
#OPTIONS:
#-i: Header-change file(s)
#-t: .dnd tree file(s)
#-U: tab-formatted Uniprotfile
#-N: Node-names, can be 'D' (Description/function), 'G', (Gene name, DEFAULT), 'T', (taxonomy), or 'A', (All, tax,descript,genename).
#-format: format, can be 'phyloxml','newick',etc.
#-p: prints trees out in STDOUT
#-c: Deletes all temporary files
#-v: prints date and time out
#
#'-i ./MGFuncrunSP/changeheaders/CL_100.headertable -t ./MGFuncrunSP/trees/CL_100.new.njtree -u ./MGFuncrunSP/parseclusterinfo/test_list2.uniprottab -n A -format newick -c -v'
#'''
