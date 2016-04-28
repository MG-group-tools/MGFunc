from __future__ import division
import argparse
from Bio import SeqIO
from datetime import datetime as dt
import time
import os
import sys
import gzip

class main:
	def __init__(self):
		self.start = time.time()
		self.d_ = dt.today()
		self.timestarted = self.d_.strftime("%d-%m-%Y %H:%M:%S")
		self.parseArgs()
	def parseArgs(self):###GETTING ARGUMENTS FROM COMMANDLINE###
		 parser = argparse.ArgumentParser(prog="swiss2tab",usage="swiss2tab.py -i <input UNIPROT> -o <output-file>",epilog="Example: python2.7 swiss2tab.py -i uniprot_sprot.dat -o uniprot_sprot.tab\n\nWritten by Kosai+Asli, OCT 2013. Last modified MAY 2014.",description="Desctription: Extracts AC,ID,DE,GN,Taxonomy,AC(cession),Organism,ncbi_taxID,GO-term,KEGG-id from STOCKHOLM-formatted file and converts it to tabular-format")
		 parser.add_argument("-i",metavar="database", help="STOCKHOLM-formatted database",nargs=1,required=True)
		 parser.add_argument("-o",metavar="OUTPUT NAME",help="output-name, put the whole output name, fx '-o uniprot.dat.tab'",nargs=1,required=True)
#		 parser.add_argument("-q","--quiet",help="Quiet-mode, suppresses all stdout output. Write \"-q\" with no arguments in commandline. Default is off.",action="store_true")
		 parser.add_argument("-v",help="Verbose. Prints out progress and details to stdout output. Write \"-v\" with no arguments in commandline. Default is off.",action="store_true")
#		 return parser.parse_args(), parser
		 self.parser = parser

	def makeTAB(self):
			fid = self.gzipopen(self.args.i[0]) #input_database
			fout = open(self.args.o[0],"w") #output_tab-file-name
			dbfile = os.popen("grep \"ID   \" "+self.args.i[0] + " | wc -l")
			ctot = dbfile.read()
			dbfile.close()
			ctot = int(ctot.split(" ")[0])
			rangelist = range(0,ctot,10000)
			timeEST = ctot*17/536489
			self.printer("Estimated time usage: "+str(round(timeEST,1))+" minutes ("+str(round(timeEST/60,1))+" hours)\n")
			input_seq_iterator = SeqIO.parse(fid, "swiss")
			fout.write("AC(name)\tID\tDE\tGN\tTaxonomy\tAccession\tOrganism\tncbi_taxID\tGO_term\tKEGG_id\n")
			rowstring = ""
			c = 0
			for record in input_seq_iterator:
					if record.name:
							rowstring += record.name+"\t"
					else:
							rowstring += "N/A\t"
					if record.id:
							rowstring += record.id+"\t"
					else:
							rowstring += "N/A\t"
					if record.description:
							rowstring += record.description+"\t"
					else:
							rowstring += "N/A\t"
					if record.annotations:
							if 'gene_name' in record.annotations:
									rowstring += str(record.annotations['gene_name'])+"\t"
							else:
									rowstring += "N/A\t"
							if "taxonomy" in record.annotations:
									rowstring += str(record.annotations["taxonomy"])+"\t"
							else:
									rowstring += "N/A\t"
							if "accessions" in record.annotations:
									rowstring += str(record.annotations['accessions'])+"\t"
							else:
									rowstring += "N/A\t"
							if "organism" in record.annotations:
									rowstring += str(record.annotations['organism'])+"\t"
							else:
									rowstring += "N/A\t"
							if "ncbi_taxid" in record.annotations:
									rowstring += str(record.annotations['ncbi_taxid'])+"\t"
							else:
									rowstring += "N/A\t"

                                        KEGG = []
                                        GO = []
                                        if record.dbxrefs:
                                                for el in record.dbxrefs:
                                                        if el[0:3] == "GO:":
#                                                               rowstring += el[3:]+";"
                                                                GO.append(el[3:])
                                                        if el[0:5] == "KEGG:":
                                                                KEGG.append(el[5:])
					if not KEGG:
#                                               rowstring += "N/A"
                                                KEGG.append("N/A")
					if not GO:
                                                GO.append("N/A")
					go = ";".join(GO)
					kegg = ";".join(KEGG)
					rowstring += go + "\t" + kegg
					fout.write(rowstring+"\n")
					rowstring = ""
					c += 1
					if c in rangelist or c==1:
							self.printer("FINISHED "+str(c)+" ENTRIES out of "+str(ctot)+"\n")
							sys.stdout.flush()
			self.printer("FINISHED "+str(c)+" ENTRIES out of "+str(ctot)+"\n")
			fid.close()
			fout.close()
			self.indextab()

	
	def printer(self,string): #surpressing output print if -q (quiet) is on
#		if not self.args.quiet:
		if self.args.v:
			print string,
	


        def indextab(self):
                fid = open(self.args.o[0],"r")
                fout = open(self.args.o[0]+".indexed","w")
		line = fid.readline()
                while 1:
                        start = fid.tell()
                        line = fid.readline()
                        if not line or not len(line):
#                               stop = fid.tell()
#                               header = line.split("\t")[0]
#                               fout.write(header + "\t" + str(start) + "," + str(stop)+"\n")
                                break
                        stop = fid.tell()
                        header = line.split("\t")[0]
                        fout.write(header + "\t" + str(start) + "," + str(stop)+"\n")
                fout.close()
                fid.close()




	def gzipopen(self,fileID):
		if fileID[-3:] == ".gz":
			return gzip.open(fileID)
		else:
			return open(fileID,"rU")


	def mainthing(self):
#		self.printer("Cluster2Fasta initialized at"+str(self.timestarted)+"\n")
		self.makeTAB()
		timeused = (time.time() - self.start) / 60
		self.printer("### Time used: "+str(round(timeused)) + " min ("+str(round(timeused/60,1))+" hours)\n")


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



#	myclass = main()
#	myclass.args = myclass.parser.parse_args(sys.argv[1:])
#	myclass.mainthing()



'''
handle=open(swissfilename, "rU")
input_seq_iterator = SeqIO.parse(handle, "swiss")
for record in input_seq_iterator:
        print record.id, record.name, record.description,record.annotations["taxonomy"],record.annotations['accessions'], record.annotations['ncbi_taxid'], record.annotations['organism'], record.annotations['gene_name']

handle.close()
'''


######################
'''
INPUT:
Extracts AC,ID,DE,GN,Taxonomy,AC(cession),Organism,ncbi_taxID,GO_term,KEGG-id
from STOCKHOLM-formatted file and converts it to tabular-format


OUTPUT:
Tabular form of a stockholm-formatted file, where each line is
an entry.


OPTIONS LIST:
-i database: STOCKHOLM-formatted database
-o OUTPUT NAME: output-name, tab-formatted
-q quiet: Quiet-mode, suppresses all stdout output. Write "-q" with noarguments in commandline. Default is off.
'''

