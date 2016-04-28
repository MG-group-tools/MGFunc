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
		parser = argparse.ArgumentParser(prog="swiss2fasta.py",usage="swiss2fasta.py -i <input UNIPROT> -o <output fasta-file>",epilog="Example: python2.7 swiss2fasta.py -i uniprot_sprot.dat -o uniprot_sprot.fasta\n\nWritten by Kosai+Asli, OCT 2013. Last modified MAR 2014.",description="Desctription: Extracts sequence and ID\'s from STOCKHOLM-formatted file and converts it to fasta-format")
		parser.add_argument("-i",metavar="database", help="STOCKHOLM-formatted database (UNIPROT)",nargs=1,required=True)
		parser.add_argument("-o",metavar="OUTPUT NAME",help="output-name, fasta",nargs=1,required=True)
#		parser.add_argument("-q","--quiet",help="Quiet-mode, suppresses all stdout output. Write \"-q\" with no arguments in commandline. Default is off.",action="store_true")
		parser.add_argument("-v",help="Verbose. Prints out progress and details to stdout output. Write \"-v\" with no arguments in commandline. Default is off.",action="store_true")
#		self.self.args = parser.parse_self.args()
		self.parser = parser

	def makeFASTA(self):
		fid = self.gzipopen(self.args.i[0]) #input_database
		fout = open(self.args.o[0],"w") #output_fasta-file-name
		dbfile = os.popen("grep \"ID   \" "+self.args.i[0] + " | wc -l")
		ctot = dbfile.read()
		dbfile.close()
		ctot = int(ctot.split(" ")[0])
		rangelist = range(0,ctot,10000)
		timeEST = ctot*17/536489
		self.printer("Estimated time usage: "+str(round(timeEST,1))+" minutes ("+str(round(timeEST/60,1))+" hours)\n")
		input_seq_iterator = SeqIO.parse(fid, "swiss")
		rowstring = ""
		seqstring = ""
		newline = "\n"
		c = 0
		for record in input_seq_iterator:
			if record.name:
				rowstring += ">" + record.name
			else:
				rowstring += ">N/A"
			fout.write(rowstring+newline)
			rowstring=""
			if record.seq:
				seqstring += record.seq
			else:
				seqstring += "N/A"
			fout.write(str(seqstring)+newline)
			seqstring = ""
			c += 1
			if c in rangelist or c==1:
				self.printer("FINISHED "+str(c)+" ENTRIES out of "+str(ctot)+"\n")
		self.printer("FINISHED "+str(c)+" ENTRIES out of "+str(ctot)+"\n")
		fid.close()
		fout.close()

        def gzipopen(self,fileID):
                if fileID[-3:] == ".gz":
                        return gzip.open(fileID)
                else:
                        return open(fileID,"rU")



	def printer(self,string): #surpressing output print if -q (quiet) is on (NEW: shows output if -v is on)
#		if not self.args.quiet:
		if self.args.v:
			print string,

	def mainthing(self):
#		self.printer("swiss2fasta initialized at"+str(self.timestarted)+"\n")
		self.makeFASTA()
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
Takes a stockholm-formatted database such as Uniprot or Swiss-prot, and output
a fasta-file with the sequence for each entry. Each entry's name is the
gene ID from the database.

OUTPUT:
Fasta-file with sequences for the uniprot.


OPTIONS LIST:
-i (database): STOCKHOLM-formatted database (UNIPROT)
-o (OUTPUT NAME): output-name, fasta
-q (quiet): Quiet-mode, suppresses all stdout output. Write "-q" with no arguments in commandline. Default is off.
'''


