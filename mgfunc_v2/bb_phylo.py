#!/usr/bin/env python
# Created: Tue Jul  6 08:49:28 1999
# Last changed: Time-stamp: <11/03/10 15:10:42 thomas>
# Thomas.Sicheritz@molbio.uu.se, http://evolution.bmc.uu.se/~thomas
# File: bb_phylo.py

import string, re
import commands, tempfile
import os, sys  # os.system, sys.argv
sys.path.insert(0, '.')
from bb_utils import myunzip, myzip
from optparse import OptionParser
#from templabels import run_paup_with__changed_labels

class phylo_tool:
    def __init__(self, alnfile=None, outfile=None, program=None, create_only = 0):
        self.logfile = None
        if alnfile:
            if alnfile[-3:] == '.gz':
                self._alnfile = alnfile
                myunzip(alnfile)
                alnfile = alnfile[:-3]
                was_zipped = 1
            else:
                was_zipped = 0
                
            #if not os.path.exists(alnfile): return
            size = os.stat(alnfile)[6]
            if size == 0: return
            self.infile = alnfile
            self.outfile = outfile
            self.program = program
            self.nr_taxas = self.count_taxas()
            print '%s has %d taxas' % (alnfile, self.nr_taxas)
            if self.nr_taxas <4: return
            self.make_tree(create_only)
            if was_zipped: myzip(alnfile)
            
    def count_taxas(self):
        fid = open(self.infile,'r')
        fid.readline()
        fid.readline()
        fid.readline()
        taxas = 0
        line = fid.readline()
        while not line.startswith(' '):
            line = fid.readline()
            if not line: break
            taxas += 1
        fid.close()
        return (taxas)
    
    def make_tree(self, create_only = 0):
        argv = string.split(self.program)
	print 'ARGV', argv
        if argv[0] == "clustalw":
            com = "clustalw " + self.infile + " -bootstrap -outputtree=phylip"
            tmpoutfile = self.infile[:-4] + ".phb"
            if not create_only:
		print com
                os.system(com)
                os.rename(tmpoutfile, self.outfile)

        elif argv[0] == "paup":
	    print 'pauping'
            basename = os.path.splitext(os.path.basename(self.infile))[0]
            nexusfile = basename + ".nexus"
            pmethod = getattr(self,'paup_neighbor')
	    print 'with', pmethod
            if len(argv) > 1: 
                if argv[1] =='jk':
                    pmethod = getattr(self,'paup_jackknife')
            tmpoutfile = pmethod()
            com = "paup " + nexusfile + " -u -n"
            print 'PAUPCOM', com
            paup_options = " -u -n"
            if not create_only:
                run_paup_with__changed_labels(nexusfile, paup_options)
                #os.system(com)
                try:
                    print >> sys.stderr, 'renaming <%s> to <%s>' % (tmpoutfile, self.outfile)
                    os.rename(tmpoutfile, self.outfile)
                except:
                    ""
                #self.FixNewickTreeStupidPaupARGHHHHHHH(basename + '.paup.tree', nexusfile)
                
        elif self.program == "phylip":
            tmpoutfile = self.phylip_tree()
            os.rename(tmpoutfile, self.outfile)
            
        elif self.program == "mrbayes":
            tmpoutfile = self.mrbayes()
            basename = os.path.splitext(os.path.basename(self.infile))[0]
            nexusfile = basename + ".nexus"
            com = 'mb3 %s' % nexusfile
            print 'MRBAYES command:', com
            if not create_only:
                os.system(com)
                try:
                    os.rename(tmpoutfile, self.outfile)
                except:
                    ""
            
    def paup_jackknife(self):
        ALNFILE = os.path.abspath(self.infile)
        basename = os.path.splitext(os.path.basename(ALNFILE))[0]
        NEXUSFILE = basename + ".nexus"
        PIRFILE = basename + ".pir"
        LOGFILE = basename + ".logjk"
        JKFILE = basename + ".jackknife.trees"
        NJFILE = basename + ".neighbour.trees"
        os.system("clustalw -convert -infile=" + ALNFILE + " -output=pir -outfile=" + PIRFILE)
        #print "clustalw " + ALNFILE + " -convert -output=pir -outfile=" + PIRFILE

        mall = '''tonexus fromfile=PIRFILE tofile = NEXUSFILE format = pir datatype = protein interleave = yes replace = yes;
        
        BEGIN paup;
        LOG file=LOGFILE REPLACE;
        
	set CRITERION=parsimony;
        SET TAXLABELS=full;
  
	JACKKNIFE
        keepall = yes
	PCTDELETE = 37
	JSEED = 13    
	NREPS = 100   
			
	[SAVE FILE=jackknifetrees
	SEarch=Heuristic/
	[HSEARCH            
	NOSTATUS            
	ADDSEQ = RANDOM     
			    
	NREPS = 1           
	RSEED = 13          
	
	NORSTATUS           
	
	SWAP = NNI          
	NCHUCK = 10
	
	CHUCKLEN = 10

        NOCOLLAPSE
	;                    [markerar avslutad kommandorad]
        contree all /le50 = yes;
        savetrees FMT=phylip FILE=JKFILE REPLACE=YES;
        [showdist;]
        describetrees all /brlens = yes;
        quit;
        endblock;
        '''

        for i in ['ALNFILE', 'PIRFILE', 'LOGFILE', 'JKFILE', 'NJFILE', 'NEXUSFILE']:
            print i, eval(i)
            mall = string.replace(mall, i, eval(i))
            
        tmpfile = tempfile.mktemp()
        fid = open(tmpfile, 'w+')
        lines = string.split(mall, '\n')
        fid.write(lines[0] + '\n')
        fid.close()
        os.system("paup %s  -u -n" % tmpfile)
        os.remove(tmpfile)
        
        txt = open(NEXUSFILE).read()
        txt = string.replace(txt, basename, "'%s'" % basename)
        fid = open(NEXUSFILE, 'w+')
        fid.write(txt)
        #fid = open(NEXUSFILE, 'a')
        fid.write("\n" + mall[len(lines[0]):] + "\n")
        fid.close()
        try:
            os.remove(PIRFILE)
        except:
            ""
        self.logfile = LOGFILE
        self.treefile = JKFILE
        return JKFILE

    def paup_neighbor(self):
        #ALNFILE = os.path.abspath(self.infile) 
	ALNFILE = self.infile
        basename = os.path.splitext(os.path.basename(ALNFILE))[0]
        NEXUSFILE = basename + ".nexus"
        PIRFILE = basename + ".pir"
        LOGFILE = basename + ".lognj"
        NJFILE = basename + ".neighbour.trees"
	print 'ALNFILE', ALNFILE
	print commands.getoutput('ls -l %s' % ALNFILE)
	print os.getcwd()
        com = "clustalw  -convert -infile=" + ALNFILE + " -output=pir -outfile=" + PIRFILE
        print com
        os.system(com)
                

        mall = '''tonexus fromfile=PIRFILE tofile = NEXUSFILE format = pir datatype = protein interleave = yes replace = yes;
        
        BEGIN paup;
        LOG file=LOGFILE REPLACE;
        SET CRITERION=distance;
        SET TAXLABELS=full;
        bootstrap search=nj keepall = yes;
        contree all /le50 = yes;
        savetrees from=1 to=1 FMT=phylip FILE=NJFILE REPLACE=YES;
        [showdist;]
        describetrees all /brlens = yes;
        quit;
        endblock;
        '''

        for i in ['ALNFILE', 'PIRFILE', 'LOGFILE', 'NJFILE', 'NEXUSFILE']:
            print i, eval(i)
            mall = string.replace(mall, i, eval(i))
            
        tmpfile = tempfile.mktemp()
        fid = open(tmpfile, 'w+')
        lines = string.split(mall, '\n')
        fid.write("#NEXUS\n" + lines[0] + '\n')
	
        fid.close()
        os.system("paup %s  -u -n" % tmpfile)
        os.remove(tmpfile)
        
        txt = open(NEXUSFILE).read()
        txt = string.replace(txt, basename, "'%s'" % basename)

        fid = open(NEXUSFILE, 'a')
        fid.write("\n" + mall[len(lines[0]):] + "\n")
	
        fid.close()
	
	os.system("paup %s  -u -n" % NEXUSFILE)
       	try:
            os.remove(PIRFILE)
        except:
            ""
        
	self.logfile = LOGFILE
        self.treefile = NJFILE
        return NJFILE

    def mrbayes(self):
        ALNFILE = os.path.abspath(self.infile)
        basename = os.path.splitext(os.path.basename(ALNFILE))[0]
        NEXUSFILE = basename + ".nexus"
        PIRFILE = basename + ".pir"
        LOGFILE = basename + ".logmb"
        MBFILE = basename + ".mrbayes"
        MBFILEP = basename + ".mrbayes.p"
        MBFILET = basename + ".mrbayes.t"
        
	print 'ALNFILE', ALNFILE
	print commands.getoutput('ls -l %s' % ALNFILE)
	print os.getcwd()
        com = "clustalw  -convert -infile=" + ALNFILE + " -output=pir -outfile=" + PIRFILE
        print com
        os.system(com)


        # remove amino acid codes which mrbayes won't accept
        fid = open(PIRFILE)
        txt = fid.read()
        fid.close()

        transtable = string.maketrans('XxZz', '----')
        fod = open(PIRFILE, 'w+')
        for line in txt.split('\n'):
            if line and line[0] != '>':
                line = line.translate(transtable)
            fod.write(line+'\n')

        fod.close()

	try:
	    NCGEN = self.NCGEN
	except:
	    NCGEN = 600000
        #NCGEN = 1000

        PRINTFREQ = int(NCGEN/100)
        SAMPLEFREQ = 100
        BURNIN = int(NCGEN/(3*SAMPLEFREQ))


        print 'NCGEN', NCGEN
        print 'PRINTFREQ', PRINTFREQ
        print 'SAMPLEFREQ', SAMPLEFREQ
        print 'BURNIN', BURNIN
        
        mall = '''tonexus fromfile=PIRFILE tofile = NEXUSFILE format = pir datatype = protein interleave = yes replace = yes;
        

        BEGIN mrbayes;
            set autoclose=yes;
            log start filename=LOGFILE replace;
            prset aamodelpr=fixed(wag);
            lset rates=invgamma Ngammacat=4; 
            mcmc ngen=%d printfreq=%d samplefreq=%d nchains=4 savebrlens=yes
            startingtree=random filename=MBFILE;

            log start filename=LOGFILE append;        
            sumt filename=MBFILET burnin=%d;

            quit;
        endblock;
        ''' % ( NCGEN, PRINTFREQ, SAMPLEFREQ, BURNIN)

        for i in ['ALNFILE', 'PIRFILE', 'LOGFILE', 'MBFILET', 'MBFILEP', 'MBFILE', 'NEXUSFILE']:
            print i, eval(i)
            mall = string.replace(mall, i, eval(i))
            
        tmpfile = tempfile.mktemp()
        print 'TMPFILE', tmpfile
        print 'MALL', mall
        fid = open(tmpfile, 'w+')
        lines = string.split(mall, '\n')
        fid.write(lines[0] + '\n')
        fid.close()
        os.system("paup %s  -u -n" % tmpfile)
        os.remove(tmpfile)
        
        txt = open(NEXUSFILE).read()
        txt = string.replace(txt, basename, "'%s'" % basename)

        fid = open(NEXUSFILE, 'a')
        fid.write("\n" + mall[len(lines[0]):] + "\n")
        fid.close()
        self.logfile = LOGFILE
        self.treefile = MBFILE

        if os.path.exists(MBFILET): os.remove(MBFILET)
        if os.path.exists(MBFILEP): os.remove(MBFILEP)
        return MBFILE
    
    def phylip_tree(self):
        tempfile.tempdir = '.'
        tmpdir = tempfile.mktemp()
        os.mkdir(tmpdir)
        
        basename = os.path.splitext(os.path.basename(self.infile))[0]
        phylipfile = tmpdir + '/infile'
        com = "clustalw  -convert " + self.infile + " -output=phylip -outfile=" + phylipfile
        os.system(com)

        os.chdir(tmpdir)

        os.system("echo y > yes")
        os.system("echo 37 >> yes")
        os.system("seqboot < yes")

        os.rename("outfile", "infile")
        os.system('cp infile infile1')
        os.system("echo p > yes")
        os.system("echo p > yes")
        os.system("echo m >> yes")
        os.system("echo d >> yes")
        os.system("echo 100 >> yes")
        os.system("echo y >> yes")
        os.system("protdist < yes")

        os.rename("outfile", "infile")
        os.system('cp infile infile2')
        os.system("echo m > yes")
        os.system("echo 100 >> yes")
        os.system("echo 37 >> yes")
        os.system("echo y >> yes")
        os.system("neighbor < yes")

        os.rename("treefile", "infile")
        os.system("echo y > yes")
        os.system("consense < yes")
        os.chdir('..')
        os.removedirs(tmpdir)
        return "%s/treefile" % tmpdir

    def remove_nodes(self, treefile, nodes):
        fid = open('retree.answer','w+')
        fid.write('0\nl\n10000000\ny\n%s\nx\nn\n' % treefile)
        fid.close()
        if os.path.exists('intree'): os.remove('intree')
        txt = commands.getoutput("retree < retree.answer")
        lines = string.split(txt,'\n')
        node_numbers = []
        for node in nodes:
            rxnode = string.replace(node," ","[ _]")
            rxnode = string.replace(rxnode,"_","[ _]")
            for line in lines:
                m = re.search("-([0-9]+):" + rxnode, line)
                print rxnode, line
                if m:
                    node_numbers.append(int(m.group(1)))

        str = '0\nl\n10000000\ny\n%s\n' % treefile
        
        for nr in node_numbers:
            str = str + 'r\n%d\n\n' % nr
        str = str + 'w\nu\nx\n\n'

        fid = open('retree.answer','w+')
        fid.write(str)
        fid.close()
        os.system("retree < retree.answer")


    def AsciiTree(self, pyphy, treefile, db_index):
        sys.path.insert(0, pyphy)
        from bb_tree import BasicTree
        if not self.logfile: return
        txt = open(self.logfile).read()
        rx = re.compile('((Bootstrap|Jackknife) 50% majority-rule consensus tree.+)Bipartitions found in one or more trees', re.S|re.M)

        ascii_tree = rx.search(txt).group(1)



        tree = BasicTree(treefile, isfile = 1)

        #print 'BBBBB', self.treefile, tree.hash.keys()

        for i in tree.hash.keys():
            line = '%s,%s,%s' % (i, db_index.Get_Organism(i), db_index.Get_Gene(i))
            ascii_tree = ascii_tree.replace(i, line)
            
        outfile = self.outfile + 'ascii'

        fid = open(outfile,'w+')
        fid.write(ascii_tree + '\n')
        fid.close()
        
    def FixNewickTreeStupidPaupARGHHHHHHH(self, newick_file, nexus_file):
        base = os.path.split(newick_file)[-1].split('.')[0]
        print >> sys.stderr, "Fixing Paup's newick export Bug/Feature of chopping long names", base,

        if os.path.exists(newick_file + '.changed'):
            print >> sys.stderr, 'already changed'
            return

        newicktree = open(newick_file).read()

        record = False
        taxa = []
        for line in open(nexus_file):
            if line.find('ntax=') > -1:
                ntax = line.split('ntax=')[1].split()[0]
            if line.find('Matrix') > -1:
                record = True
                continue

            if record:
                if line.find(';') > -1:
                    record = False
                    break

                fields = line.split()
                if len(fields) >1:
                    taxa.append(fields[0])

        for name in taxa:
            if len(name) > 10:
                shortname = name[:10]
                print >> sys.stderr, 'changing <%s> to <%s>' % (shortname, name)
                newicktree = newicktree.replace(shortname, name)

        os.rename(newick_file, newick_file + '.changed')
        fid = open(newick_file, 'w+')
        print >> fid,  newicktree
        fid.close()
        

        
if __name__ == '__main__':
    usage = "%prog [options] file (or - for stdin)\n"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default = 0)
    parser.add_option("-c", "--create_only", action="store_true", dest="create_only", default = False)
    parser.add_option("-t", "--tool", action="store", dest="tool", default = 'mrbayes')
    (options, args) = parser.parse_args()
    verbose = options.verbose
    files = args

    for file in files:
        outfile = os.path.splitext(os.path.split(file)[-1])[0] + '.phy'
        test= phylo_tool(file, outfile, options.tool, create_only = options.create_only)
        #test= phylo_tool(file, "dum.phy", options.tool, create_only = options.create_only)
#     from getgene import DB_Index
#     db_index = DB_Index()
#     test = phylo_tool()
#     test.logfile = 'logfile'
#     test.outfile = 'out'
#     test.AsciiTree(os.environ['PYPHY'], file, db_index)

    
