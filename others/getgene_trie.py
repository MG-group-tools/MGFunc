#!/usr/bin/env python
# Created: Sun Oct 15 16:16:20 2000
# Last changed: Time-stamp: <08/05/21 12:24:44 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: getgene.py

# SPyPhy imporved variant of getgene.py


import string, re
import os, sys, glob
from commands import getoutput
import gdbm, cPickle
import getopt

MACHINE = os.uname()[-1]

class DB_Index:
    temp_entries = {}

    def __init__(self, open = 1, tryoldids = 0, debug = 0, pyphydb=None):
        self.debug = debug
        
        return

    def GetValues(self, id):
        com = '/home/people/thomas/projects/AutoTree/Ftrie_src3/getgene_%s -p up -u /home/people/thomas/projects/AutoTree/pyphynr.triendex %s' % (MACHINE, id)
        if self.debug: print com
        res = getoutput(com)
        return ('hhs_ftrie', res)
            
    def Get(self, id):
        try:
            # try to lookup in temporary application specific dictionary first
            # usefull for temporarily adding values, fake entries or overriding
            #print >> sys.stderr, 'checking temp_entries'
            #print id, len(self.temp_entries.keys())
            txt = self.temp_entries[id]
            return txt
        except:
            pass


        #print >> sys.stderr, 'DBG 1 trying to get', id
        
        try:
            which_db, values = self.GetValues(id)
            #print >> sys.stderr, 'DBG 2 got', values
        except:
            if self.debug:
                print >> sys.stderr, 'Cannot retrieve', id
            return None
        return values
    
    def Get_OrganismKW(self, id):
        entry = self.Get(id)
        if not entry: return None
        OS = ''
        KW = ''
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = OS + string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
            elif line[0:5] == 'KW   ':
                KW = KW + string.strip(line[5:])
                if KW[-1] ==".": KW = KW[0:-1]
            if line[0:2] =="//": break
        if not len(OS): return id
        return '%s\t->\t%s' % (OS, KW)
    
    def Get_OrganismOC(self, id):
        entry = self.Get(id)
        if not entry: return None
        OS = ''
        OC = ''
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = OS + string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
            elif line[0:5] == 'OC   ':
                OC = OC + string.strip(line[5:])
                if OC[-1] ==".": OC = OC[0:-1]
            if line[0:2] =="//": break
        if not len(OS): return id
        return '%s\t->\t%s' % (OS, OC)


    def Get_OX(self, id):
        entry = self.Get(id)
        if not entry: return None
        OS = ''
        OC = ''
	OX = ''
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = OS + string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
            elif line[0:5] == 'OC   ':
                OC = OC + string.strip(line[5:])
                if OC[-1] ==".": OC = OC[0:-1]
            elif line[0:5] == 'OX   ':
                OX = OX + string.strip(line[16:])
                if OX[-1] ==";": OX = OX[0:-1]
	    
	    if line[0:2] =="//": break
        if not len(OS): return id
        return (OS, OC, OX)
     
    def Get_Organism(self, id):
        entry = self.Get(id)
        if not entry: return None
        OS = ''
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = OS + string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
                return OS
            if line[0:2] =="//": break
        return OS

    def FixOC(self, OC):
        D = {

            'alpha subdivision':'Alphaproteobacteria',
            'beta subdivision':'Betaproteobacteria',
            'gamma subdivision':'Gammaproteobacteria',
            
#             'delta subdivision':'Deltaproteobacteria',
#             'epsilon subdivision':'Epsilonproteobacteria',
            }
        for old, new in D.items():
            OC = OC.replace(old,new)
        return OC
        
    def Get_TaxonomyKW(self, id):
        entry = self.Get(id)
        if not entry: return None
        OC = ""
        KW = ""
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OC   ':
                OC = OC + ' ' + string.strip(line[5:])
                if OC[-1] ==".": OC = OC[0:-1]
            if line[0:5] == 'KW   ':
                KW = KW + string.strip(line[5:])
                if KW[-1] ==".": KW = KW[0:-1]
            if line[0:2] =="//": break
        OC = self.FixOC(OC)
        return '%s %s' % (OC, KW)

    def FixOS(self, os):
        os = string.split(os,',')[0]
        os = string.split(os,'(')[0]
        return string.strip(os)
    
    def Get_Taxonomy(self, geneid):
        entry = self.Get(geneid)
        if not entry: return None
        OC = []

        #print >> sys.stderr, 'TTTT', OC, len(OC), entry

        for line in entry.splitlines():
            #print >> sys.stderr, 'OCCCCCCCCCCC <%s>' % line
            if line[0:5] == 'OC   ':
                OC.append(line[5:].strip())
            elif line[0:2] =="//": break

        OC = ' '.join(OC)
        #print >> sys.stderr, 'ttttttttttttt', geneid, OC, len(OC)
        try:
            if OC[-1] ==".": OC = OC[0:-1]
            OC = self.FixOC(OC)
        except:
            print >> sys.stderr, 'Error in',  geneid, OC, 'ENTRY <<<%s>>>' % entry
            print >> sys.stderr, self.GetValues(geneid)
            print >> sys.stderr, 'ENTRY1 <%s>' % self.temp_entries.get(geneid,'NADA')
            print >> sys.stderr, 'ENTRY2 <%s>' % self.Get(geneid)
            sys.exit(1)
            return None
        return OC
        
    def Get_Kingdom(self, geneid):
        res = self.Get_Taxonomy(geneid)
        #print id, res
        if not res: return "U"
        kd = string.strip(string.split(res,";")[0])
        if kd == "Eubacteria" or kd == "Prokaryota" or kd == "Bacteria": return "B"
        elif kd == "Eukaryota" or kd =="Eukaryotae": return "E"
        elif kd == "Archaebacteria" or kd == "Archaea": return "A"
        elif kd == "Viridae" or kd == "Viruses": return "V"
        else:
            print kd, "UNKNOWN"
            return "U"
        
    def Get_Gene(self, id):
        entry = self.Get(id)
        if not entry: return None
        GN = ''
        for line in string.split(entry, '\n'):
            if line[0:5] == 'GN   ':
                GN = string.strip(line[5:])
                if GN[-1] ==".": GN = GN[0:-1]
                return GN
            if line[0:2] =="//": break
        return GN

    def Get_EC(self, id):
        OS, GC, DE = self.Get_OS_GN_DE(id)
        if not DE: return None
        m = re.search('(EC=([0-9\.]+))',DE)
        if not m: return None
        return m.group(2)
    
    def Get_OS_OC_GN_EC(self, id):
        OS, OC, GN = self.Get_OS_OC_GN(id)
        EC = self.Get_EC(id)
        return (OS, OC, GN, EC)
    
    def Get_OS_OC_GN(self, id):
        entry = self.Get(id)
        if not entry: return None, None, None
        OS, OC, GN = "","",""
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = OS + string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
            if line[0:5] == 'OC   ':
                OC = OC + ' ' + string.strip(line[5:])
                if OC[-1] ==".": OC = OC[0:-1]
            if line[0:5] == 'GN   ':
                GN = string.strip(line[5:])
                if GN[-1] ==".": GN = GN[0:-1]
            if line[0:2] =="//": break
        OC = self.FixOC(OC)
        return OS, OC, GN
    
    def Get_OS_OC_DE(self, id):
        entry = self.Get(id)
        if not entry: return None, None, None
        OS, OC, DE = "","",""
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = OS + string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
            if line[0:5] == 'OC   ':
                OC = OC + ' ' + string.strip(line[5:])
                if OC[-1] ==".": OC = OC[0:-1]
            if line[0:5] == 'DE   ':
                DE = DE + string.strip(line[5:])
                if DE and DE[-1] ==".": DE = DE[0:-1]

            if line[0:2] =="//": break
        OC = self.FixOC(OC)
        return OS, OC, DE
    
    def Get_OS_OC_OG(self, id):
        entry = self.Get(id)
        if not entry: return None, None, None
        OS, OC, OG = "","",""
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = OS + string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
            if line[0:5] == 'OC   ':
                OC = OC + ' ' + string.strip(line[5:])
                if OC[-1] ==".": OC = OC[0:-1]
            if line[0:5] == 'OG   ':
                OG = string.strip(line[5:])
                if OG[-1] ==".": OG = OG[0:-1]
            if line[0:2] =="//": break
        OC = self.FixOC(OC)            
        return OS, OC, OG
    
    def Get_OS_GN_DE(self, id):
        entry = self.Get(id)
        if not entry: return None, None, None
        OS, GN, DE = "","",""
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = OS + string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
            if line[0:5] == 'GN   ':
                GN = string.strip(line[5:])
                if GN[-1] ==".": GN = GN[0:-1]
            if line[0:2] =="//": break
            if line[0:5] == 'DE   ':
                DE = DE + string.strip(line[5:])
                if DE and DE[-1] ==".": DE = DE[0:-1]
            if line[0:2] =="//": break
        return OS, GN, DE

    def Get_SQ(self, id, fasta = 1):
        entry = self.Get(id)
        if not entry: return ""
        SQ = ""
        record = 0
        for line in string.split(entry, '\n'):
            if record: SQ = SQ + string.strip(line[5:])
            if line[0:5] == 'SQ   ': record = 1
            if line[0:2] =="//": break

        SQ = re.sub('[ \n]','',SQ)
        if fasta: SQ = '>%s\n%s' % (id, re.sub('(.{60})','\\1\n',SQ))
        return SQ

    def Get_XX(self, id, xx):
        entry = self.Get(id)
        if not entry: return ""
        XX = ""
        for line in string.split(entry, '\n'):
            if line[0:5] == '%s   ' % xx:
                XX = XX + string.strip(line[5:])
                if XX[-1] ==".": XX = XX[0:-1]
            if line[0:2] =="//": break
        return XX
    
    def Get_KEGG(self, id):
        entry = self.Get(id)
        if not entry: return ""
        res = []
        for line in entry.splitlines():
            if line[:10] == 'DR   KEGG;':
                res.append(line.split(';')[1].strip())
            if line[0:2] =="//": break
        return res
    
    def Get_GO(self, id):
        entry = self.Get(id)
        if not entry: return ""
        res = [id]
        for line in entry.splitlines():
            if line.startswith('DR   GO;'):
                res.append(line.split(';')[1].strip())
            if line[0:2] =="//": break
        return res
        
    def Get_Keywords(self, id):
        entry = self.Get(id)
        if not entry: return []
        keywords = []
        for line in string.split(entry, '\n'):
            if line[0:5] == 'KW   ':
                for i in string.split(string.strip(line[5:]),';'):
                    kw = string.strip(i)
                    if len(kw) < 2: continue
                    if kw[-1] == '.': kw = kw[:-1]
                    keywords.append(kw)
            if line[0:2] =="//": break
        return keywords
    
    def GetAnnotations(self, swiss_id):
        entry = self.Get(swiss_id)
        if not entry: return None, None, None
        GN, KW, DE = "","",""
        for line in string.split(entry, '\n'):
            if line[0:5] == 'GN   ':
                GN = GN + string.strip(line[5:])
                if GN[-1] ==".": GN = GN[0:-1]
            if line[0:5] == '	KW   ':
                KW = string.strip(line[5:])
                if KW[-1] ==".": KW = KW[0:-1]
            if line[0:2] =="//": break
            if line[0:5] == 'DE   ':
                DE = string.strip(line[5:])
                if DE[-1] ==".": DE = DE[0:-1]
            if line[0:2] =="//": break
        return (GN, KW, DE)
        
    def GetHomologues(self, swiss_id):
        try:
            gene, taxa = swiss_id.split('_')
        except:
            return []
        
        import commands
        
        # check if we use pyphy database files in $PYPHY or somewhere else
        pyphydb = os.environ.get('PYPHYDB',None)
        if not pyphydb: pyphydb = os.environ['PYPHY']

        nrid = os.path.join(pyphydb,'pyphynr.id')
        com = "grep \"^%s_\" %s" % (gene, nrid)
        id_list = [id.strip() for id in commands.getoutput(com).split('\n')]
        return id_list


    def Name2NiceName(self, swissid):
        OS, GN, DE = self.Get_OS_GN_DE(swissid)
        if self.debug:
            print >> sys.stderr, 'debugging!!!'
            print >> sys.stderr, 'Name2NiceName OS', swissid, OS
            print >> sys.stderr, 'Name2NiceName GN', swissid, GN
            print >> sys.stderr, 'Name2NiceName ID', swissid, DE

        if not GN:
            GN = DE
        if not GN:
            return swissid
        
        if GN.find('Name=') >-1:
            gene = GN.split('Name=')[1].split(';')[0]
        elif GN.find('ORFNames') > -1:
            gene = GN.split('ORFNames=')[1].split(';')[0].split(',')[0].split()[0]
        elif GN.find('OrderedLocusNames=') > -1:
            gene = GN.split('OrderedLocusNames=')[1].split(';')[0].split(',')[0].split()[0]
        else:
            gene = re.sub('\(.+\)','', GN).strip().replace(' ','_')
            
        OS = re.sub('\(.+\)','',OS).strip().replace(' ','_')

        name = '%s_%s' % (OS, gene)
        return name

    def Translate(self, file):
        rx = re.compile('[A-Z0-9]+_[A-Z0-9]+')
        if file == '-':
            fid = sys.stdin
        else:
            fid = open(file)
        for line in fid.readlines():
            line = line.strip()
            for hit in rx.findall(line):
                taxa = self.Get_Taxonomy(hit)
                if taxa:
                    line = line + '    ' + taxa
            print line
            

        
def help(exit = 0):
    name = os.path.basename(sys.argv[0])
    print 'Usage: %s <db> <gene ID>' % name
    print '  or   %s --index <db.dat>' % name
    if exit: sys.exit(0)

if __name__ == '__main__':
    # check if we use pyphy database files in $PYPHY or somewhere else
    if len(sys.argv) == 1: help(exit = 1)
    debug = 0
    ids = None
    optlist, args = getopt.getopt(sys.argv[1:], 'hg', [ 'index=', 'index2=', 'prefix=', 'get=', 'db=', 'debug', 'file=', 'translate='])
    db_index = DB_Index(debug = debug)   
    func = db_index.Get

    for arg, value in optlist:
        if arg == '--prefix':
            prefix = value
        elif arg == '--get':
            func = getattr(db_index, value)
        elif arg == '--db':
            db = value
        elif arg in ['--debug', '-g']:
            debug = 1
        elif arg == '--file':
            if value == '-':
                fid = sys.stdin
            else:
                fid = open(value)
            ids = [x.strip() for x in fid.readlines()]
        elif arg == '--translate':
            translate_file = value
        elif arg == '-h' or arg == '--help': help(exit = 1)


    if not ids:
        ids = args
    db_index.debug = debug
    for id in ids:
        #print db_index.Get(id)
        res = func(id.strip())
        if res: print res
        


