#!/usr/bin/env python
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen (asli@cbs.dtu.dk)
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)


from _version import __version__

import sys, gzip
import re, string
import argparse
import os
from operator import itemgetter, attrgetter
import pickle
import networkx as nx
from os.path import basename
import commands, tempfile
import time
from datetime import datetime as dt
import fasta2index
import traceback
import glob
import Queue
import threading

from mgfunc import MGFunc

if __name__ == "__main__":
    try:
        myobj = MGFunc()
        myobj.opts = myobj.parser.parse_args(sys.argv[1:])
        myobj.printer("\n### "+sys.argv[0]+" initialized at "+ myobj.timestarted + "\n")
        myobj.printer("### OPTIONS: "+str(myobj.opts)+"\n")
        myobj.mainthing()
    #except IOError as i:
     #   print "I/O error({0}): {1}".format(i.errno, i.strerror)
    except Exception,e:
        print str(e)
        traceback.print_exc()

