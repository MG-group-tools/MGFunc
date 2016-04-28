!/usr/bin/env python
# This version is created at Mon Mar 17 12:54:44 CET 2014
# Author: Asli I. Ozen (asli@cbs.dtu.dk)
# License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)
import sys, gzip
import re, string
import argparse
import os
from Bio.Blast import NCBIStandalone
from operator import itemgetter, attrgetter
import pickle
import shutil
sys.path.append('/home/projects5/pr_53035/people/asli/bin/lib/python2.7/site-packages')
import networkx as nx
from os.path import basename
import time
from datetime import datetime as dt
import matplotlib.pyplot as plt

usage = sys.argv[0] + " [-h] -p P [-b B] [-a] [-o] [-m M] [-v]\n" 
helpstr = '''
description: This script clusters genes from a given blastdecide file using reciprocal best hits or triple-hits graphs.
'''

epi="Author: Asli I. Ozen (asli@cbs.dtu.dk)" 
