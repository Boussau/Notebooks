#! /usr/bin/env python
# Computes a consensus chronogram and saves it as a NHX file

import sys
from ete3 import Tree, TreeStyle, NodeStyle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import re

if (len(sys.argv)<3):
    print("Usage: python computeConsensusChronogram.py treefile output.nhx")
    exit(-1)

print (sys.argv)

argv=sys.argv[1:]
file=argv[0]
out=argv[1]

try:
    f=open(file, 'r')
except IOError:
    print ("Unknown file: ",file)
    sys.exit()


#File where I store useful functions
exec (open("/home/boussau/Programming/Notebooks/code/functions.py").read ())

allTrees = list()
for l in f:
    l2 = l.strip()
    # removing anything within square brackets
    if "[" in l2:
        l2 = re.sub('\[[^\]]+\]\s*','', l2)
    print(l2)
    t = Tree( l2 )
    createNameToLeavesLink( t )
    allTrees.append( t )
f.close()


id2Heights=list()
for t in allTrees:
    node2Height,id2Height = getNodeHeights( t )
    print(id2Height)
    id2Heights.append(id2Height)

#print(len(id2Heights))

# Creating a uniform weight vector
weights = [1] * len(id2Heights)
outputsWeightedChronogram (allTrees[0].copy(), id2Heights, out, weights)
