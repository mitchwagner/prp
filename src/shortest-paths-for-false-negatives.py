#!/usr/bin/python

## Modified from /data/poirel/research/signaling-pathways/hop-analysis/hop-output/compute-NetPath-hops.py

from utilsPoirel import *
import sys
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from random import choice
import os.path

## TODO replace these with arguments.
NET = '/data/annaritz/datasets/svn-data/interactomes/human/pathlinker-signaling-children-reg-weighted.txt'
CANONDIR = '/data/annaritz/datasets/svn-data/interactions/netpath/pathways/%s-nodes.txt'

# Read the pathways to plot
PATHWAYS = readItemList('data/netpath-analyzed-pathways.txt')

edges = readColumns(NET,1,2)
nodes = set([u for u,v in edges]).union([v for u,v in edges])

G=nx.DiGraph()
G.add_nodes_from(nodes)
G.add_edges_from(edges)
print 'Background interactome has %d nodes and %d edges' % (G.number_of_nodes(),G.number_of_edges())

for p in PATHWAYS:
    print p

    ## get pathway nodes
    canonnodes = [i for i in readItemSet(CANONDIR % p,1) if '-' not in i]
    print '%s nodes' % (len(canonnodes))
    outfile = 'data/shortest-paths-for-false-negatives/%s-dist.txt' % (p)
    if os.path.isfile(outfile):
        continue
    out = open(outfile,'w')
    out.write('#node1\tnode2\tshortest-path-length\n')
    for i in range(len(canonnodes)):
        print i,'...'
        n1 = canonnodes[i]
        #print '  Computing shortest-paths...'
        paths = nx.shortest_path(G.to_undirected(),source=n1)
        for j in range(i+1,len(canonnodes)):
            n2 = canonnodes[j]
            out.write('%s\t%s\t%d\n' %(n1,n2,len(paths[n2])-1))
    out.close()
    print '  wrote to %s' % (outfile)
    
print 'done'


