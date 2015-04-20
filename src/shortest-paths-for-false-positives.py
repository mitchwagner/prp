#!/usr/bin/python

## Modified from /data/poirel/research/signaling-pathways/hop-analysis/hop-output/compute-NetPath-hops.py
from utilsPoirel import *
import sys
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from random import choice
from math import log
from optparse import OptionParser

############################################################
def neglog(w):
    return -log(max([0.000000001, float(w)]))/log(10)

############################################################
def main(args):
    usage = 'shortest-paths-for-false-positives.py [options]'
    parser = OptionParser(usage=usage)

    parser.add_option('','--ppidir',type='string',metavar='STR',\
                      help='Directory of pathway-specific interactomes. Third column includes weights between 0 and 1.  Required.')
    parser.add_option('','--datadir',type='string',metavar='STR',\
                      help='Prefix of curated node files. Required.')
    parser.add_option('','--pathway',type='string',metavar='STR',\
                      help='Pathway.')
    (opts, args) = parser.parse_args()
    if not opts.ppidir or not opts.datadir or not opts.pathway:
        parser.print_help()
        sys.exit('ERROR: argument is missing. Exiting.')

    ## get interactome nodes/edges
    ppifile = '%s/%s-interactome.txt' % (opts.ppidir,opts.pathway)
    lines = readColumns(ppifile,1,2,3)
    nodes = set([u for u,v,w in lines]).union([v for u,v,w in lines])
    ## take neg log (base 10) of weights. Add it to attr dictionary keyed by 'weight'.
    edges = [(u,v,{'weight':float(w),'neglogweight':neglog(float(w))}) for u,v,w in lines]

    ## make networkx object
    G=nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    print 'Background interactome has %d nodes and %d edges' % (G.number_of_nodes(),G.number_of_edges())

    ## get pathway nodes
    nodefile = '%s/%s-nodes.txt' % (opts.datadir,opts.pathway)
    canonnodes = readItemSet(nodefile,1)
    edgefile = '%s/%s-edges.txt' % (opts.datadir,opts.pathway)
    canonedges = readColumns(edgefile,1,2)
    
    # contract nodes in pathway to single source node.
    pathwaynode = 'PATHWAYNODE'
    G.add_node(pathwaynode)
    for n in canonnodes:
        # contract out-edges from n
        outedges = G.out_edges(n)
        for e in outedges:
            G.add_edge(pathwaynode,e[1],G.get_edge_data(e[0],e[1]))
            G.remove_edge(e[0],e[1])
        # contract in-edges to n
        inedges = G.in_edges(n)
        for e in inedges:
            G.add_edge(e[0],pathwaynode,G.get_edge_data(e[0],e[1]))
            G.remove_edge(e[0],e[1])
        # remove n
        if n in G.nodes():
            G.remove_node(n)
    print '  Contracted interactome has %d nodes and %d edges' % (G.number_of_nodes(),G.number_of_edges())

    print '  Computing DIRECTED shortest-paths...'
    paths = nx.shortest_path(G,source='PATHWAYNODE')
    print '    %d paths computed (max path length = %d)' % (len(paths),max([getpathlength(paths[n]) for n in paths]))
    print '  Computing min-cost paths...'
    paths = nx.shortest_path(G,source='PATHWAYNODE',weight='neglogweight')
    print '    %d paths computed (max path length = %d)' % (len(paths),max([getpathlength(paths[n]) for n in paths]))
    outfile = 'data/shortest-paths-for-false-positives/%s-node-dists-directed.txt' % (opts.pathway)
    out = open(outfile,'w')
    out.write('#node\tdist(edges)\tpath_weight\tneglog_path_weight\n')
    for n in canonnodes: # True Positive: 0 edges away (min) and pathweight of 1 (max)
        out.write('%s\t%d\t%d\t%d\n' % (n,0,1,neglog(1)))
    for n in nodes:
        if n in paths: # shortest path exists
            pathweight = getpathweight(paths[n],G)
            out.write('%s\t%d\t%.2e\t%.2e\n' % (n,getpathlength(paths[n]),pathweight,neglog(pathweight)))
        else: # not connected
            out.write('%s\t%d\t%d\t%d\n' % (n,-1,-1,-1))
    out.close()
    print '  wrote to %s' % (outfile)

    outfile = 'data/shortest-paths-for-false-positives/%s-edge-dists-directed.txt' % (opts.pathway)
    out = open(outfile,'w')
    out.write('#tail\thead\tdist(edges)\tpath_weight\tneglog_path_weight\n')
    for u,v,attrdict in edges: # do this for ALL edges, not just ones that are left after contraction.
        if (u,v) in canonedges: # True Positive, 0 edges away (min) and pathwweight of 1 (max)
            out.write('%s\t%s\t%d\t%d\t%d\n' % (u,v,0,1,neglog(1)))
        else: # compute shortest path to u or v + 1
            if u in canonnodes: 
                ## the tail is in the pathway; the head is not.
                ## Simply add the weight of edge (u,v).
                pathdist = 1
                pathweight = attrdict['weight']
                out.write('%s\t%s\t%d\t%.2e\t%.2e\n' % (u,v,pathdist,pathweight,neglog(pathweight)))
            elif u in paths: # shortest path to u exists: add (u,v) as last edge.
                ## note that this is the distance even if the head (v) is in the pathway.
                pathdist = getpathlength(paths[u])+1
                pathweight = getpathweight(paths[u],G)*attrdict['weight']
                out.write('%s\t%s\t%d\t%.2e\t%.2e\n' % (u,v,pathdist,pathweight,neglog(pathweight)))
            else: # path doesn't exist to (u,v) ending at v. 
                out.write('%s\t%s\t%d\t%d\t%d\n' % (u,v,-1,-1,-1))
    out.close()
    print '  wrote to %s' % (outfile)

    G_undirected = G.to_undirected()

    print '  Computing UNDIRECTED shortest-paths...'
    paths = nx.shortest_path(G_undirected,source='PATHWAYNODE')
    print '    %d paths computed (max path length = %d)' % (len(paths),max([getpathlength(paths[n]) for n in paths]))
    print '  Computing min-cost paths...'
    paths = nx.shortest_path(G_undirected,source='PATHWAYNODE',weight='neglogweight')
    print '    %d paths computed (max path length = %d)' % (len(paths),max([getpathlength(paths[n]) for n in paths]))
    outfile = 'data/shortest-paths-for-false-positives/%s-node-dists-undirected.txt' % (opts.pathway)
    out = open(outfile,'w')
    out.write('#node\tdist(edges)\tpath_weight\tneglog_path_weight\n')
    for n in canonnodes: # True Positive: 0 edges away (min) and pathweight of 1 (max)
        out.write('%s\t%d\t%d\t%d\n' % (n,0,1,neglog(1)))
    for n in nodes:
        if n in paths: # shortest path exists
            pathdist = getpathlength(paths[n])
            pathweight = getpathweight(paths[n],G_undirected)
            out.write('%s\t%d\t%.2e\t%.2e\n' % (n,pathdist,pathweight,neglog(pathweight)))
        else: # not connected
            out.write('%s\t%d\t%d\t%d\n' % (n,-1,-1,-1))
    out.close()
    print '  wrote to %s' % (outfile)

    outfile = 'data/shortest-paths-for-false-positives/%s-edge-dists-undirected.txt' % (opts.pathway)
    out = open(outfile,'w')
    out.write('#tail\thead\tdist(edges)\tpath_weight\tneglog_path_weight\n')
    seen = set()
    for u,v,attrdict in edges: # do this for ALL edges, not just ones that are left after contraction.
        sortede = tuple(sorted([u,v]))
        u = sortede[0]
        v = sortede[1]
        if sortede in seen:
            continue
        seen.add(sortede)
        if (u,v) in canonedges or (v,u) in canonedges: 
            # True Positive, 0 edges away (min) and pathwweight of 1 (max)
            out.write('%s\t%s\t%d\t%d\t%d\n' % (u,v,0,1,neglog(1)))
            out.write('%s\t%s\t%d\t%d\t%d\n' % (v,u,0,1,neglog(1)))
        else: # compute shortest path to u or v + 1
            if u in canonnodes or v in canonnodes: 
                ## if either the tail or the head is in canonnodes, add it.
                ## Simply add the weight of edge (u,v).
                pathdist = 1
                pathweight = attrdict['weight']
                out.write('%s\t%s\t%d\t%.2e\t%.2e\n' % (u,v,pathdist,pathweight,neglog(pathweight)))
                out.write('%s\t%s\t%d\t%.2e\t%.2e\n' % (v,u,pathdist,pathweight,neglog(pathweight)))
            elif u in paths and v in paths:
                ## shortest path to u and to v exists.  Take min of these and add edge.
                pathdist = min(getpathlength(paths[u]),getpathlength(paths[v]))+1
                pathweight = min(getpathweight(paths[u],G_undirected),
                                 getpathweight(paths[v],G_undirected))*attrdict['weight']
                out.write('%s\t%s\t%d\t%.2e\t%.2e\n' % (u,v,pathdist,pathweight,neglog(pathweight)))
                out.write('%s\t%s\t%d\t%.2e\t%.2e\n' % (v,u,pathdist,pathweight,neglog(pathweight)))
            elif u in paths: # shortest path to u exists. Use that.
                pathdist = getpathlength(paths[u])+1
                pathweight = getpathweight(paths[u],G_undirected)*attrdict['weight']
                out.write('%s\t%s\t%d\t%.2e\t%.2e\n' % (u,v,pathdist,pathweight,neglog(pathweight)))
                out.write('%s\t%s\t%d\t%.2e\t%.2e\n' % (v,u,pathdist,pathweight,neglog(pathweight)))
            elif u in paths: # shortest path to v exists. Use that.
                pathdist = getpathlength(paths[v])+1
                pathweight = getpathweight(paths[v],G_undirected)*attrdict['weight']
                out.write('%s\t%s\t%d\t%.2e\t%.2e\n' % (u,v,pathdist,pathweight,neglog(pathweight)))
                out.write('%s\t%s\t%d\t%.2e\t%.2e\n' % (v,u,pathdist,pathweight,neglog(pathweight)))
            else: # path doesn't exist to u or v
                out.write('%s\t%s\t%d\t%d\t%d\n' % (u,v,-1,-1,-1))
                out.write('%s\t%s\t%d\t%d\t%d\n' % (v,u,-1,-1,-1))
    out.close()
    print '  wrote to %s' % (outfile)
    return

############################################################
def getpathweight(path,G):
    pathweight = 1
    for i in range(len(path)-1):
        #print (path[i],path[i+1],G[path[i]][path[i+1]])
        pathweight*=G[path[i]][path[i+1]]['weight']
    return pathweight

############################################################
def getpathlength(path):
    return len(path)-1
        
############################################################
if __name__=='__main__':
    main(sys.argv)
