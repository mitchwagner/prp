#!/usr/bin/python

## Modified from /data/poirel/research/signaling-pathways/viz/precision-recall-faster/compute-all-precision-recall.py

import os
import sys
from optparse import OptionParser, OptionGroup
from utilsPoirel import *
from scipy.stats.mstats import rankdata
import matplotlib.pyplot as plt
import os.path
import random


############################################################
def getPredictions(edgefile,edgecol,nodefile,nodecol,descending):
    ## Always read these predictions from the specified files.

    ## ANNA CHANGE: this reads both edge and node predictions FROM THE SAME FILE.
    ## Old version read node predictions separately, sometimes from different files.

    #print 'Retrieving edge and node predictions...'
    if edgecol == None: 
        # read entire file as items tied with rank 1
        lines = readColumns(edgefile,1,2)
        cols = [(t,h,'1') for t,h in lines]
    else:
        cols = readColumns(edgefile,1,2,edgecol)

    epreds = {}
    npreds = {}
    for t,h,v in cols:
        # ignore if t,h,v are weird values.
        if h=='-' or t=='-':
            continue

        val = float(v)
        if descending: # make negative if descending.
            val = -val

        # Edge key is the sorted edge (thus making it undirected)
        key = tuple(sorted([t,h]))
        ## we may have seen this key before. If so, keep the smallest or largest 
        ## value according to the descending flag
        if key not in epreds or val < epreds[key]:
            epreds[key]=val

        if nodefile == None: # read nodes as we read edges.
            ## Node keys: add both nodes as keys in preds dict.
            for key in [t,h]:
                ## we may have seen this key before. If so, keep the smallest or largest 
                ## value according to the descending flag
                if key not in npreds or \
                   (val < npreds[key] and not descending) or \
                   (val > npreds[key] and descending):
                    npreds[key]=val

    # nodes haven't been read. Read from separate file.
    if nodefile != None:
        if nodecol == None: 
            # read entire file as items tied with rank 1
            lines = readItemSet(nodefile,1)
            cols = [(n,'1') for n in lines]
        else:
            cols = readColumns(nodefile,1,nodecol)

        for n,v in cols:
            # ignore if n,v are weird values.
            if n=='-':
                continue

            val = float(v)
            if descending: # make negative if descending.
                val = -val

            key = n
            ## we may have seen this key before. If so, keep the smallest or largest 
            ## value according to the descending flag
            if key not in npreds or \
               (val < npreds[key] and not descending) or \
               (val > npreds[key] and descending):
                npreds[key]=val

    ## DON"T NEED TO SORT!!! ComputePR() does this.
    #sortedEdgePreds = sorted(epreds.items(), key=lambda x: x[1])
    #sortedNodePreds = sorted(npreds.items(), key=lambda x: x[1])

    #print '\t%d edge predictions and %d node predictions read' %(len(sortedEdgePreds),len(sortedNodePreds))
    return epreds.items(),npreds.items()

############################################################
def getPosNeg(ppifile,edgefile,nodefile,outprefix,negtype,negfactor,force,ignorednodefile,ignorededgefile):

    #print 'Reading positive node set from %s' % (nodefile)
    nodelines = readColumns(nodefile,1,2) # first column name, second column type
    posNodes = set([n for n,t in nodelines])
    receptors = set([n for n,t in nodelines if t == 'receptor'])
    tfs = set([n for n,t in nodelines if t == 'tf'])
    
    #print 'Reading positive edge set from %s' % (edgefile)
    # Need to keep directionality for when we remove some edges (e.g., incoming to receptor)
    posEdges = set(readColumns(edgefile,1,2))

    #print 'Reading PPI edges and nodes'
    # Need to keep directionality for when we remove some edges (e.g., incoming to receptor)
    ppiEdges = set([e for e in readColumns(ppifile,1,2)])
    ppiNodes = set([t for t,h in ppiEdges])
    ppiNodes.update([h for t,h in ppiEdges])

    ## It is possible that some positive edges are not in the
    ## network because (1) they are not part of the largest
    ## weakly connected component or (2) they were connected
    ## to a large hub in the network that was removed.
    if len(posEdges.difference(ppiEdges)) != 0:
        #print 'WARNING! positive edges are not a subset of the PPI.'
        posEdges.intersection_update(ppiEdges)
    if len(posNodes.difference(ppiNodes)) != 0:
        #print 'WARNING! positive nodes are not a subset of the PPI.'        
        posNodes.intersection_update(ppiNodes)

    ## Remove edges from posEdges if they have a tf in the tail or a receptor in the head.
    ## These are removed from the PPI when running PathLinker, thus they cannot be recovered.
    #print '  removing edges that have a tf in the tail or a receptor in the head'
    removedEdges = set([(tail,head) for tail,head in ppiEdges if tail in tfs or head in receptors])
    #print 'posEdges before removing: %d' % (len(posEdges))
    posEdges.difference_update(removedEdges)
    #print 'posEdges after removing: %d ' % (len(posEdges))

    negEdgeOutfile = '%s-exclude_%s-%dX-negative-edges.txt' % (outprefix,negtype,negfactor)
    negNodeOutfile = '%s-exclude_%s-%dX-negative-nodes.txt' % (outprefix,negtype,negfactor)
    if force or not os.path.isfile(negEdgeOutfile) or not os.path.isfile(negNodeOutfile):
        print 'Constructing all negative edges and nodes'
        ## need to remove (u,v) and (v,u) for all (u,v) in posEdges.  
        ## note that we cannot sort posEdges and allNegEdges because
        ## we need to remove edges that have a tf in the tail or a 
        ## receptor in the head.
        allNegEdges = ppiEdges.difference(set([(u,v) for u,v in posEdges]).union(set([(v,u) for u,v in posEdges])))
        allNegNodes = ppiNodes.difference(posNodes)
        
        #print '  removing edges that have a tf in the tail or a receptor in the head'
        allNegEdges.difference_update(removedEdges)
        before = len(allNegEdges)
        if negtype == 'adjacent': # remove edges where one node is a positive.
            #print '  removing pathway-adjacent negatives'
            adjacentEdges = set([(tail,head) for tail,head in ppiEdges if tail in posNodes or head in posNodes])
            allNegEdges.difference_update(adjacentEdges)

            adjacentNodes = set([t for t,h in adjacentEdges])
            adjacentNodes.update(set([h for t,h in adjacentEdges]))
            allNegNodes.difference_update(adjacentNodes)

        elif negtype == 'file': # read edges/nodes from ignored files and remove these.
            #print '  removing negatives from ignored files'
            ignoredEdges = readColumns(ignorededgefile,1,2)
            allNegEdges.difference_update(ignoredEdges)

            ignoredNodes = readItemSet(ignorednodefile,1)
            allNegNodes.difference_update(ignoredNodes)
        
        ## Now, we can sort the positive and the negative edges.
        ## ANNA CHANGE: in the old scripts, we subsampled with directed edges, then 
        ## converted to undirected right before we set. Now we sort before subsampling.
        posEdges = set([tuple(sorted([t,h])) for t,h in posEdges])
        allNegEdges = set([tuple(sorted([t,h])) for t,h in allNegEdges])

        if len(posEdges)*negfactor > len(allNegEdges):
            negEdges = allNegEdges
        else:
            #print '  subsampling negative edge set...'
            negEdges = set(random.sample(allNegEdges, len(posEdges)*negfactor))

        if len(posNodes)*negfactor > len(allNegNodes):
            negNodes = allNegNodes
        else:
            #print '  subsampling negative node set...'
            negNodes = set(random.sample(allNegNodes, len(posNodes)*negfactor))
   
        # write output files
        out = open(negEdgeOutfile, 'w')
        out.write('#node1\tnode2\n')
        for (u,v) in negEdges:
            out.write('%s\t%s\n' % (u,v))
        out.close()
        #print 'Wrote to %s' % (negEdgeOutfile)

        out = open(negNodeOutfile, 'w')
        out.write('#node\n')
        for n in negNodes:
            out.write('%s\n' % (n))
        out.close()
        #print 'Wrote to %s' % (negNodeOutfile)

    else:
        print 'Reading negative edge set from %s' % (negEdgeOutfile)
        negEdges = set(readColumns(negEdgeOutfile,1,2))

        print 'Reading negative node set from %s' % (negNodeOutfile)
        negNodes = readItemSet(negNodeOutfile,1)

        ## have to sort positive edges here; negative edges are already sorted.
        posEdges = set([tuple(sorted([t,h])) for t,h in posEdges])

    if len(posEdges.intersection(negEdges))!=0:
        sys.exit("ERROR: there is a positive edge in the negative edge set. Exiting.")

    return posEdges,negEdges,posNodes,negNodes

def main(args):
    global VALID_METHODS

    usage = '''compute-precision-recall.py [options]

Computes precision and recall for ranked nodes and ranked edges; outputs these files to a user-specified location.  When sampling positives and negatives, first checks to see if these have already been sampled (unless --force is specified) and will read from those files..
'''
    parser = OptionParser(usage=usage)

    parser.add_option('--force', action='store_true',\
                      help='Overwrite existing files.')

    parser.add_option('--outprefix', type='string',metavar='STR',\
                      help='A string to prepend to all output files. Required.')

    group = OptionGroup(parser,'Prediction Files and Arguments')
    group.add_option('--edgefile',type='string',metavar='STR',\
                      help='File used for computing precision/recall of edges. Required.')
    group.add_option('--edgecol',type='int',metavar='INT',\
                      help='Column of edgefile to rank predictions. If not specified, takes entire file.')
    group.add_option('--nodefile',type='string',metavar='STR',\
                     help='File used for computing precision/recall of nodes. If not specified, --edgefile is used.')
    group.add_option('--nodecol',type='int',metavar='INT',\
                     help='Column of nodefile to rank predictions.  Required if --nodefile is specified.')
    group.add_option('--descending',action='store_true',default=False,\
                      help='If specified, ranks values in decreasing order. Default is to rank values in increasing order.')
    parser.add_option_group(group)

    group = OptionGroup(parser,'Ground Truth Files and Arguments')
    group.add_option('--trueedgefile',type='string',metavar='STR',\
                     help='File of positive edges. Directed edges are specified as (u,v) in columns 1 and 2. Required.')
    group.add_option('--truenodefile',type='string',metavar='STR',\
                     help='File of positive nodes in column 1. Required.')
    group.add_option('--sampledoutprefix',type='string',metavar='STR',\
                     help='Prefix of sampled output; will read from these files if they are available; will write to this directory if files don\'t exist. This directory will be empty before you run any instances of this script. Required.')
    parser.add_option_group(group)

    group = OptionGroup(parser,'Interactome and Sampling Arguments')
    group.add_option('--ppi',type='string',metavar='STR',\
                     help='PPI to sample negatives from.  Predictions must be a subset of the edges/nodes from this PPI.  Required.')
    group.add_option('--negtype', type='string', metavar='STR',\
                      help='Specifies which edges to exclude from the set of negatives: none, adjacent, file. Required.')
    group.add_option('--ignorededgefile',type='string',metavar='STR',\
                     help='File of edges to ignore. Directed edges are specified as (uv) in columns 1 and 2. Required if "--negtype file" is specified.')
    group.add_option('--ignorednodefile',type='string',metavar='STR',\
                     help='File of nodes to ignore in column 1. Required if "--negtype file" is specified.')
    group.add_option('--neg-factor', type='int', default='10', metavar='INT',\
                     help='Select f*|Positives| pairs randomly as negatives, where f is the positive integer value provided with this option.  Default is 10.')
    parser.add_option('','--randseed',type='int',metavar='INT',default=1234567,\
                      help='Seed for random number generator Optional; default = 1234567.')
    parser.add_option_group(group)
    
    # parse the command line arguments
    (opts, args) = parser.parse_args()

    if opts.edgefile == None:
        sys.exit('\nERROR: edgefile required.')
    if opts.edgecol == None:
        print 'WARNING: No edge column specified.  Taking entire file as single point.'
    if opts.nodefile != None and opts.nodecol == None:
        print 'WARNING: Node file specified but no node column specified.  Taking entire file as single point.'
    if opts.trueedgefile == None or opts.truenodefile == None or opts.sampledoutprefix == None:
        sys.exit('\nERROR: true edge file, true node file, and sampled output prefix name are required.')
    if opts.ppi == None:
        sys.exit('\nERROR: PPI must be specified.')
    if opts.outprefix == None:
        sys.exit('\nERROR: output file prefix is required.')
    if opts.negtype not in ['none', 'adjacent','file']:
        sys.exit('\nERROR: --negtype must be one of "none", "adjacent", or "file"')
    if opts.negtype == 'file' and (opts.ignorednodefile==None or opts.ignorededgefile==None):
        sys.exit('\nERROR: if --negtype is "file", then --ignorednodefile and --ignorededgefile must be specified.')

    # seed the random number generator
    random.seed(opts.randseed)
    
    print '\nOPTIONS ARE',opts
    
    # Get the pathway predictions
    ## these contain: tuple ((u,v), value), sorted in ascending or descending order.
    predEdges,predNodes = getPredictions(opts.edgefile,opts.edgecol,opts.nodefile,opts.nodecol,opts.descending)
    print '%d predicted edges and %d predicted nodes' % (len(predEdges),len(predNodes))

    # Get the set of positives and negatives.
    posEdges,negEdges,posNodes,negNodes = getPosNeg(opts.ppi,opts.trueedgefile,opts.truenodefile,opts.sampledoutprefix,opts.negtype,opts.neg_factor,opts.force,opts.ignorednodefile,opts.ignorededgefile)
    print '%d positive edges and %d positive nodes' % (len(posEdges),len(posNodes))
    print '%d negative edges (%.2fX positives) and %d negative nodes (%.2fX positives)'% (len(negEdges),len(negEdges)/float(len(posEdges)),len(negNodes),len(negNodes)/float(len(posNodes)))
    if len(negEdges)/float(len(posEdges)) > opts.neg_factor:
        sys.exit('ERROR: # of negative edges is more than %dX times the number of positives. Exiting.'   % (opts.neg_factor))
    if len(negNodes)/float(len(posNodes)) > opts.neg_factor:
        sys.exit('ERROR: # of negative nodes is more than %dX times the number of positives. Exiting.'   % (opts.neg_factor))

    # compute precision and recall
    print 'Computing precision and recall for edges...'
    PREdge = computePR(posEdges, negEdges, predEdges,compressed=False)
    
    outfile = '%s-exclude_%s-sample_%dX-edge-precision-recall.txt' % (opts.outprefix,opts.negtype,opts.neg_factor)
    out = open(outfile,'w')
    out.write('#node1_sorted\tnode2_sorted\tvalue\tpos/neg/ignore\tprecision\trecall\n')
    # write ranked values
    for item,val,itemtype,prec,rec in PREdge:
        out.write('%s\t%s\t%0.5e\t%s\t%0.5e\t%0.5e\n' % (item[0],item[1],val,itemtype,prec,rec))
 
    if len(PREdge)==0:   
        # write unranked positives.
        for item in posEdges.difference(set([e for e,t in predEdges])):
            out.write('%s\t%s\tInf\tpos\t%0.5e\t%0.5e\n' % (item[0],item[1],0,0))
        # write unranked negatives.
        for item in negEdges.difference(set([e for e,t in predEdges])):
            out.write('%s\t%s\tInf\tneg\t%0.5e\t%0.5e\n' % (item[0],item[1],0,0))
    else:
        # write unranked positives.
        for item in posEdges.difference(set([e for e,t in predEdges])):
            out.write('%s\t%s\tInf\tpos\t%0.5e\t%0.5e\n' % (item[0],item[1],PREdge[-1][3],PREdge[-1][4]))
        # write unranked negatives.
        for item in negEdges.difference(set([e for e,t in predEdges])):
            out.write('%s\t%s\tInf\tneg\t%0.5e\t%0.5e\n' % (item[0],item[1],PREdge[-1][3],PREdge[-1][4]))
    print 'Wrote to %s' % (outfile)

    print 'Computing precision and recall for nodes...'
    PRNode = computePR(posNodes, negNodes, predNodes,compressed=False)

    outfile = '%s-exclude_%s-sample_%dX-node-precision-recall.txt' % (opts.outprefix,opts.negtype,opts.neg_factor)
    out = open(outfile,'w')
    out.write('#node\tvalue\tpos/neg/ignore\tprecision\trecall\n')
    # write ranked values
    for item,val,itemtype,prec,rec in PRNode:
        out.write('%s\t%0.5e\t%s\t%0.5e\t%0.5e\n' % (item,val,itemtype,prec,rec))
    if len(PRNode)==0:
        # write unranked positives.
        for item in posNodes.difference(set([n for n,t in predNodes])):
            out.write('%s\tInf\tpos\t%0.5e\t%0.5e\n' % (item,0,0))
        # write unranked negatives.
        for item in negNodes.difference(set([n for n,t in predNodes])):
            out.write('%s\tInf\tneg\t%0.5e\t%0.5e\n' % (item,0,0))
    else:
        # write unranked positives.
        for item in posNodes.difference(set([n for n,t in predNodes])):
            out.write('%s\tInf\tpos\t%0.5e\t%0.5e\n' % (item,PRNode[-1][3],PRNode[-1][4]))
        # write unranked negatives.
        for item in negNodes.difference(set([n for n,t in predNodes])):
            out.write('%s\tInf\tneg\t%0.5e\t%0.5e\n' % (item,PRNode[-1][3],PRNode[-1][4]))
        
    print 'Wrote to %s' % (outfile)  

    print 'DONE'

if __name__=='__main__':
    main(sys.argv)
