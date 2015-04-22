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
def getNetPathPathways(onlykeggoverlap):
    if onlykeggoverlap:
        analyzedpathwayfile = 'data/netpath-dbcompare-pathways.txt'
    else:
        analyzedpathwayfile = 'data/netpath-analyzed-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,1)]
    return pathways

###########################################################
def getKEGGPathways(onlynetpathoverlap):
    ## TODO: update this with kegg-dbcompare-pathways.txt 
    ## now, the only KEGG pathways are the 6 that overlap with NetPath
    if onlynetpathoverlap:
        analyzedpathwayfile = 'data/kegg-dbcompare-pathways.txt'
    else:
        analyzedpathwayfile = 'data/kegg-analyzed-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,2)]
    # dictionary of keggnames to netpath names.
    kegg2netpath = readDict(analyzedpathwayfile,2,1)
    return pathways,kegg2netpath

############################################################
def readFiles(indir,param,negtype,subsample,pathways,descending,debug):
    # append pathway to tuple of columns.
    negEdges = set()
    posEdges = set()
    predEdges = set()
    negNodes = set()
    posNodes = set()
    predNodes = set()
    for pathway in pathways:
        # parse edges
        if not debug:
            if not param:
                edgefile = '%s/%s-exclude_%s-sample_%dX-edge-precision-recall.txt' % \
                           (indir,pathway,negtype,subsample)
            else:
                edgefile = '%s/%s-%s-exclude_%s-sample_%dX-edge-precision-recall.txt' % \
                           (indir,pathway,param,negtype,subsample)
        else:
            # old file (for debugging)
            # param doesn't matter here.
            edgefile = '/data/annaritz/signaling/2014-11-weighted-interactome/weighted-viz/precrecfiles-sample-once-per-pathway/varyparams-exclude_%s-KSP-0.5-%s-ranked-edges-pos_neg.txt' % (negtype,pathway)

        if not os.path.isfile(edgefile):
            sys.exit('ERROR: Edge file is not found: %s' % (edgefile))
        if not debug:
            lines = readColumns(edgefile,1,2,3,4)
        else: # old version has an extra column.
            lines = readColumns(edgefile,2,3,4,5)

        for u,v,val,etype in lines:
            if val != 'Inf': ## if it's a prediction (not Inf), add it to predEdges
                predEdges.add( ((pathway,(u,v)), float(val)) )
            if 'pos' in etype: ## if pos or positive, add it to posEdges
                posEdges.add( (pathway,(u,v)) )
            if 'neg' in etype: ## if neg or negative, add it to negEdges
                negEdges.add( (pathway,(u,v)) )
                
        # parse nodes
        if not debug:
            if not param:
                nodefile = '%s/%s-exclude_%s-sample_%dX-node-precision-recall.txt' % \
                           (indir,pathway,negtype,subsample)
            else:
                nodefile = '%s/%s-%s-exclude_%s-sample_%dX-node-precision-recall.txt' % \
                           (indir,pathway,param,negtype,subsample)
        else:
            # old file (for debugging)
            nodefile = '/data/annaritz/signaling/2014-11-weighted-interactome/weighted-viz/precrecfiles-sample-once-per-pathway/varyparams-exclude_%s-KSP-0.5-%s-ranked-nodes-pos_neg.txt' % (negtype,pathway)

        if not os.path.isfile(nodefile):
            sys.exit('ERROR: Node file is not found: %s' % (nodefile))
        if not debug:
            lines = readColumns(nodefile,1,2,3)
        else:# old version has an extra column
            lines = readColumns(nodefile,2,3,4)
        for n,val,ntype in lines:
            if val != 'Inf': ## if it's a prediction (not Inf), add it to predNodes
                predNodes.add( ((pathway,n), float(val)) )
            if 'pos' in ntype: ## if pos or positive, add it to posNodes
                posNodes.add( (pathway,n) )
            if 'neg' in ntype: ## if neg or negative, add it to negNodes
                negNodes.add( (pathway,n) )

    # sort predicted nodes and edges
    if descending:
        predEdges = sorted(predEdges, key=lambda x: x[1],reverse=True)
        predNodes = sorted(predNodes, key=lambda x: x[1],reverse=True)
    else:
        predEdges = sorted(predEdges, key=lambda x: x[1])
        predNodes = sorted(predNodes, key=lambda x: x[1])

    return posEdges,negEdges,predEdges,posNodes,negNodes,predNodes

############################################################
def main(args):

    usage = '''compute-aggregate-precision-recall.py [options]

Computes aggregate precision and recall from already-computed individual files.
'''
    parser = OptionParser(usage=usage)
    parser.add_option('--inputdir', type='string',metavar='STR',\
                      help='Directory of already-computed individual files. Aggregate files will be placed here. Required.')
    parser.add_option('--netpath',action='store_true',default=False,\
                     help='If specified, computes netpath pathways.  Written in anticipation of adding other pathway DBs.')
    parser.add_option('--kegg',action='store_true',default=False,\
                      help='If specified, computes kegg pathways.  Written in anticipation of adding other pathway DBs.')
    parser.add_option('--negtype', type='string', metavar='STR',\
                      help='Specifies which edges to exclude from the set of negatives: none, adjacent, file. Required.')
    parser.add_option('--neg-factor', type='int', default='50', metavar='INT',\
                     help='Select f*|Positives| pairs randomly as negatives, where f is the positive integer value provided with this option.  Default is 50.')
    parser.add_option('--descending',action='store_true',default=False,\
                     help='If specified, ranks values in decreasing order. Default is to rank values in increasing order.')
    parser.add_option('--param',type='string',metavar='STR',\
                      help='Specifies parameter infix of <inputdir>/<pathway>-<param>.  E.g., q_0.50 for pagerank.')
    parser.add_option('--ignorekegg',action='store_true',default=False,\
                      help='Aggregate NetPath pathways over 6 pathways shared by KEGG. Default = False.')
    parser.add_option('--ignorenetpath',action='store_true',default=False,\
                      help='Aggregate KEGG pathways over 6 pathways shared by NetPath. Default = False.')
    parser.add_option('--union',action='store_true',default=False,\
                      help='Only consider 6 pathways shared by KEGG/Netpath.')
    parser.add_option('--debug',action='store_true',default=False,\
                      help='Read from older directory of precision and recall.  Still write to input directory.')

    
    # parse the command line arguments
    (opts, args) = parser.parse_args()

    if opts.inputdir == None:
        sys.exit('\nERROR: input directory required.  Must also have both read/write permissions.')
    if opts.negtype not in ['none', 'adjacent','file']:
        sys.exit('\nERROR: --negtype must be one of "none", "adjacent", "file"')

    print '\nOPTIONS ARE', opts
    
    # get pathways
    pathways = set()
    if opts.netpath:
        pathways = getNetPathPathways(opts.ignorekegg or opts.union)
    elif opts.kegg: # only aggregate 6 KEGG pathways
        pathways,keggnames = getKEGGPathways(opts.ignorenetpath or opts.union)
        

    # read files
    posEdges,negEdges,predEdges,posNodes,negNodes,predNodes = readFiles(opts.inputdir,opts.param,opts.negtype,opts.neg_factor,pathways,opts.descending,opts.debug)
    print '%d predicted edges and %d predicted nodes' % (len(predEdges),len(predNodes))
    print '%d positive edges and %d positive nodes' % (len(posEdges),len(posNodes))
    print '%d negative edges (%.2fX positives) and %d negative nodes (%.2fX positives)'% (len(negEdges),len(negEdges)/float(len(posEdges)),len(negNodes),len(negNodes)/float(len(posNodes)))
    if len(negEdges)/float(len(posEdges)) > opts.neg_factor:
        sys.exit('Error: # of negative edges is more than %dX times the number of positives. Exiting.'   % (opts.neg_factor))
    if len(negNodes)/float(len(posNodes)) > opts.neg_factor:
        sys.exit('Error: # of negative nodes is more than %dX times the number of positives. Exiting.'   % (opts.neg_factor))    

    # compute precision and recall
    print 'Computing precision and recall for edges...'
    PREdge = computePR(posEdges, negEdges, predEdges,compressed=False)
    
    if opts.param == None:
        if opts.ignorekegg:
            outfile = '%s/aggregate-pathways_shared_with_kegg-exclude_%s-sample_%dX-edge-precision-recall.txt' % (opts.inputdir,opts.negtype,opts.neg_factor)
        elif opts.ignorenetpath:
             outfile = '%s/aggregate-pathways_shared_with_netpath-exclude_%s-sample_%dX-edge-precision-recall.txt' % (opts.inputdir,opts.negtype,opts.neg_factor)
        else:
            outfile = '%s/aggregate-exclude_%s-sample_%dX-edge-precision-recall.txt' % (opts.inputdir,opts.negtype,opts.neg_factor)
    else:
        if opts.ignorekegg:
            outfile = '%s/aggregate-pathways_shared_with_kegg-%s-exclude_%s-sample_%dX-edge-precision-recall.txt' % (opts.inputdir,opts.param,opts.negtype,opts.neg_factor)
        elif opts.ignorenetpath:
            outfile = '%s/aggregate-pathways_shared_with_netpath-%s-exclude_%s-sample_%dX-edge-precision-recall.txt' % (opts.inputdir,opts.param,opts.negtype,opts.neg_factor)
        else:
            outfile = '%s/aggregate-%s-exclude_%s-sample_%dX-edge-precision-recall.txt' % (opts.inputdir,opts.param,opts.negtype,opts.neg_factor)
    out = open(outfile,'w')
    out.write('#pathway\tnode1_sorted\tnode2_sorted\tvalue\tpos/neg/ignore\tprecision\trecall\n')
    # write ranked values
    for item,val,itemtype,prec,rec in PREdge:
        out.write('%s\t%s\t%s\t%0.5e\t%s\t%0.5e\t%0.5e\n' % (item[0],item[1][0],item[1][1],val,itemtype,prec,rec))
    # write unranked positives.
    for item in posEdges.difference(set([e for e,t in predEdges])):
        out.write('%s\t%s\t%s\tInf\tpos\t%0.5e\t%0.5e\n' % (item[0],item[1][0],item[1][1],PREdge[-1][3],PREdge[-1][4]))
    # write unranked negatives.
    for item in negEdges.difference(set([e for e,t in predEdges])):
        out.write('%s\t%s\t%s\tInf\tneg\t%0.5e\t%0.5e\n' % (item[0],item[1][0],item[1][1],PREdge[-1][3],PREdge[-1][4]))

    print 'Wrote to %s' % (outfile)

    print 'Computing precision and recall for nodes...'
    PRNode = computePR(posNodes, negNodes, predNodes,compressed=False)

    if opts.param == None:
        if opts.ignorekegg:
            outfile = '%s/aggregate-pathways_shared_with_kegg-exclude_%s-sample_%dX-node-precision-recall.txt' % (opts.inputdir,opts.negtype,opts.neg_factor)
        elif opts.ignorenetpath:
            outfile = '%s/aggregate-pathways_shared_with_netpath-exclude_%s-sample_%dX-node-precision-recall.txt' % (opts.inputdir,opts.negtype,opts.neg_factor)
        else:
            outfile = '%s/aggregate-exclude_%s-sample_%dX-node-precision-recall.txt' % (opts.inputdir,opts.negtype,opts.neg_factor)
    else:
        if opts.ignorekegg:
            outfile = '%s/aggregate-pathways_shared_with_kegg-%s-exclude_%s-sample_%dX-node-precision-recall.txt' % (opts.inputdir,opts.param,opts.negtype,opts.neg_factor)
        elif opts.ignorenetpath:
            outfile = '%s/aggregate-pathways_shared_with_netpath-%s-exclude_%s-sample_%dX-node-precision-recall.txt' % (opts.inputdir,opts.param,opts.negtype,opts.neg_factor)
        else:
            outfile = '%s/aggregate-%s-exclude_%s-sample_%dX-node-precision-recall.txt' % (opts.inputdir,opts.param,opts.negtype,opts.neg_factor)
    out = open(outfile,'w')
    out.write('#pathway\tnode\tvalue\tpos/neg/ignore\tprecision\trecall\n')
    # write ranked values
    for item,val,itemtype,prec,rec in PRNode:
        out.write('%s\t%s\t%0.5e\t%s\t%0.5e\t%0.5e\n' % (item[0],item[1],val,itemtype,prec,rec))
    # write unranked positives.
    for item in posNodes.difference(set([n for n,t in predNodes])):
        out.write('%s\t%s\tInf\tpos\t%0.5e\t%0.5e\n' % (item[0],item[1],PRNode[-1][3],PRNode[-1][4]))
    # write unranked negatives.
    for item in negNodes.difference(set([n for n,t in predNodes])):
        out.write('%s\t%s\tInf\tneg\t%0.5e\t%0.5e\n' % (item[0],item[1],PRNode[-1][3],PRNode[-1][4]))
        
    print 'Wrote to %s' % (outfile)  
    print 'DONE'

if __name__=='__main__':
    main(sys.argv)
