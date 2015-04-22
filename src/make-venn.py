#!/usr/bin/python

from matplotlib_venn import venn2,venn3,venn2_circles,venn3_circles
from utilsPoirel import *
import matplotlib.pyplot as plt
import sys
import os
from optparse import OptionParser

DBS = ['netpath','kegg']

NAMES = {'EGFR1':'EGFR',
         'TCR':'T Cell Receptor',
         'TGF_beta_Receptor':'TGF Beta Receptor',
         'BCR':'B Cell Receptor'}

############################################################
## getKEGGPathways reads the pathways for KEGG and returns
## them as a list. It also returns a dictionary of {keggid:netpathname}
def getKEGGPathways():
    analyzedpathwayfile = 'data/kegg-analyzed-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,2)]
    # dictionary of keggnames to netpath names.
    kegg2netpath = readDict(analyzedpathwayfile,2,1)
    return pathways,kegg2netpath

############################################################
def readPositives(pathways,datadir):
    keggpathways,kegg2netpath = getKEGGPathways()
    netpath2kegg = {kegg2netpath[k]:k for k in kegg2netpath}

    nodes = {pathway:{db:set() for db in DBS} for pathway in pathways}
    edges = {pathway:{db:set() for db in DBS} for pathway in pathways}
    for pathway in pathways:
        ## READ FILES
        
        for db in DBS:
            if db=='netpath':
                prefix = '%s/interactions/netpath/pathways/%s' % (datadir,pathway)
            else:
                prefix = '%s/interactions/kegg/2015-03-23/hsa/edge-files/%s' % (datadir,netpath2kegg[pathway])
            nodes[pathway][db] = readItemSet('%s-nodes.txt' % (prefix),1)

            ## edges are sorted (direction doesn't matter here)
            edgelines = readColumns('%s-edges.txt' % (prefix),1,2)
            edgelist = [tuple(sorted(list(tup))) for tup in edgelines]
            edges[pathway][db] = set(edgelist)
        for db in DBS:
            print '%s %s: %d nodes' % (pathway,db,len(nodes[pathway][db]))
            print '%s %s: %d edges' % (pathway,db,len(edges[pathway][db]))

    return nodes,edges

############################################################
def readPathLinker(pathways,inputdir,topk,):
    prednodes = {pathway:set() for pathway in pathways}
    prededges = {pathway:set() for pathway in pathways}
    for pathway in pathways:
        ## READ FILES
        prefix = '%s/%s-exclude_none-sample_50X' % (inputdir,pathway)
        prednodes[pathway] = set([n for n,val,t in readColumns('%s-node-precision-recall.txt' % (prefix),1,2,3) if val != 'Inf' and float(val) <= topk and t != 'ignore'])

        ## edges are sorted (direction doesn't matter here)
        edgelines = [(u,v) for u,v,val,t in readColumns('%s-edge-precision-recall.txt' % (prefix),1,2,3,4) if val != 'Inf' and float(val) <= topk and t != 'ignore']
        edgelist = [tuple(sorted(list(tup))) for tup in edgelines]
        prededges[pathway] = set(edgelist)

    return prednodes,prededges

############################################################
## http://matthiaseisen.com/pp/patterns/p0144/ for ideas
def twowayvenn(items,itemtype,prefix,pathways):

    ## NODES
    fig = plt.figure(figsize=(10,8))
    subplot = 1
    pathways = sorted(pathways)
    for i in range(len(pathways)):
        pathway = pathways[i]
        ax = fig.add_subplot(2,3,i)
    
        venn2([items[pathway][db] for db in DBS], set_labels=('NetPath','KEGG'),ax=ax)
        c = venn2_circles(subsets=[items[pathway][db] for db in DBS], linestyle='solid')
        ax.set_title('%s %s' % (NAMES.get(pathway,pathway),itemtype),fontsize='18')
  
    figname = '%s.png' %(prefix)
    print '\twriting to %s' % (figname)
    plt.savefig(figname)
    figname = '%s.pdf' %(prefix)
    print '\twriting to %s' % (figname)
    plt.savefig(figname)
    os.system('pdfcrop %s %s' %(figname, figname))

    return 

############################################################
## http://matthiaseisen.com/pp/patterns/p0144/ for ideas
def threewayvenn(items,itemtype,prefix,pathways):
    toplot = DBS + ['pathlinker']

    ## NODES
    fig = plt.figure(figsize=(10,8))
    subplot = 1
    pathways = sorted(pathways)
    for i in range(len(pathways)):
        pathway = pathways[i]
        ax = fig.add_subplot(2,3,i)
    
        venn3([items[pathway][k] for k in toplot], set_labels=('NetPath','KEGG','PathLinker'),ax=ax,normalize_to=5.0)
        c = venn3_circles(subsets=[items[pathway][k] for k in toplot], linestyle='solid',normalize_to=5.0)
        c[2].set_ls('dashed')  # Line style
        c[0].set_lw(0.5)  # Line style
        c[1].set_lw(0.5)  # Line style
        c[2].set_lw(0.5)  # Line style
        ax.set_title('%s %s' % (NAMES.get(pathway,pathway),itemtype),fontsize='16')
  
    figname = '%s.png' %(prefix)
    print '\twriting to %s' % (figname)
    plt.savefig(figname)
    figname = '%s.pdf' %(prefix)
    print '\twriting to %s' % (figname)
    plt.savefig(figname)
    os.system('pdfcrop %s %s' %(figname, figname))

    return 

############################################################
def main(args):

    usage = '''compute-aggregate-precision-recall.py [options]
    
    Computes aggregate precision and recall from already-computed individual files.
    '''
    parser = OptionParser(usage=usage)
    
    parser.add_option('','--datadir',type='str',metavar='STR',\
                      help='Data directory from SVN. Required.')
    parser.add_option('--inputdir', type='string',metavar='STR',\
                      help='Directory PathLinker precision recal values. If specified, 3-way .')
    parser.add_option('','--topk',type='int',default=200,metavar='INT',\
                     help='# of paths to visualize (Default is 200).')
    
    # parse the command line arguments
    (opts, args) = parser.parse_args()

    if opts.datadir == None:
        sys.exit('ERROR: data directory must be specified. Exiting.')

    ## get the pathways that KEGG & NetPath have in common.
    pathways = readItemSet('data/netpath-dbcompare-pathways.txt',1)
    print '%d pathways:' % (len(pathways))
    print pathways
 
    ## get positives for each pathway
    nodes,edges = readPositives(pathways,opts.datadir)

    ## make 2-way venn:
    twowayvenn(nodes,'Proteins','viz/venn/positives-nodes',pathways)
    twowayvenn(edges,'Interactions','viz/venn/positives-edges',pathways)

    if opts.inputdir:
        pass
        ## read all PathLinker
        prednodes,prededges = readPathLinker(pathways,opts.inputdir,opts.topk)

        for pathway in pathways:
            nodes[pathway]['pathlinker'] = prednodes[pathway]
            edges[pathway]['pathlinker'] = prededges[pathway]

        ## make 3D venn
        threewayvenn(nodes,'Proteins','viz/venn/predicted-nodes',pathways)
        threewayvenn(edges,'Interactions','viz/venn/predicted-edges',pathways)
    return
        

if __name__=='__main__':
    main(sys.argv)
