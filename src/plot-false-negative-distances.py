#!/usr/bin/python

## Modified from /data/poirel/research/signaling-pathways/hop-analysis/plot-recall-value.py

import os
import sys
from optparse import OptionParser, OptionGroup
from utilsPoirel import *
from scipy.stats.mstats import rankdata
import matplotlib 
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import os.path
import numpy as np
import random        

BULLSEYE_COLORS = [[0,0,0],[.2,0,0],[.4,0,0],[.6,0,0],[.8,0,0],[1,0,0]]

###############################################################
def plotFalseNegatives(rankedlists,hopdict,recalls,pdf,outprefix,edges=False,aggregate=False):
    maxhops = 6

    # write output to files.
    if edges:
        out = open(outprefix+'-edges-counts.txt','w')
        outnorm = open(outprefix+'-edges-counts-normalized.txt','w')
    else:
        out = open(outprefix+'-nodes-counts.txt','w')
        outnorm = open(outprefix+'-nodes-counts-normalized.txt','w')
    out.write('#Recall\tSubset')
    for i in range(maxhops+1):
        out.write('\tNumber-Dist%d' % (i))
    out.write('\n')

    outnorm.write('#Recall\tSubset')
    for i in range(maxhops+1):
        outnorm.write('\tProportion-Dist%d' % (i))
    outnorm.write('\n')

    subplot = 1
    barvals = {r:[] for r in recalls}
    maxval = 0
    if not edges:
        outfile = '%snodes' % (outprefix)
        title ='Predicted Nodes' 
    else:
        outfile = '%sedges' % (outprefix)
        title ='Predicted Edges' 
    for recall in sorted(recalls):
        fig = plt.figure(figsize=(8,8))
        for m in rankedlists:
            

        # set of positives
        allpositives = set([r for r in rankedlists[m] if r[-1]=='positive'])
        numpositives = float(len(allpositives))

        
            counts = [0]*(maxhops+1)
            # start with entire list:
            newlist = sorted([r for r in rankedlists[m] if r[-2] != 'Inf' and r[-1] != 'ignored'],key=lambda x: float(x[-2]))
            truepositives = set()

            # get set just from recall.
            if not edges:  ### NODE PREDICTIONS
                # get set of nodes at fixed recall
                j = 0
                while j < len(newlist) and  len(truepositives)/numpositives < recall:
                    if newlist[j][-1] == 'positive':
                        truepositives.add(newlist[j])
                    j+=1

            else:  ### EDGE PREDICTIONS
                # get set of edges at fixed recall
                j = 0
                while j<len(newlist) and len(truepositives)/numpositives < recall:
                    if newlist[j][-1] == 'positive':
                        truepositives.add(newlist[j])
                    j+=1
                    
            falsenegatives = allpositives.difference(truepositives)
            print 'recall=%.2f: %d positives: %d true positives and %d false negatives' % (recall,len(allpositives),len(falsenegatives),len(truepositives))

            # get set just from recall.
            if not edges:  ### NODE PREDICTIONS
                # get hops to each node
                for line in falsenegatives:
                    
                    if aggregate:
                        minval = maxhops
                        for p in [a for a in truepositives if a[0] == line[0]]:
                            minval = min(minval,hopdict[line[0]][line[1]][p[1]])
                        counts[minval]+=1
                    else:
                        minval = maxhops
                        for p in truepositives:
                            minval = min(minval,hopdict[line[1]][p[1]])
                        counts[minval]+=1
            else:
                for line in falsenegatives:
                    #print line
                    if aggregate:
                        n1minval = maxhops
                        n2minval = maxhops
                        for p in [a for a in truepositives if a[0] == line[0]]:
                            #print line,p
                            for j in [1,2]:
                                if line[1]==p[j]:
                                    n1minval = 0
                                else:
                                    n1minval = min(n1minval,hopdict[line[0]][line[1]][p[j]])
                                if line[2] ==p[j]:
                                    n2minval = 0
                                else:
                                    n2minval = min(n2minval,hopdict[line[0]][line[2]][p[j]])
                        num = min(n1minval,n2minval)+1
                    else:
                        n1minval = maxhops
                        n2minval = maxhops
                        for p in truepositives:
                            #print line,p
                            for j in [1,2]:
                                if line[1]==p[j]:
                                    n1minval = 0
                                else:
                                    n1minval = min(n1minval,hopdict[line[1]][p[j]])
                                if line[2] ==p[j]:
                                    n2minval = 0
                                else:
                                    n2minval = min(n2minval,hopdict[line[2]][p[j]])
                        num = min(n1minval,n2minval)+1
                    counts[num]+=1
            
            #counts[0] = len(truepositives)
            maxval=max([i for i in range(len(counts)) if counts[i]>0]+[maxval])
            barvals[recall] = counts
        print barvals
        print maxval
        counts = barvals[recall] # HACK to fit with plot-bullseye code pasted below.

        out.write('%f\t%s' % (recall,m))
        for i in range(maxval+1):
            out.write('\t%d' % (counts[i]))
        out.write('\n')

        #print counts
        # normalize and make cumulative counts
        for c in range(len(counts)):
            counts[c] = counts[c]/float(len(falsenegatives))
            if c>0:
                counts[c]+=counts[c-1]
        outnorm.write('%f\t%s' % (recall,m))
        for i in range(maxval+1):
            if i == 0:
                outnorm.write('\t%f' % (counts[i]))
            else:
                outnorm.write('\t%f' % (counts[i]-counts[i-1]))
        outnorm.write('\n')

        ## Plot Bar Chart         
        ax = fig.add_subplot(1,1,1)
        for l in reversed(range(maxval+1)):
            if l == 0:
                ax.bar(0, counts[l], width=2 * np.pi, bottom=0,
                       color=BULLSEYE_COLORS[l], edgecolor = BULLSEYE_COLORS[l],label='%d (%.2f)' % (l,counts[l]))
            elif abs(counts[l] - counts[l-1]) > 0.001:
                ax.bar(0, counts[l], width=2 * np.pi, bottom=0,
                       color=BULLSEYE_COLORS[l], edgecolor = BULLSEYE_COLORS[l],label='%d (%.2f)' % (l,counts[l]-counts[l-1]))
            else:
                ax.bar(0, counts[l], width=2 * np.pi, bottom=0,
                       color=BULLSEYE_COLORS[l], edgecolor = BULLSEYE_COLORS[l],label='_nolegend_')
        plt.ylim(0,1)
        #ax.legend(title='SPL',loc='lower right',bbox_to_anchor= (.4,-.1), prop={'size':10},)
        ax.set_yticks([])
        ax.set_xticks([])
        #ax.set_title(title)

        plt.tight_layout()
        figname = '%s.png' %(outfile)
        print '\twriting to %s' % (figname)
        plt.savefig(figname)
        if pdf:
            figname = '%s.pdf' %(outfile)
            print '\twriting to %s' % (figname)
            plt.savefig(figname)
            os.system('pdfcrop %s %s' %(figname, figname))

    out.close()
    outnorm.close()
    return

###############################################################
def readRankedLists(outputPrefix): 
    EdgeRankedList = dict([(p,{}) for p in PATHWAYS])
    NodeRankedList = dict([(p,{}) for p in PATHWAYS])
    for m in METHODS:
        for p in PATHWAYS:
            if p == 'aggregate':
                EdgeRankedList[p][m] = []
                NodeRankedList[p][m] = []
                for p2 in VALID_PATHWAYS:
                    edgeFile = '%s-%s-%s-ranked-edges-pos_neg.txt' % (outputPrefix,m,p2)
                    nodeFile = '%s-%s-%s-ranked-nodes-pos_neg.txt' % (outputPrefix,m,p2)
                    EdgeRankedList[p][m] += readColumns(edgeFile,1,2,3,4,5) #p,u,v,val,posneg
                    NodeRankedList[p][m] += readColumns(nodeFile,1,2,3,4) #p,n,val,posneg
                ## sort!
                EdgeRankedList[p][m] = sorted(EdgeRankedList[p][m],key=lambda x: float(x[3]))
                NodeRankedList[p][m] = sorted(NodeRankedList[p][m],key=lambda x: float(x[2]))
            else:
                edgeFile = '%s-%s-%s-ranked-edges-pos_neg.txt' % (outputPrefix,m,p)
                nodeFile = '%s-%s-%s-ranked-nodes-pos_neg.txt' % (outputPrefix,m,p)                
                EdgeRankedList[p][m] = readColumns(edgeFile,1,2,3,4,5) # (p,u,v,val,posneg)
                NodeRankedList[p][m] = readColumns(nodeFile,1,2,3,4) # (p,n,val,posneg)
                
    return  EdgeRankedList, NodeRankedList

###############################################################
def readHops():
    hopdict = {p:{} for p in PATHWAYS}
    for p in PATHWAYS:
        if p == 'aggregate':
            hopdict[p] = {p2:{} for p2 in VALID_PATHWAYS}
            for p2 in VALID_PATHWAYS:
                hopfile = 'data/shortest-paths-for-false-negatives/%s-dist.txt' % (p2)
                lines = readColumns(hopfile,1,2,3)
                nodes = set([u for u,v,val in lines]).union(set([v for u,v,val in lines]))
                hopdict[p][p2] = {n:{} for n in nodes}
                for u,v,val in lines:
                    hopdict[p][p2][u][v] = int(val)
                    hopdict[p][p2][v][u] = int(val)
                
        else:
            hopfile = 'data/shortest-paths-for-false-negatives/%s-dist.txt' % (p)
            lines = readColumns(hopfile,1,2,3)
            nodes = set([u for u,v,val in lines]).union(set([v for u,v,val in lines]))
            hopdict[p] = {n:{} for n in nodes}
            for u,v,val in lines:
                hopdict[p][u][v] = int(val)
                hopdict[p][v][u] = int(val)
    return hopdict

###############################################################
def main(args):
    global METHODS
    global PATHWAYS

    usage = '''prec-rec.py [options]
'''

    parser = OptionParser(usage=usage)

    # General Options

    parser.add_option('-o', '--out', type='string', default='bullseye-output/', metavar='STR',\
        help='A string to prepend to all output files. Default is "bullseye-output/"')

    parser.add_option('-i', '--prefix', type='string', default='out', metavar='STR',\
                          help='Prefix of outfile from compute-all-precision-recall.py. Ex: "all-exclude_none" ')
    
    parser.add_option('-p', '--pathway', action='append', type='str', default=[],\
                          help='A pathway to plot; this option can be provided multiple times. "aggregate" is also allowed')
    
    parser.add_option('-a', '--alg', action='append', type='str', default=[],\
                          help='An algorithm to plot; this option can be provided multiple times.')
    
    parser.add_option('', '--pdf', action='store_true', default=False,\
        help='Output images in PDF as well as PNG.')

    parser.add_option('-r','--recall',action='append',type='float',default=[],metavar='FLOAT',\
                          help='Plot bullseye plots at recall; this option can be provided multiple times')

    parser.add_option('','--edges',action='store_true',default=False,\
                          help='Plot edge figures')

    parser.add_option('','--nodes',action='store_true',default=False,\
                          help='Plot node figures')

   
    # parse the command line arguments
    (opts, args) = parser.parse_args()

    if opts.edges == False and opts.nodes == False:
        sys.exit('\nERROR: need to specify --edges, --nodes, or both. Exiting')

    if len(opts.alg)==0:
        sys.exit('\nERROR: no algorithms were provided (--alg).')

    for a in opts.alg:
        if a not in VALID_METHODS:
            sys.exit('\nERROR: "%s" is not a valid method (given to --alg).' %(a))
        METHODS.append(a)
    if len(METHODS) > 5:
        sys.exit('\nERROR: there is only room to plot 5 bulls-eye plots for 5 methods (--alg).  %d methods were passed.\n' % (len(METHODS)))

    if opts.pathway == ['all']:
        PATHWAYS = VALID_PATHWAYS
    else:
        for p in opts.pathway:
            if p not in VALID_PATHWAYS and p != 'aggregate':
                sys.exit('\nERROR: "%s" is not a valid pathway (given to --pathway).' %(p))
            PATHWAYS.append(p)

    if opts.recall == []:
        sys.exit('\nERROR: must enter a valid --recall.')
    
    print '\nOPTIONS ARE', opts
    print

    # READ FILES
    print 'reading all files...'
    #EdgePR, NodePR, = readPRs(opts.prefix)
    EdgeRankedList, NodeRankedList, = readRankedLists(opts.prefix)
    hopdict = readHops()


    ## EDGES
    if opts.edges:
        print '  Drawing edges precision/recall plots:'
        titleStub = 'Predicted Edges' 
        for recall in opts.recall:
            out = opts.out+'recall=%.2f-' % (recall)
            for p in EdgeRankedList:
                if p == 'aggregate':
                    plotFalseNegatives(EdgeRankedList[p],hopdict[p],recall,opts.pdf,out,edges=True,aggregate=True)
                else:
                    plotFalseNegatives(EdgeRankedList[p],hopdict[p],recall,opts.pdf,out,edges=True)

    ##NODES
    if opts.nodes:
        print '  Drawing nodes precision/recall plots:'
        titleStub = 'Predicted Nodes'
        for recall in opts.recall:
            out = opts.out+'recall=%.2f-' % (recall)
            for p in NodeRankedList:
                if p == 'aggregate':
                    plotFalseNegatives(NodeRankedList[p],hopdict[p],recall,opts.pdf,out,aggregate=True)
                else:
                    plotFalseNegatives(NodeRankedList[p],hopdict[p],recall,opts.pdf,out)


if __name__=='__main__':
    main(sys.argv)
