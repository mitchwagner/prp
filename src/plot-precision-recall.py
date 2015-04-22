#!/usr/bin/python

## Modified from /data/poirel/research/signaling-pathways/viz/precision-recall-faster/plot-precision-recall.py

import os
import sys
from optparse import OptionParser, OptionGroup
from utilsPoirel import *
from scipy.stats.mstats import rankdata
import matplotlib 
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import os.path
import random

## To make make a list unique while preserving order:
## http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order
from collections import OrderedDict
  
# dictionary of {name:file location}
FILELOCATIONS = {
    'singleparam': {  # "default" parameters
        'pathlinker':'%s/pathlinker/%s-exclude_%s-sample_50X-%s-precision-recall.txt',
        'pagerank':  '%s/pagerank/%s-q_0.50-exclude_%s-sample_50X-%s-precision-recall.txt',
        'shortestpaths': '%s/shortestpaths/%s-exclude_%s-sample_50X-%s-precision-recall.txt',
        'inducedsubgraph':'%s/inducedsubgraph/%s-exclude_%s-sample_50X-%s-precision-recall.txt',
        'eqed':'%s/eqed/%s-exclude_%s-sample_50X-%s-precision-recall.txt',
        'responsenet':'%s/responsenet/%s-gamma_20-exclude_%s-sample_50X-%s-precision-recall.txt',
        'pcsf':'%s/pcsf/%s-prize5-omega0.01-exclude_%s-sample_50X-%s-precision-recall.txt',
        'anat':'%s/anat/%s-alpha0.00-exclude_%s-sample_50X-%s-precision-recall.txt',
        'ipa':'%s/ipa/%s-nmax35-exclude_%s-sample_50X-%s-precision-recall.txt',
    },
    'varyparams': { # variable parameters 
        'pathlinker':'%s/pathlinker/%s-exclude_%s-sample_50X-%s-precision-recall.txt',
        'pagerank':  '%s/pagerank/%s-q_%.2f-exclude_%s-sample_50X-%s-precision-recall.txt',
        'shortestpaths': '%s/shortestpaths/%s-exclude_%s-sample_50X-%s-precision-recall.txt',
        'inducedsubgraph':'%s/inducedsubgraph/%s-exclude_%s-sample_50X-%s-precision-recall.txt',
        'eqed':'%s/eqed/%s-exclude_%s-sample_50X-%s-precision-recall.txt',
        'responsenet':'%s/responsenet/%s-gamma_%d-exclude_%s-sample_50X-%s-precision-recall.txt',
        'pcsf':'%s/pcsf/%s-prize%d-omega%.2f-exclude_%s-sample_50X-%s-precision-recall.txt',
        'anat':'%s/anat/%s-alpha%.2f-exclude_%s-sample_50X-%s-precision-recall.txt',
        'ipa':'%s/ipa/%s-nmax%d-exclude_%s-sample_50X-%s-precision-recall.txt',
    },
}

## VARYPARAMS
VARYPARAMS = {'q': [0.1, 0.25, 0.5, 0.75, 0.9],
              'gamma':[10, 15, 20, 25, 30],
              'omega':[0, 0.01, 0.1],
              'prize':[1, 3, 5, 7, 9],
              'alpha':[0, 0.1, 0.25, 0.4, 0.5],
              'nmax':[5,10,15, 25, 35, 50, 75, 100, 200, 500],
          }

# Set the color of the plot for each algorithm.
# TODO: Write the names of these colors down, or find comparable named HTML
# colors and replace the cryptic HTML color codes.
COLORS = {
    'pagerank' : '#99CCFF', # light blue
    'pathlinker' : '#009933',  # medium green
    'eqed' : '#0404B4', 
    'pcsf' : '#CC3300', # red
    'anat' : 'b', # blue
    'responsenet' : '#FE9A2E', # orange
     'inducedsubgraph': '#66CCFF',
    'shortestpaths':'#F4FA58', # yellow
    'ipa':'#F78181', # pink
}

SHAPES = {
    'shortestpaths': 's',
    'responsenet':'^',
    'pcsf':'v',
    'anat':'d',
    'ipa':'o',
    'inducedsubgraph':'v'
}

NAMES = {
    'pagerank' : 'RWR',
    'pathlinker' : 'PathLinker', 
    'eqed' : 'eQED',
    'pcsf' : 'PCSF',
    'anat' : 'ANAT',
    'responsenet' : 'ResponseNet',
     'inducedsubgraph': 'InducedSubgraph',
    'shortestpaths':'ShortestPaths',
    'ipa':'IPA',
}

##################################################################
def getTitle(name,ignorekegg,ignorenetpath,numpathways):
    if name != 'aggregate':
        if ignorekegg:
            titles = {
                'node': {
                    'none': 'Proteins in the %s Reconstruction' % (name),
                    'adjacent':'Proteins in the %s Reconstruction\nIgnoring Pathway-Adjacent Negatives' % (name),
                    'file':'Proteins in the NetPath %s Reconstruction\nIgnoring Proteins in KEGG' % (name)
                },
                'edge': {
                    'none': 'Interactions in the %s Reconstruction' % (name),
                    'adjacent':'Interactions in the %s Reconstruction\nIgnoring Pathway-Adjacent Negatives' % (name),
                    'file':'Interactions in the  NetPath%s Reconstruction\nIgnoring Interactions in KEGG' % (name)
                },
            }
        else:
            titles = {
                'node': {
                    'none': 'Proteins in the %s Reconstruction' % (name),
                    'adjacent':'Proteins in the %s Reconstruction\nIgnoring Pathway-Adjacent Negatives' % (name),
                    'file':'Proteins in the KEGG %s Reconstruction\nIgnoring Proteins in NetPath' % (name)
                },
                'edge': {
                    'none': 'Interactions in the %s Reconstruction' % (name),
                    'adjacent':'Interactions in the %s Reconstruction\nIgnoring Pathway-Adjacent Negatives' % (name),
                    'file':'Interactions in the KEGG %s Reconstruction\nIgnoring Interactions in NetPath' % (name)
                },
            }
    else:
        if ignorekegg:
            titles = {
                'node': {
                    'none': 'Proteins Aggregated over \n%d Reconstructions' % (numpathways),
                    'adjacent':'Proteins Aggregated over %s Reconstructions\nIgnoring Pathway-Adjacent Negatives' % (numpathways),
                    'file':'Proteins Aggregated over %s Reconstructions\nIgnoring Proteins in KEGG' % (numpathways)
            },
                'edge': {
                    'none': 'Interactions Aggregated over \n%d Reconstructions'% (numpathways),
                    'adjacent':'Interactions Aggregated over %d Reconstructions\nIgnoring Pathway-Adjacent Negatives'% (numpathways),
                    'file':'Interactions Aggregated over %d Reconstructions\nIgnoring Interactions in KEGG'% (numpathways)
                },
            }
        else:
            titles = {
                'node': {
                    'none': 'Proteins Aggregated over \n%d Reconstructions' % (numpathways),
                    'adjacent':'Proteins Aggregated over %d Reconstructions\nIgnoring Pathway-Adjacent Negatives'% (numpathways),
                    'file':'Proteins Aggregated over %d Reconstructions\nIgnoring Proteins in NetPath'% (numpathways)
            },
                'edge': {
                    'none': 'Interactions Aggregated over \n%d Reconstructions'% (numpathways),
                    'adjacent':'Interactions Aggregated over %d Reconstructions\nIgnoring Pathway-Adjacent Negatives'% (numpathways),
                    'file':'Interactions Aggregated over %d Reconstructions\nIgnoring Interactions in NetPath'% (numpathways)
                },
            }
    return titles

##################################################################
def plotPR(methodlist,precrecs,f,name,pdf,negtypes,ignorekegg,ignorenetpath,numpathways):
    fig = plt.figure(figsize=(12,12))

    subplot = 1
    titles = getTitle(name,ignorekegg,ignorenetpath,numpathways)

    for display in ['node','edge']:
        for negtype in negtypes:
            ax = fig.add_subplot(2,2,subplot)
            subplot+=1
            ax.set_xlabel('Recall', size=14)
            ax.set_ylabel('Precision', size=14)
            ax.set_title(titles[display][negtype],size=16)
            increments = {alg:0.0 for alg in COLORS.keys()}

            ## if negtype == 'file' or negtype == 'adjacent', add original ('none') lines/points with transparency
            #if negtype == 'file':
            for alg in methodlist:
                color,linewidth,markersize,markerstyle = getstyle(alg,increments,len(methodlist))
                pr = precrecs[display]['none'][alg]

                if len(pr)==1: # plot single point
                    ax.plot([r for p,r in pr], [p for p,r in pr],markerstyle,ms=markersize,color=color,label='_nolegend_',alpha=0.4)
                elif alg == 'ipa': # plot points AND line
                    ax.plot([r for p,r in pr],[p for p,r in pr],'--'+markerstyle,ms=markersize,\
                            lw=1,color=color,label='_nolegend_',alpha=0.4)
                else: # plot line only.
                    ax.plot([r for p,r in pr], [p for p,r in pr],lw=1,color=color,label='_nolegend_',alpha=0.4)

            ## add opaque lines/points for this negtype
            for alg in methodlist:
                color,linewidth,markersize,markerstyle = getstyle(alg,increments,len(methodlist))
                pr = precrecs[display][negtype][alg]

                if len(pr)==1: # plot single point
                    ax.plot([r for p,r in pr], [p for p,r in pr],markerstyle,ms=markersize,color=color,label=NAMES.get(alg,alg))
                elif alg == 'ipa': # plot points AND line
                    ax.plot([r for p,r in pr],[p for p,r in pr],'--'+markerstyle,ms=markersize,\
                            lw=linewidth,color=color,label=NAMES.get(alg,alg))
                else: # plot line only.
                    ax.plot([r for p,r in pr], [p for p,r in pr],lw=linewidth,color=color,label=NAMES.get(alg,alg))
            if display == 'node':
                ax.legend(loc='lower left', ncol=2,frameon=True, prop={'size':12}, numpoints=1)
            else:
                ax.legend(loc='upper right', ncol=2,frameon=True, prop={'size':12}, numpoints=1)
            ax.set_ylim([-.02,1.02])
            ax.set_xlim([-.02,1.02])

    if len(methodlist)<20:
        plt.tight_layout()
    
    figname = '%s.png' %(f)
    print '\twriting to %s' % (figname)
    plt.savefig(figname)
    if pdf:
        figname = '%s.pdf' %(f)
        print '\twriting to %s' % (figname)
        plt.savefig(figname)
        os.system('pdfcrop %s %s' %(figname, figname))

    return

#######################################################################
def getstyle(alg,increments,nummethods):
    if alg in COLORS:
        color = COLORS[alg]
        linewidth = 2
        markersize = 8
        markerstyle=SHAPES.get(alg,'-')
    else:
        ## get original algorithm (withouth parameter annotations)
        key = [a for a in increments if a in alg][0]

        ## shade of blue
        color = [increments[key]/nummethods,increments[key]/nummethods,1]
        color = [min(1,c) for c in color]
        increments[key] += 1
        linewidth = 1+increments[key]
        markersize = 2+10*increments[key]/nummethods
        markerstyle = SHAPES.get(key,'-')
    return color,linewidth,markersize,markerstyle

#######################################################################
# to remove dupliates; use OrderedDict.fromkeys() function.
def readFiles(varyparams,algs,indir,pathway,negtypes,ignorekegg,ignorenetpath):
    if varyparams:
        filepatterns = FILELOCATIONS['varyparams']
    else:
        filepatterns = FILELOCATIONS['singleparam']

    # shift variable shifts the column if it's aggregate.
    if pathway == 'aggregate': # first column is pathway
        shift = 1
    else: 
        shift = 0

    if ignorekegg and pathway == 'aggregate':
        ## only get 6 pathways
        pathway += '-pathways_shared_with_kegg'
    elif ignorenetpath and pathway == 'aggregate':
        pathway += '-pathways_shared_with_netpath'

    nodeprecrec = {}
    edgeprecrec = {}
    methodlist = []
    for negtype in negtypes:
        nodeprecrec[negtype] = {}
        edgeprecrec[negtype] = {}
        for alg in algs:
            print 'Reading %s exclude %s' % (alg,negtype)
            if alg not in filepatterns:
                sys.exit('ERROR: add %s to FILELOCATIONS global variable.'% (alg))

            ## read nodes and edges for each parameter
            if not varyparams: # add one line

                if alg == 'ipa':
                    ## read in all params; concatenate into one line.
                    nodeprecrec[negtype][alg] = []
                    edgeprecrec[negtype][alg] = []
                    for param in VARYPARAMS['nmax']:
                        nodelist,edgelist = populateprecrecs(FILELOCATIONS['varyparams'][alg],indir,pathway,negtype,param,shift)
                        nodeprecrec[negtype][alg] += nodelist
                        edgeprecrec[negtype][alg] += edgelist
                    if negtype == 'none':
                        methodlist.append(alg)
                    continue # skip rest of block.

                ## for all other algorithms, read in node/edge file and get list of precision/recall values.
                nodefile = filepatterns[alg] % (indir,pathway,negtype,'node')
                edgefile = filepatterns[alg] % (indir,pathway,negtype,'edge')
                nodeprecrec[negtype][alg] = [(float(p),float(r)) for t,p,r in \
                                             readColumns(nodefile,3+shift,4+shift,5+shift) if t != 'ignore']
                edgeprecrec[negtype][alg] = [(float(p),float(r)) for t,p,r in \
                                             readColumns(edgefile,4+shift,5+shift,6+shift) if t != 'ignore']
                if negtype == 'none':
                    methodlist.append(alg)
                    if pathway != 'aggregate':
                        numpathways = 1
                    else:
                        numpathways = len(set([p for p in readItemSet(nodefile,1)]))
            else: # add lines for each parameter
                if alg == 'pagerank':
                    for param in VARYPARAMS['q']:
                        nodelist,edgelist = populateprecrecs(filepatterns[alg],indir,pathway,negtype,param,shift)
                        nodeprecrec[negtype][alg+'-q_%.2f' % (param)] = nodelist
                        edgeprecrec[negtype][alg+'-q_%.2f' % (param)] = edgelist
                        if negtype == 'none':
                            methodlist.append(alg+'-q_%.2f' % (param))
                elif alg == 'responsenet':
                    for param in VARYPARAMS['gamma']:
                        nodelist,edgelist = populateprecrecs(filepatterns[alg],indir,pathway,negtype,param,shift)
                        nodeprecrec[negtype][alg+'-gamma_%d' % (param)] = nodelist
                        edgeprecrec[negtype][alg+'-gamma_%d' % (param)] = edgelist
                        if negtype == 'none':
                            methodlist.append(alg+'-gamma_%d' % (param))
                elif alg == 'pcsf':
                    for param1 in VARYPARAMS['prize']:
                        for param2 in VARYPARAMS['omega']:

                            nodelist,edgelist = populateprecrecs(filepatterns[alg],indir,pathway,\
                                                                 negtype,param1,shift,param2=param2)
                            nodeprecrec[negtype][alg+'-prize_%d-omega_%.2f' % (param1,param2)] = nodelist
                            edgeprecrec[negtype][alg+'-prize_%d-omega_%.2f' % (param1,param2)] = edgelist
                            if negtype == 'none':
                                methodlist.append(alg+'-prize_%d-omega_%.2f' % (param1,param2))
                elif alg == 'anat':
                    for param in VARYPARAMS['alpha']:
                        nodelist,edgelist = populateprecrecs(filepatterns[alg],indir,pathway,negtype,param,shift)
                        nodeprecrec[negtype][alg+'-alpha_%.2f' % (param)] = nodelist
                        edgeprecrec[negtype][alg+'-alpha_%.2f' % (param)] = edgelist
                        if negtype == 'none':
                            methodlist.append(alg+'-alpha_%.2f' % (param))
                elif alg == 'ipa':
                    for param in VARYPARAMS['nmax']:
                        nodelist,edgelist = populateprecrecs(filepatterns[alg],indir,pathway,negtype,param,shift)
                        nodeprecrec[negtype][alg+'-nmax_%d' % (param)] = nodelist
                        edgeprecrec[negtype][alg+'-nmax_%d' % (param)] = edgelist
                        if negtype == 'none':
                            methodlist.append(alg+'-nmax_%d' % (param))
                else: # no parameters. 
                    nodefile = filepatterns[alg] % (indir,pathway,negtype,'node')
                    edgefile = filepatterns[alg] % (indir,pathway,negtype,'edge')
                    nodeprecrec[negtype][alg] = [(float(p),float(r)) for t,p,r in \
                                                 readColumns(nodefile,3+shift,4+shift,5+shift) if t != 'ignore']
                    edgeprecrec[negtype][alg] = [(float(p),float(r)) for t,p,r in \
                                                 readColumns(edgefile,4+shift,5+shift,6+shift) if t != 'ignore']
                    if negtype == 'none':
                         methodlist.append(alg)
                         if pathway != 'aggregate':
                             numpathways = 1
                         else:
                             numpathways = len(set([p for p in readItemSet(nodefile,1)]))

    ## make all lists unique, but preserve order.
    for negtype in negtypes:
        for key in nodeprecrec[negtype].keys():
            nodeprecrec[negtype][key] = list(OrderedDict.fromkeys(nodeprecrec[negtype][key]))
            edgeprecrec[negtype][key] = list(OrderedDict.fromkeys(edgeprecrec[negtype][key]))

    return nodeprecrec,edgeprecrec,methodlist,numpathways

#######################################################################
def populateprecrecs(filepattern,indir,pathway,negtype,param1,shift,param2=None):
    ## build node file and edge files using param information.
    if param2==None:
        nodefile = filepattern % (indir,pathway,param1,negtype,'node')
        edgefile = filepattern % (indir,pathway,param1,negtype,'edge')
    else:
        nodefile = filepattern % (indir,pathway,param1,param2,negtype,'node')
        edgefile = filepattern % (indir,pathway,param1,param2,negtype,'edge')

    nodelist = [(float(p),float(r)) for t,p,r in readColumns(nodefile,3+shift,4+shift,5+shift) if t != 'ignore']
    edgelist = [(float(p),float(r)) for t,p,r in readColumns(edgefile,4+shift,5+shift,6+shift) if t != 'ignore']

    return nodelist,edgelist
    
#######################################################################
def main(args):
    usage = '''plot-precision-recally.py [options]
'''
    parser = OptionParser(usage=usage)

    # General Options
    parser.add_option('','--indir',type='string',metavar='STR',\
                      help='input directory containing precision recall files. Required.')
    parser.add_option('-o', '--outprefix', type='string', metavar='STR',\
        help='A string to prepend to all output files. Required.')
    parser.add_option('-p','--pathway',type='string',metavar='STR',\
                      help='Pathway to run or "aggregate". Required.')
    parser.add_option('-a', '--alg', action='append', type='str', default=[],\
        help='An algorithm to plot; this option can be provided multiple times.')
    parser.add_option('', '--pdf', action='store_true', default=False,\
        help='Output images in PDF as well as PNG.')
    parser.add_option('','--varyparams',action='store_true',default=False,\
                      help='Plot different parameters if specified.')
    parser.add_option('','--seed',type='int',default=123456,\
                      help='seed for selecting colors.  Default is 123456.')
    parser.add_option('','--ignorekegg',action='store_true',default=False,\
                      help='Instead of plotting exclude-adjacent, plot exclude-file ignoring KEGG positives.')
    parser.add_option('','--ignorenetpath',action='store_true',default=False,\
                      help='Instead of plotting exclude-adjacent, plot exclude-file ignoring NetPath positives.')

    # parse the command line arguments
    (opts, args) = parser.parse_args()
    if opts.indir == None:
        sys.exit('ERROR: input directory required.')
    if opts.outprefix == None:
        sys.exit('ERROR: output directory required.')
    if len(opts.alg) == 0:
        sys.exit('must specify at least one algorithm with --alg')

    print '\nOPTIONS ARE', opts
    print

    ## set seed
    random.seed(opts.seed)

    ##
    if opts.ignorekegg or opts.ignorenetpath:
        negtypes = ['none','file']
    else:
        negtypes = ['none','adjacent']

    ## Read files for nodes and edges
    precrecs = {}
    precrecs['node'],precrecs['edge'],methodlist,numpathways = readFiles(opts.varyparams,opts.alg,opts.indir,opts.pathway,negtypes,opts.ignorekegg,opts.ignorenetpath)

    if not opts.varyparams:
        ## plot PR
        plotPR(methodlist,precrecs,opts.outprefix,opts.pathway,opts.pdf,negtypes,opts.ignorekegg,opts.ignorenetpath,numpathways)
    else:
        ## plot PR for every algorithm specified.
        for alg in opts.alg:
            outprefix = '%s-varyparams-%s' % (opts.outprefix,alg)
            parammethods = [m for m in methodlist if alg in m]
            parammethods.reverse()
            pr = {'node':{},'edge':{}}
            for negtype in opts.negtypes:
                pr['node'][negtype] = {key:precrecs['node'][negtype][key] for key in parammethods}
                pr['edge'][negtype] = {key:precrecs['edge'][negtype][key] for key in parammethods}
            plotPR(parammethods,pr,outprefix,opts.pathway,opts.pdf,negtypes,opts.ignorekegg,opts.ignorenetpath,numpathways)

if __name__=='__main__':
    main(sys.argv)
