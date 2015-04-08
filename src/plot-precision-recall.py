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
              'nmax':[15, 25, 35, 50, 75, 100, 200, 500],
          }

COLORS = {
    'pathlinker':'c',
    'pagerank':  'r',
    'shortestpaths': [1,.8,.8],
    'inducedsubgraph':'k',
    'eqed':'k',
    'responsenet':'g',
    'pcsf':'b',
    'anat':'m',
    'ipa':[.8,.8,1],
}

SHAPES = {
    'shortestpaths': 's',
    'responsenet':'^',
    'pcsf':'v',
    'anat':'d',
    'ipa':'o',
    'inducedsubgraph':'v'
}

##################################################################
def plotPR(methodlist,precrecs,f,name,pdf):
    fig = plt.figure(figsize=(12,12))

    subplot = 1
    titles = {
        'node': {
            'none': '%s Pathway Nodes' % (name),
            'adjacent':'%s Pathway Nodes,\nExclude Pathway-Adjacent Negatives' % (name)
        },
        'edge': {
            'none': '%s Pathway Edges' % (name),
            'adjacent':'%s Pathway Edges,\nExclude Pathway-Adjacent Negatives' % (name)
        },
    }

    for display in ['node','edge']:
        for negtype in ['none','adjacent']:
            ax = fig.add_subplot(2,2,subplot)
            subplot+=1
            ax.set_xlabel('Recall', size=12)
            ax.set_ylabel('Precision', size=12)
            ax.set_title(titles[display][negtype],size=12)
            increments = {alg:0.0 for alg in COLORS.keys()}
            for alg in methodlist:
                color,linewidth,markersize,markerstyle = getstyle(alg,increments,len(methodlist))
                pr = precrecs[display][negtype][alg]

                if len(pr)==1: # plot single point
                    ax.plot([r for p,r in pr], [p for p,r in pr],markerstyle,ms=markersize,color=color,label=alg)
                elif alg == 'ipa': # plot points AND line
                    ax.plot([r for p,r in pr],[p for p,r in pr],'--'+markerstyle,ms=markersize,lw=linewidth,color=color,label=alg)
                else: # plot line only.
                    ax.plot([r for p,r in pr], [p for p,r in pr],lw=linewidth,color=color,label=alg)

            if len < 20:
                    # put legend outside plot
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width, box.height * 0.7])
                if negtype=='none':
                    ax.legend(loc='center left', ncol=5, frameon=False, prop={'size':8}, bbox_to_anchor=(0, 1.32), numpoints=1)
            else:
                ax.legend(loc='lower left', ncol=3,frameon=False, prop={'size':8}, numpoints=1)

            ax.set_ylim([-.01,1.01])
            ax.set_xlim([-.01,1.01])

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
        markersize = 6
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
def readFiles(varyparams,algs,indir,pathway):
    if varyparams:
        filepatterns = FILELOCATIONS['varyparams']
    else:
        filepatterns = FILELOCATIONS['singleparam']

    # shift variable shifts the column if it's aggregate.
    if pathway == 'aggregate': # first column is pathway
        shift = 1
    else: 
        shift = 0

    nodeprecrec = {}
    edgeprecrec = {}
    methodlist = []
    for negtype in ['none','adjacent']:
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

    ## make all lists unique, but preserve order.
    for negtype in ['none','adjacent']:
        for key in nodeprecrec[negtype].keys():
            nodeprecrec[negtype][key] = list(OrderedDict.fromkeys(nodeprecrec[negtype][key]))
            edgeprecrec[negtype][key] = list(OrderedDict.fromkeys(edgeprecrec[negtype][key]))

    return nodeprecrec,edgeprecrec,methodlist

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

    ## Read files for nodes and edges
    precrecs = {}
    precrecs['node'],precrecs['edge'],methodlist = readFiles(opts.varyparams,opts.alg,opts.indir,opts.pathway)

    if not opts.varyparams:
        ## plot PR
        plotPR(methodlist,precrecs,opts.outprefix,opts.pathway,opts.pdf)
    else:
        ## plot PR for every algorithm specified.
        for alg in opts.alg:
            outprefix = '%s-varyparams-%s' % (opts.outprefix,alg)
            parammethods = [m for m in methodlist if alg in m]
            parammethods.reverse()
            pr = {'node':{},'edge':{}}
            pr['node']['none'] = {key:precrecs['node']['none'][key] for key in parammethods}
            pr['node']['adjacent'] = {key:precrecs['node']['adjacent'][key] for key in parammethods}
            pr['edge']['none'] = {key:precrecs['edge']['none'][key] for key in parammethods}
            pr['edge']['adjacent'] = {key:precrecs['edge']['adjacent'][key] for key in parammethods}
            plotPR(parammethods,pr,outprefix,opts.pathway,opts.pdf)

if __name__=='__main__':
    main(sys.argv)
