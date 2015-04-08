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

VERSION_NAMES = {
    '2015pathlinker': '"regulation of signal transduction"',
    'pathlinker-signaling': '"signal transduction"', 
    'pathlinker-signaling-children': 'Children of "signal transduction"',   
    'pathlinker-signaling-children-reg': 'Children of "signal transduction"\nplus "regulation of signal transduction"',   
}

PRDIR = 'results/%s/weighted/netpath/precision-recall/pathlinker/Wnt-exclude_%s-sample_50X-%s-precision-recall.txt'
ORIGDIR = '/data/annaritz/signaling/2014-11-weighted-interactome/weighted-viz/precrecfiles-sample-once-per-pathway/varyparams-exclude_%s-KSP-Wnt-ranked-%ss-pr.txt'
PAPERDIR = '/data/annaritz/signaling/2014-06-linker/viz/precrecfiles-sample-once-per-pathway/precrec-exclude_%s-PRflux+KSP-Wnt-ranked-%ss-pr.txt'

AGGPRDIR = 'results/%s/weighted/netpath/precision-recall/pathlinker/aggregate-exclude_%s-sample_50X-%s-precision-recall.txt'
AGGORIGDIR = '/data/annaritz/signaling/2014-11-weighted-interactome/weighted-viz/precrecfiles-sample-once-per-pathway/varyparams-exclude_%s-KSP-0.5-aggregate-%ss-pr.txt'
AGGORIDNEWRUNDIR = 'results/old-weighted-ksp//aggregate-exclude_%s-sample_50X-%s-precision-recall.txt'
AGGPAPERDIR = '/data/annaritz/signaling/2014-06-linker/viz/precrecfiles-sample-once-per-pathway/precrec-exclude_%s-PRflux+KSP-aggregate-%ss-pr.txt'

COLORS =  {
    '2015pathlinker': 'r',
    'pathlinker-signaling': 'g', 
    'pathlinker-signaling-children': 'b',   
    'pathlinker-signaling-children-reg':  'c',   
}

def plotPR():
    f = 'viz/sanity-checks/wnt-weighted-versions-of-ppi-pathlinker'
    fig = plt.figure(figsize=(8,8))
    
    subplot = 1
    for display in ['node','edge']:
        for negtype in ['none','adjacent']:
            ax = fig.add_subplot(2,2,subplot)
            subplot+=1
            ax.set_xlabel('Recall', size=12)
            ax.set_ylabel('Precision', size=12)
            if negtype == 'adjacent':
                ax.set_title('PathLinker %s (k=20000)\nExclude Pathway-Adj. Negs' % (display),size=12)
            else:
                ax.set_title('PathLinker %s (k=20000)' % (display),size=12)
            for version in VERSION_NAMES:
                if display == 'node':
                    pr = readColumns(PRDIR % (version,negtype,display),4,5)
                else:
                    pr = readColumns(PRDIR % (version,negtype,display),5,6)
                ax.plot([float(r) for p,r in pr],[float(p) for p,r in pr],COLORS[version],label=VERSION_NAMES[version])

            # # plot originals 
            pr = readColumns(ORIGDIR % (negtype,display),1,2)
            ax.plot([float(r) for p,r in pr],[float(p) for p,r in pr],color=[.8,.8,.8],lw=2,\
                    label='KSP on old weighted network,\nOld PR calculation')

            pr = readColumns(PAPERDIR % (negtype,display),1,2)
            ax.plot([float(r) for p,r in pr],[float(p) for p,r in pr],'k',lw=2,label='PRFlux+KSP original paper')
                    
            if subplot==3:
                ax.legend(loc='lower left', prop={'size':8}, numpoints=1)
            #else:
            #    ax.legend(loc='upper right', prop={'size':8}, numpoints=1)
            ax.set_xlim([-.01,1.01])
            ax.set_ylim([-.01,1.01])

    plt.tight_layout()
    figname = '%s.png' %(f)
    print '\twriting to %s' % (figname)
    plt.savefig(figname)

    return

def plotAggPR():
    f = 'viz/sanity-checks/aggregate-weighted-versions-of-ppi-pathlinker'
    fig = plt.figure(figsize=(8,8))
    
    subplot = 1
    for display in ['node','edge']:
        for negtype in ['none','adjacent']:
            ax = fig.add_subplot(2,2,subplot)
            subplot+=1
            ax.set_xlabel('Recall', size=12)
            ax.set_ylabel('Precision', size=12)
            if negtype == 'adjacent':
                ax.set_title('PathLinker %s (k=20000)\nExclude Pathway-Adj. Negs' % (display),size=12)
            else:
                ax.set_title('PathLinker %s (k=20000)' % (display),size=12)
            for version in VERSION_NAMES:
                if display == 'node':
                    pr = readColumns(AGGPRDIR % (version,negtype,display),5,6)
                else:
                    pr = readColumns(AGGPRDIR % (version,negtype,display),6,7)
                ax.plot([float(r) for p,r in pr],[float(p) for p,r in pr],COLORS[version],label=VERSION_NAMES[version])

            # # plot originals 
            pr = readColumns(AGGORIGDIR % (negtype,display),1,2)
            ax.plot([float(r) for p,r in pr],[float(p) for p,r in pr],color=[.8,.8,.8],lw=2,\
                    label='KSP on old weighted network,\nOld PR calculation')

            if display == 'node':
                pr = readColumns(AGGORIDNEWRUNDIR % (negtype,display),5,6)
            else:
                pr = readColumns(AGGORIDNEWRUNDIR % (negtype,display),6,7)
            ax.plot([float(r) for p,r in pr],[float(p) for p,r in pr],'--',color=[.5,.5,.5],lw=2,label='KSP on old weighted network,\nNew PR Calculation')

            pr = readColumns(AGGPAPERDIR % (negtype,display),1,2)
            ax.plot([float(r) for p,r in pr],[float(p) for p,r in pr],'k',lw=2,label='PRFlux+KSP original paper')
                    
            if subplot==3:
                ax.legend(loc='lower left', prop={'size':8}, numpoints=1)
            #else:
            #    ax.legend(loc='upper right', prop={'size':8}, numpoints=1)
            ax.set_xlim([-.01,1.01])
            ax.set_ylim([-.01,1.01])

    plt.tight_layout()
    figname = '%s.png' %(f)
    print '\twriting to %s' % (figname)
    plt.savefig(figname)

    return

def main(args):
    plotPR()
    plotAggPR()

if __name__=='__main__':
    main(sys.argv)
