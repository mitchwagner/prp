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
    
#######################################################################
def main(args):
    usage = '''plot-precision-recall.py file1prefix file2prefix outname file1name file2name
'''
    if len(args) != 6:
        print usage
        sys.exit('%d arguments.'%  (len(args)))

    file1prefix = args[1]
    file2prefix = args[2]
    figname = args[3]
    file1label = args[4]
    file2label = args[5]

    fig = plt.figure(figsize=(12,12,))
    subplot = 1
    for display in ['node','edge']:
        for negtype in ['none','adjacent']:
            file1 = '%s-exclude_%s-sample_50X-%s-precision-recall.txt' %  (file1prefix,negtype,display)
            file2 = '%s-exclude_%s-sample_50X-%s-precision-recall.txt' %  (file2prefix,negtype,display)
            if display == 'node':
                try:
                    list1 = [(float(p),float(r)) for t,p,r in \
                             readColumns(file1,3,4,5) if t != 'ignore']
                    list2 = [(float(p),float(r)) for t,p,r in \
                             readColumns(file2,3,4,5) if t != 'ignore']
                except ValueError: # maybe this is aggregate - shift by one column
                    list1 = [(float(p),float(r)) for t,p,r in \
                             readColumns(file1,4,5,6) if t != 'ignore']
                    list2 = [(float(p),float(r)) for t,p,r in \
                             readColumns(file2,4,5,6) if t != 'ignore']
            else:
                try:
                    list1 = [(float(p),float(r)) for t,p,r in \
                             readColumns(file1,4,5,6) if t != 'ignore']
                    list2 = [(float(p),float(r)) for t,p,r in \
                             readColumns(file2,4,5,6) if t != 'ignore']
                except ValueError: # maybe this is aggregate - shift by one column
                    list1 = [(float(p),float(r)) for t,p,r in \
                             readColumns(file1,5,6,7) if t != 'ignore']
                    list2 = [(float(p),float(r)) for t,p,r in \
                             readColumns(file2,5,6,7) if t != 'ignore']
            ax = fig.add_subplot(2,2,subplot)
            subplot+=1
            ax.set_xlabel('Recall', size=12)
            ax.set_ylabel('Precision', size=12)
            ax.set_title('%s exclude %s' % (display,negtype),size=12)
            if len(list1)==1: 
                ax.plot([r for p,r in list1], [p for p,r in list1],'sr',ms=6,label=file1label)
            else: # plot line only.
                ax.plot([r for p,r in list1], [p for p,r in list1],'r',lw=2,label=file1label)
            if len(list2)==1: 
                ax.plot([r for p,r in list2], [p for p,r in list2],'sb',ms=6,label=file2label)
            else: # plot line only.
                ax.plot([r for p,r in list2], [p for p,r in list2],'b',lw=2,label=file2label)

            ax.legend(loc='lower left',frameon=False, prop={'size':10}, numpoints=1)

            ax.set_ylim([-.01,1.01])
            ax.set_xlim([-.01,1.01])

    plt.tight_layout()
    
    print '\twriting to %s' % (figname)
    plt.savefig(figname)

if __name__=='__main__':
    main(sys.argv)
