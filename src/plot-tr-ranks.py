from optparse import OptionParser
import sys
from utilsPoirel import *
import networkx as nx
import matplotlib 
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

## http://matplotlib.org/examples/pylab_examples/broken_axis.html

############################################################
## getNetPathPathways reads the pathways for NetPath and returns
## them as a list.
## onlynetpathwnt: if True, only return Wnt
## dbcompare: if True, only return the 6 pathways in common with KEGG
def getNetPathPathways(onlynetpathwnt,overlapwithkegg):
    if onlynetpathwnt:
        return ['Wnt']
    if overlapwithkegg: # only return the 6 pathways in common with KEGG.
        analyzedpathwayfile = 'data/netpath-dbcompare-pathways.txt'
    else:
        analyzedpathwayfile = 'data/netpath-analyzed-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,1)]
    return pathways

############################################################
## getKEGGPathways reads the pathways for KEGG and returns
## them as a list. It also returns a dictionary of {keggid:netpathname}
def getKEGGPathways(overlapwithnetpath):
    if overlapwithnetpath:
        analyzedpathwayfile = 'data/kegg-dbcompare-pathways.txt'
    else:
        analyzedpathwayfile = 'data/kegg-analyzed-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,2)]
    # dictionary of keggnames to netpath names.
    kegg2netpath = readDict(analyzedpathwayfile,2,1)
    return pathways,kegg2netpath

#######################################################################
def readNodes(datadir,pathway,kegg):
    if pathway != 'aggregate':
        nodefile = '%s/%s-nodes.txt' % (datadir,pathway)
        tfs = {u:n for u,t,n in readColumns(nodefile,1,2,3) if t == 'tf'}
        receptors = {u:n for u,t,n in readColumns(nodefile,1,2,3) if t == 'receptor'}
    else:
        tfs = {}
        receptors = {}
        if kegg:
            pathways,kegg2netpath = getKEGGPathways(False)
        else:
            pathways = getNetPathPathways(False,False)
        for p in pathways:
            nodefile = '%s/%s-nodes.txt' % (datadir,p)
            tfs.update({(p,u):n for u,t,n in readColumns(nodefile,1,2,3) if t == 'tf'})
            receptors.update({(p,u):n for u,t,n in readColumns(nodefile,1,2,3) if t == 'receptor'})
    print '%d tfs and %d receptors.' % (len(tfs),len(receptors))
    return tfs,receptors

#######################################################################
def readPathLinkerEdges(tfs,receptors,indir,pathway,kegg):
    pl_tfs = {}
    pl_receptors = {}
    if pathway != 'aggregate':
        infile = '%s/pathlinker/%s-k_20000-ranked-edges.txt' % (indir,pathway)
        lines = readColumns(infile,1,2,3)
        for u,v,ksp in lines:
            if u in tfs and u not in pl_tfs:
                pl_tfs[u] = int(ksp)
            if v in tfs and v not in pl_tfs:
                pl_tfs[v] = int(ksp)
            if u in receptors and u not in pl_receptors:
                pl_receptors[u] = int(ksp)
            if v in receptors and v not in pl_receptors:
                pl_receptors[v] = int(ksp)
    else:
        if kegg:
            pathways,kegg2netpath = getKEGGPathways(False)
        else:
            pathways = getNetPathPathways(False,False)
        for p in pathways:
            infile = '%s/pathlinker/%s-k_20000-ranked-edges.txt' % (indir,p)
            lines = readColumns(infile,1,2,3)
            for u,v,ksp in lines:
                if (p,u) in tfs and (p,u) not in pl_tfs:
                    pl_tfs[(p,u)] = int(ksp)
                if (p,v) in tfs and (p,v) not in pl_tfs:
                    pl_tfs[(p,v)] = int(ksp)
                if (p,u) in receptors and (p,u) not in pl_receptors:
                    pl_receptors[(p,u)] = int(ksp)
                if (p,v) in receptors and (p,v) not in pl_receptors:
                    pl_receptors[(p,v)] = int(ksp)
    print '%d pathlinker TFs and %d pathlinker receptors' % (len(pl_tfs),len(pl_receptors))
    return pl_tfs,pl_receptors

#######################################################################
def readPageRankNodes(tfs,receptors,indir,pathway,kegg):
    pr_tfs = {}
    pr_receptors = {}
    if pathway != 'aggregate':
        infile = '%s/pagerank/%s-q_0.50-node-pagerank.txt' % (indir,pathway)
        lines = readColumns(infile,1,2)
        i=1
        for n,p in sorted(lines,key=lambda x:float(x[1]),reverse=True):
            if n in tfs and n not in pr_tfs:
                pr_tfs[n] = i
            if n in receptors and n not in pr_receptors:
                pr_receptors[n] = i
            i+=1
    else:
        if kegg:
            pathways,kegg2netpath = getKEGGPathways(False)
        else:
            pathways = getNetPathPathways(False,False)
        for p in pathways:
            infile = '%s/pagerank/%s-q_0.50-node-pagerank.txt' % (indir,p)
            lines = readColumns(infile,1,2)
            i=1
            for n,v in sorted(lines,key=lambda x:float(x[1]),reverse=True):
                if (p,n) in tfs and (p,n) not in pr_tfs:
                    pr_tfs[(p,n)] = i
                if (p,n) in receptors and (p,n) not in pr_receptors:
                    pr_receptors[(p,n)] = i
                i+=1
    print '%d pagerank TFs and %d pagerank receptors' % (len(pr_tfs),len(pr_receptors))
    return pr_tfs,pr_receptors
  
#######################################################################
# no longer used.  Plots ranks
# def plotRanks(tfs,pathlinker,pagerank,ipa,outprefix):
#     iparanks = [pagerank[k] for k in ipa]
#     maxval = max(max(pagerank.values()),max(pathlinker.values()))
#     ply = .3
#     pry = 0
#     ipay = -.1
#     # If we were to simply plot pts, we'd lose most of the interesting
#     # details due to the outliers. So let's 'break' or 'cut-out' the y-axis
#     # into two portions - use the top (ax) for the outliers, and the bottom
#     # (ax2) for the details of the majority of our data
#     fig,(ax,ax2) = plt.subplots(1,2,sharey=True,figsize=(18,3))
#     # plot the same data on both axes
#     ax.plot(sorted(pagerank.values()),[pry]*len(pagerank),'s',ms=8)
#     ax.plot(sorted(pathlinker.values()),[ply]*len(pathlinker),'s',ms=8)
#     ax.plot(sorted(iparanks),[ipay]*len(iparanks),'r*',ms=8)
#     for tf in pathlinker.keys():
#         ax.plot([pagerank[tf],pathlinker[tf]],[pry+0.01,ply-0.01],'-ok',ms=4)
#     ax2.plot(sorted(pagerank.values()),[pry]*len(pagerank),'s',ms=8)
#     ax2.plot(sorted(pathlinker.values()),[ply]*len(pathlinker),'s',ms=8)
#     ax2.plot(sorted(iparanks),[ipay]*len(iparanks),'r*',ms=8)
#     for tf in pathlinker.keys():
#         ax2.plot([pagerank[tf],pathlinker[tf]],[pry+0.01,ply-0.01],'-ok',ms=4)
#     # zoom-in / limit the view to different portions of the data
#     ax.set_xlim(0,1300)
#     ax.set_ylim(ipay-0.05,ply+0.05)
#     ax2.set_xlim(maxval-1000,maxval+50) # outliers only
#     ax2.set_ylim(ipay-0.05,ply+0.05)
#     # hide the spines between ax and ax2
#     ax.spines['right'].set_visible(False)
#     ax.yaxis.tick_left()
#     ax.yaxis.set_ticks([ipay,pry,ply])
#     ax.yaxis.set_ticklabels(['IPA','PageRank','PathLinker'])
#     ax.tick_params(labelright='off') # don't put tick labels at the top
#     ax2.spines['left'].set_visible(False)
#     # This looks pretty good, and was fairly painless, but you can get that
#     # cut-out diagonal lines look with just a bit more work. The important
#     # thing to know here is that in axes coordinates, which are always
#     # between 0-1, spine endpoints are at these locations (0,0), (0,1),
#     # (1,0), and (1,1).  Thus, we just need to put the diagonals in the
#     # appropriate corners of each of our axes, and so long as we use the
#     # right transform and disable clipping.
#     d = .015 # how big to make the diagonal lines in axes coordinates
#     # arguments to pass plot, just so we don't keep repeating them
#     kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
#     ax.plot((1-d,1+d),(1-d,1+d), **kwargs)      # top-left diagonal
#     ax.plot((1-d,1+d),(-d,+d), **kwargs)    # top-right diagonal
#     kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
#     ax2.plot((-d,+d),(1-d,1+d), **kwargs)   # bottom-left diagonal
#     ax2.plot((-d,+d),(-d,+d), **kwargs) # bottom-right diagonal
#     #plt.tight_layout()
#     plt.savefig(outprefix+'.png')
#     print 'Wrote to '+outprefix+'.png'
#     return

#######################################################################
def plotDistribution(tfs,receptors,pl_tfs,pl_receptors,pr_tfs,pr_receptors,outprefix,pathway, truncate):
    

    ## from plot-precision-recall.py
    plcolor = '#009933'  # medium green
    prcolor = '#99CCFF'  # light blue
    receptorcolor = '#8CE1DC' #blue
    tfcolor = '#FFFF60' #yellow

    fig = plt.figure(figsize=(5,5))

    xmax = max(pl_receptors.values() + pr_receptors.values() + pl_tfs.values() + pr_tfs.values())
    print 'max X value is %d' % (xmax)
    ## plot receptors
    ax = plt.subplot(2,1,1)
    ax.plot(sorted(pl_receptors.values()),range(1,len(pl_receptors)+1),'-d',color=plcolor,
            label='PathLinker Receptors',lw=2)
    ax.plot(sorted(pr_receptors.values()),range(1,len(pr_receptors)+1),'-d',color=prcolor,
            label='RWR Receptors',lw=2)   
    ax.set_title('Receptors in the %s reconstruction' % (pathway))
    ax.set_xlabel('Protein Ranking',size=10)
    ax.set_ylabel('# Receptors',size=10)
    ax.legend(loc='lower right', prop={'size':8}, numpoints=1)


    ax2 = plt.subplot(2,1,2)
    ax2.plot(sorted(pl_tfs.values()),range(1,len(pl_tfs)+1),'-s',color=plcolor,label='PathLinker TRs',lw=2)
    ax2.plot(sorted(pr_tfs.values()),range(1,len(pr_tfs)+1),'-s',color=prcolor,label='RWR TRs',lw=2)
    ax2.set_title('TRs in the %s reconstruction' % (pathway))
    ax2.set_xlabel('Protein Ranking',size=10)
    ax2.set_ylabel('# TRs',size=10)
    ax2.legend(loc='lower right', prop={'size':8}, numpoints=1)

    ax.set_ylim([1,len(receptors)+1])
    ax.set_yticks([1,len(receptors)])
    ax2.set_ylim([1,len(tfs)+1])
    ax2.set_yticks([1,len(tfs)])

    plt.tight_layout()
    plt.savefig(outprefix+'-distribution.png')
    print 'Wrote to '+outprefix+'-distribution.png'
    plt.savefig(outprefix+'-distribution.pdf')
    print 'Wrote to '+outprefix+'-distribution.pdf'

    if not truncate:
        ax.set_xlim([-10,xmax+10])
        ax2.set_xlim([-10,xmax+10])  
        plt.savefig(outprefix+'-distribution-fixed-xaxis.png')
        print 'Wrote to '+outprefix+'-distribution-fixed-xaxis.png'
        plt.savefig(outprefix+'-distribution-fixed-xaxis.pdf')
        print 'Wrote to '+outprefix+'-distribution-fixed-xaxis.pdf'
    else:
        ax.set_xlim([-10,truncate+10])
        ax2.set_xlim([-10,truncate+10])
        plt.savefig(outprefix+'-distribution-fixed-xaxis-%d.png' % (truncate))
        print 'Wrote to '+outprefix+'-distribution-fixed-xaxis-%d.png'  % (truncate)
        plt.savefig(outprefix+'-distribution-fixed-xaxis-%d.pdf' % (truncate))
        print 'Wrote to '+outprefix+'-distribution-fixed-xaxis-%d.pdf' % (truncate)
    return

#######################################################################
def main(args):
    usage = '''plot-tr-ranks.py [options]
'''
    parser = OptionParser(usage=usage)

    # General Options
    parser.add_option('','--datadir',type='string',metavar='STR',\
                      help='Data directory for node files.')
    parser.add_option('-o', '--outprefix', type='string', metavar='STR',\
        help='A string to prepend to all output files. Required.')
    parser.add_option('','--indir',type='string',metavar='STR',\
                      help='Input directory.')
    parser.add_option('','--pathway',type='string',metavar='STR',\
                      help='pathway or "aggregate"')
    parser.add_option('','--kegg',action='store_true',default=False,\
                      help='run kegg instead of netpath.')
    parser.add_option('','--truncate',type='int',metavar='INT',\
                      help='Truncate x-axis to this value.')

    # parse the command line arguments
    (opts, args) = parser.parse_args()
    if opts.indir == None:
        sys.exit('ERROR: input directory required.')
    if opts.outprefix == None:
        sys.exit('ERROR: output prefix required.')
    if opts.datadir == None:
        sys.exit('ERROR: data directory required.')
    if opts.pathway == None:
        sys.exit('ERROR: must specify pathway or "aggregate".')

    print '\nOPTIONS ARE', opts
    
    # read nodes
    tfs,receptors = readNodes(opts.datadir,opts.pathway,opts.kegg)
  
    out = open(opts.outprefix+'.txt','w')
    out.write('#%d tfs and %d receptors\n' % (len(tfs),len(receptors)))
    if opts.pathway!='aggregate':
        out.write('#id\tname\ttf/receptor\tpathlinker/pagerank\tval\n')
    else:
        out.write('#pathway\tid\tname\ttf/receptor\tpathlinker/pagerank\tval\n')

    # read pathlinker edges
    pl_tfs,pl_receptors = readPathLinkerEdges(tfs,receptors,opts.indir,opts.pathway,opts.kegg)
    if opts.pathway!='aggregate':
        for tf in pl_tfs:
            out.write('%s\t%s\ttf\tpathlinker\t%d\n' % (tf,tfs[tf],pl_tfs[tf]))
        for receptor in pl_receptors:
            out.write('%s\t%s\treceptor\tpathlinker\t%d\n' % (receptor,receptors[receptor],pl_receptors[receptor]))
    else:
        for p,tf in pl_tfs:
            out.write('%s\t%s\t%s\ttf\tpathlinker\t%d\n' % (p,tf,tfs[(p,tf)],pl_tfs[(p,tf)]))
        for p,receptor in pl_receptors:
            out.write('%s\t%s\t%s\treceptor\tpathlinker\t%d\n' % (p,receptor,receptors[(p,receptor)],pl_receptors[(p,receptor)]))

    # read pagerank edges
    pr_tfs,pr_receptors = readPageRankNodes(tfs,receptors,opts.indir,opts.pathway,opts.kegg)
    if opts.pathway!='aggregate':
        for tf in pr_tfs:
            out.write('%s\t%s\ttf\tpagerank\t%d\n' % (tf,tfs[tf],pr_tfs[tf]))
        for receptor in pl_receptors:
            out.write('%s\t%s\treceptor\tpagerank\t%d\n' % (receptor,receptors[receptor],pr_receptors[receptor]))
    else:
        for p,tf in pr_tfs:
            out.write('%s\t%s\t%s\ttf\tpagerank\t%d\n' % (p,tf,tfs[(p,tf)],pr_tfs[(p,tf)]))
        for p,receptor in pl_receptors:
            out.write('%s\t%s\t%s\treceptor\tpagerank\t%d\n' % (p,receptor,receptors[(p,receptor)],pr_receptors[(p,receptor)]))
        
    out.close()
    print 'wrote to %s.txt' % (opts.outprefix)

    ## only get receptors/tfs that appear in pathlinker results
    if opts.pathway!='aggregate':
        tfs = {t:tfs[t] for t in pr_tfs.keys()}
        receptors = {r:receptors[r] for r in pr_receptors.keys()}
    else:
        tfs = {t:tfs[t] for t in pr_tfs.keys()}
        receptors = {r:receptors[r] for r in pr_receptors.keys()}
      
    #plotRanks(tfs,pathlinkertfs,pageranktfs,ipa,opts.outprefix)
    plotDistribution(tfs,receptors,pl_tfs,pl_receptors,pr_tfs,pr_receptors,opts.outprefix,opts.pathway,opts.truncate)

    
if __name__=='__main__':
    main(sys.argv)
