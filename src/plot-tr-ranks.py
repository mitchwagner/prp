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
        lines = set([(u,v,int(k)) for u,v,k in readColumns(infile,1,2,3)])
        pl_tfs,pl_receptors = walkDownRankings(lines,tfs,receptors,reverse=False,pathway=None)
    else:
        if kegg:
            pathways,kegg2netpath = getKEGGPathways(False)
        else:
            pathways = getNetPathPathways(False,False)
        for p in pathways:
            infile = '%s/pathlinker/%s-k_20000-ranked-edges.txt' % (indir,p)
            lines = set([(u,v,int(k)) for u,v,k in readColumns(infile,1,2,3)])
            ranked_tfs,ranked_receptors = walkDownRankings(lines,tfs,receptors,reverse=False,pathway=p)
            pl_tfs.update(ranked_tfs)
            pl_receptors.update(ranked_receptors)
    print '%d pathlinker TFs and %d pathlinker receptors' % (len(pl_tfs),len(pl_receptors))
    return pl_tfs,pl_receptors

#######################################################################
def readPageRankEdges(tfs,receptors,indir,pathway,kegg):
    pr_tfs = {}
    pr_receptors = {}
    if pathway != 'aggregate':
        infile = '%s/pagerank/%s-q_0.50-edge-fluxes.txt' % (indir,pathway)
        lines = set([(u,v,float(k)) for u,v,k in readColumns(infile,1,2,3)])
        pr_tfs,pr_receptors = walkDownRankings(lines,tfs,receptors,reverse=True,pathway=None)
    else:
        if kegg:
            pathways,kegg2netpath = getKEGGPathways(False)
        else:
            pathways = getNetPathPathways(False,False)
        for p in pathways:
            infile = '%s/pagerank/%s-q_0.50-edge-fluxes.txt' % (indir,p)
            lines = set([(u,v,float(k)) for u,v,k in readColumns(infile,1,2,3)])
            ranked_tfs,ranked_receptors = walkDownRankings(lines,tfs,receptors,reverse=True,pathway=p)
            pr_tfs.update(ranked_tfs)
            pr_receptors.update(ranked_receptors)
    print '%d pagerank TFs and %d pagerank receptors' % (len(pr_tfs),len(pr_receptors))
    return pr_tfs,pr_receptors
  
#######################################################################
def walkDownRankings(lines,tfs,receptors,reverse=True,pathway=None):
    ranked_tfs = {}
    ranked_receptors = {}
    i=1
    prevval = -1
    if pathway == None:
        for u,v,k in sorted(lines,key=lambda x:float(x[2]),reverse=reverse):
            # only increment i if there isn't a tie.
            if prevval != k and prevval != -1:
                i+=1
            prevval=k

            if u in tfs and u not in ranked_tfs:
                ranked_tfs[u] = i
            if v in tfs and v not in ranked_tfs:
                ranked_tfs[v] = i
            if u in receptors and u not in ranked_receptors:
                ranked_receptors[u] = i
            if v in receptors and v not in ranked_receptors:
                ranked_receptors[v] = i

    else: # append pathway to items (for aggregate)
        for u,v,k in sorted(lines,key=lambda x:float(x[2]),reverse=reverse):
            # only increment i if there isn't a tie.
            if prevval != k and prevval != -1:
                i+=1
            prevval=k
            if (pathway,u) in tfs and (pathway,u) not in ranked_tfs:
                ranked_tfs[(pathway,u)] = i
            if (pathway,v) in tfs and (pathway,v) not in ranked_tfs:
                ranked_tfs[(pathway,v)] = i
            if (pathway,u) in receptors and (pathway,u) not in ranked_receptors:
                ranked_receptors[(pathway,u)] = i
            if (pathway,v) in receptors and v not in ranked_receptors:
                ranked_receptors[(pathway,v)] = i
            i+=1

    return ranked_tfs,ranked_receptors

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
    ax.plot(sorted(pl_receptors.values()),range(1,len(pl_receptors)+1),'--d',color=plcolor,
            label='PathLinker Receptors',lw=1)
    ax.plot(sorted(pr_receptors.values()),range(1,len(pr_receptors)+1),'--d',color=prcolor,
            label='RWR Receptors',lw=1)   
    ax.set_title('Receptors in the %s reconstruction' % (pathway))
    ax.set_xlabel('Ranked Interaction',size=12)
    ax.set_ylabel('# Receptors',size=12)
    ax.legend(loc='lower right', prop={'size':10}, numpoints=1)


    ax2 = plt.subplot(2,1,2)
    ax2.plot(sorted(pl_tfs.values()),range(1,len(pl_tfs)+1),'--s',color=plcolor,label='PathLinker TRs',lw=1)
    ax2.plot(sorted(pr_tfs.values()),range(1,len(pr_tfs)+1),'--s',color=prcolor,label='RWR TRs',lw=1)
    ax2.set_title('TRs in the %s reconstruction' % (pathway))
    ax2.set_xlabel('Ranked Interaction',size=12)
    ax2.set_ylabel('# TRs',size=12)
    ax2.legend(loc='lower right', prop={'size':10}, numpoints=1)

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
    pr_tfs,pr_receptors = readPageRankEdges(tfs,receptors,opts.indir,opts.pathway,opts.kegg)
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
