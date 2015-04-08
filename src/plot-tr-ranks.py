from optparse import OptionParser
import sys
from utilsPoirel import *
import networkx as nx
import matplotlib 
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

## http://matplotlib.org/examples/pylab_examples/broken_axis.html

#######################################################################
def plotRanks(tfs,pathlinker,pagerank,ipa,outprefix):
    iparanks = [pagerank[k] for k in ipa]
    maxval = max(max(pagerank.values()),max(pathlinker.values()))
    ply = .3
    pry = 0
    ipay = -.1
    # If we were to simply plot pts, we'd lose most of the interesting
    # details due to the outliers. So let's 'break' or 'cut-out' the y-axis
    # into two portions - use the top (ax) for the outliers, and the bottom
    # (ax2) for the details of the majority of our data
    fig,(ax,ax2) = plt.subplots(1,2,sharey=True,figsize=(18,3))

    # plot the same data on both axes
    ax.plot(sorted(pagerank.values()),[pry]*len(pagerank),'s',ms=8)
    ax.plot(sorted(pathlinker.values()),[ply]*len(pathlinker),'s',ms=8)
    ax.plot(sorted(iparanks),[ipay]*len(iparanks),'r*',ms=8)
    for tf in pathlinker.keys():
        ax.plot([pagerank[tf],pathlinker[tf]],[pry+0.01,ply-0.01],'-ok',ms=4)

    ax2.plot(sorted(pagerank.values()),[pry]*len(pagerank),'s',ms=8)
    ax2.plot(sorted(pathlinker.values()),[ply]*len(pathlinker),'s',ms=8)
    ax2.plot(sorted(iparanks),[ipay]*len(iparanks),'r*',ms=8)
    for tf in pathlinker.keys():
        ax2.plot([pagerank[tf],pathlinker[tf]],[pry+0.01,ply-0.01],'-ok',ms=4)
    
    # zoom-in / limit the view to different portions of the data
    ax.set_xlim(0,1300)
    ax.set_ylim(ipay-0.05,ply+0.05)
    ax2.set_xlim(maxval-1000,maxval+50) # outliers only
    ax2.set_ylim(ipay-0.05,ply+0.05)

    # hide the spines between ax and ax2
    ax.spines['right'].set_visible(False)
    ax.yaxis.tick_left()
    ax.yaxis.set_ticks([ipay,pry,ply])
    ax.yaxis.set_ticklabels(['IPA','PageRank','PathLinker'])
    ax.tick_params(labelright='off') # don't put tick labels at the top
    ax2.spines['left'].set_visible(False)

    # This looks pretty good, and was fairly painless, but you can get that
    # cut-out diagonal lines look with just a bit more work. The important
    # thing to know here is that in axes coordinates, which are always
    # between 0-1, spine endpoints are at these locations (0,0), (0,1),
    # (1,0), and (1,1).  Thus, we just need to put the diagonals in the
    # appropriate corners of each of our axes, and so long as we use the
    # right transform and disable clipping.

    d = .015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1-d,1+d),(1-d,1+d), **kwargs)      # top-left diagonal
    ax.plot((1-d,1+d),(-d,+d), **kwargs)    # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d),(1-d,1+d), **kwargs)   # bottom-left diagonal
    ax2.plot((-d,+d),(-d,+d), **kwargs) # bottom-right diagonal

    #plt.tight_layout()
    plt.savefig(outprefix+'.png')
    print 'Wrote to '+outprefix+'.png'
    return

#######################################################################
def main(args):
    usage = '''plot-tr-ranks.py [options]
'''
    parser = OptionParser(usage=usage)

    # General Options
    parser.add_option('','--nodefile',type='string',metavar='STR',\
                      help='Node file with "tf" labeled TRs.')
    parser.add_option('-o', '--outprefix', type='string', metavar='STR',\
        help='A string to prepend to all output files. Required.')
    parser.add_option('','--indir',type='string',metavar='STR',\
                      help='Input directory.')
    parser.add_option('', '--pdf', action='store_true', default=False,\
        help='Output images in PDF as well as PNG.')

    # parse the command line arguments
    (opts, args) = parser.parse_args()
    if opts.indir == None:
        sys.exit('ERROR: input directory required.')
    if opts.outprefix == None:
        sys.exit('ERROR: output directory required.')
    if opts.nodefile == None:
        sys.exit('ERROR: must specify node file.')

    print '\nOPTIONS ARE', opts
    
    # read nodes
    tfs = {u:n for u,t,n in readColumns(opts.nodefile,1,2,3) if t == 'tf'}
    receptors = {u:n for u,t,n in readColumns(opts.nodefile,1,2,3) if t == 'receptor'}
    
    # read pathlinker edges
    infile = '%s/pathlinker/Wnt-k_20000-ranked-edges.txt' % (opts.indir)
    lines = readColumns(infile,1,2,3)
    pathlinker = {}
    for u,v,ksp in lines:
        if u in tfs and u not in pathlinker:
            pathlinker[u] = int(ksp)
        if v in tfs and v not in pathlinker:
            pathlinker[v] = int(ksp)
    
    # read pagerank edges
    infile = '%s/pagerank/Wnt-q_0.50-node-pagerank.txt' % (opts.indir)
    lines = readColumns(infile,1,2)
    i=1
    pagerank = {}
    for n,p in sorted(lines,key=lambda x:float(x[1]),reverse=True):
        if n in tfs:
            pagerank[n] = i
        i+=1
    
    # read IPA edges
    infile = '%s/ipa/Wnt-nmax35.out' % (opts.indir)
    G = nx.Graph()
    G.add_edges_from(readColumns(infile,1,2))
    G.add_node('SOURCE')
    for r in receptors:
        G.add_edge('SOURCE',r)
    ipa = set([tf for tf in tfs if nx.has_path(G,'SOURCE',tf)])

    for tf,val in sorted(pagerank.items(), key=lambda x:x[1]):
        print tf,tfs[tf],val,pathlinker.get(tf,'NA'),tf in ipa

    plotRanks(tfs,pathlinker,pagerank,ipa,opts.outprefix)
    

    
if __name__=='__main__':
    main(sys.argv)
