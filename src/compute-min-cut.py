import glob
from utilsPoirel import *
import networkx as nx
from networkx import DiGraph 
from optparse import OptionParser

# mincut variables
SOURCE='source'
SINK='sink'
MINCUTTHRES = 0

## Computes the minimum cut, after connecting all the sources to a dummy
## node and all the sinks to another dummy node.
## Input: list of tuples (edges), set of nodes (sources), set of nodes (sinks)
def getMinCut(nodefile, ppidir, mapper, out):
    pathway = nodefile.split('/')[-1]
    pathway = pathway.replace('-nodes.txt','')
    
    # read PPI nodes
    ppiedges = readColumns('%s/%s-interactome.txt' % (ppidir,pathway),1,2)
    ppinodes = set([u for u,v in ppiedges]).union(set([v for u,v in ppiedges]))

    # read nodes
    nodes = readColumns(nodefile,1,2)
    receptors = set([n for n,ntype in nodes if ntype =='receptor'])
    tfs = set([n for n,ntype in nodes if ntype =='tf'])

    edgefile = nodefile.replace('nodes','edges')
    edges = set(readColumns(edgefile,1,2))

    ## get dead-ends
    tails = set([u for u,v in edges])
    deadends = set([v for u,v in edges if v not in tails])
    nontfdeadends = deadends.difference(tfs).difference(receptors)
    if len(nontfdeadends)>0:
        print 'WARNING: pathway %s (%s) contains non-TF dead ends' % (pathway,mapper.get(pathway,pathway)) 
        print ' --> ', edgefile
        print ' --> ', nontfdeadends

    ## get mincut
    if len(receptors)>0 and len(tfs)>0:
        # make super source and super sink
        for s in receptors:
            edges.add((SOURCE,s))
        for t in tfs:
            edges.add((t,SINK))

        # make networkX digraph
        G = DiGraph()
        for t,h in edges:
            if t == SOURCE or h == SINK:
                G.add_edge(t,h,capacity=10000000000.0)
            else:
                G.add_edge(t,h,capacity=1.0)
        mincut = nx.minimum_cut(G,SOURCE,SINK)
    else:
        mincut = -1       

    print '%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d' % \
        (pathway, mapper.get(pathway,pathway), len(nodes), len(edges),\
         mincut[0], len(receptors), len(tfs),\
         len(receptors.intersection(ppinodes)),\
         len(tfs.intersection(ppinodes)), len(nontfdeadends))

    out.write('%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % \
        (pathway, mapper.get(pathway,pathway), len(nodes), len(edges),\
        mincut[0], len(receptors), len(tfs),\
        len(receptors.intersection(ppinodes)),\
        len(tfs.intersection(ppinodes)), len(nontfdeadends)))

    return

#######################################################################
def main(args):
    usage = '''compute-min-cut.py [options]
'''
    parser = OptionParser(usage=usage)

    # General Options
    parser.add_option('','--datadir',type='string',metavar='STR',\
                      help='directory containing curated pathways')
    parser.add_option('','--ppidir',type='string',metavar='STR',\
                      help='location of PPI directory.')
    parser.add_option('-o', '--outfile', type='string', metavar='STR',\
        help='output file. Required.')
    parser.add_option('','--mapfile',type='string',metavar='STR',\
                      help='if specified, prints an additional column')

    # parse the command line arguments
    (opts, args) = parser.parse_args()
    if opts.datadir == None:
        sys.exit('ERROR: data directory required.')
    if opts.ppidir == None:
        sys.exit('ERROR: ppi directory required.')
    if opts.outfile == None:
        sys.exit('ERROR: outfile required.')

    print '\nOPTIONS ARE', opts
    print

    if opts.mapfile == None:
        mapper = {}
    else:
        mapper = readDict(opts.mapfile,1,2)    
    
    nodefiles = glob.glob('%s/*-nodes.txt' % (opts.datadir))
    out = open(opts.outfile,'w')
    out.write('#ID\tName\tnumnodes\tnumedges\tmincut\tnum_receptors\tnum_tfs\tnum_recoverable_receptors\tnum_recoverable_tfs\tnum_nontfs_deadends\n')
    print '#ID\tName\tnumnodes\tnumedges\tmincut\tnum_receptors\tnum_tfs\tnum_recoverable_receptors\tnum_recoverable_tfs\tnum_nontfs_deadends'
    for nodefile in nodefiles:
        getMinCut(nodefile, opts.ppidir, mapper, out)
    out.close()
    print 'Wrote to %s' % (opts.outfile)
        

if __name__=='__main__':
    main(sys.argv)
