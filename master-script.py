#!/usr/bin/python

## Anna Ritz, March 2015
## Master Script for running PathLinker for Nat. Methods submisssion.

from optparse import OptionParser,OptionGroup
from utilsPoirel import *
import sys
import os
import os.path
import subprocess



## INTERACTOME VARIABLES (will be filled in in main())
PPIVERSION = ''
ALLOWEDVERSIONS = [#'2013linker', # Chris's version w/ vinayagam, etc.\
                   #'2015pathlinker', # new version, reg weighting\
                   #'pathlinker-old-kegg-buggy', # debug version\
                   #'pathlinker-old-kegg-fixed', # debug version\
                   #'pathlinker-old-kegg-new-parse', # debug version\
                   #'pathlinker-old-kegg-buggy-oldannotations', # debug version\
                   #'pathlinker-signaling',  # only weighted by sig-trans     \
                   #'pathlinker-signaling-children',  # weighting\
                   'pathlinker-signaling-children-reg',  # children and reg \
               ]
DATADIR = '/data/annaritz/datasets/svn-data/'
PPIDIR = '/data/annaritz/projects/2015-03-pathlinker/data/pathway-specific-interactomes/'

## MAPPING VARIABLES
MAPPINGFILE = '%s/namespace-mappers/human-gene-map.txt' % (DATADIR)
MAPFROMCOL = 6
MAPTOCOL = 1

## PATHWAY DIRECTORIES
NETPATHDIR='/data/annaritz/datasets/svn-data/interactions/netpath/pathways/'
KEGGDIR='/data/annaritz/datasets/svn-data/interactions/kegg/2015-03-23/hsa/edge-files/'

## VARYPARAMS
VARYPARAMS = {'q': [0.1, 0.25, 0.5, 0.75, 0.9],
              'gamma':[10, 15, 20, 25, 30],
              'omega':[0, 0.01, 0.1],
              'prize':[1, 3, 5, 7, 9],
              'alpha':[0, 0.1, 0.25, 0.4, 0.5],
              'nmax':[5,10,15, 25, 35, 50, 75, 100, 200, 500],
          }
    
def main(args):
    global PPIVERSION, PPIDIR
    usage = 'master-script.py [options]\n'
    parser = OptionParser(usage=usage)

    parser.add_option('','--forcealg',action='store_true',default=False,\
                          help='Run algorithms even if they will overwrite existing files.  Default is that algorithms are not run if output files exist.')
    parser.add_option('','--forceprecrec',action='store_true',default=False,\
                      help='Run computations of precision and recall, even if they will overwrite existing files. Default is that precision and recall are not run if output files exist.')
    parser.add_option('','--printonly',action='store_true',default=False,\
                          help='Print the commands to stdout, but do not execute.')

    group = OptionGroup(parser,'Datasets')
    group.add_option('','--ppiversion',type='string',default='2015pathlinker',\
                     help='Version of the PPI to run.  Options are %s. Default is "2015pathlinker."' % (', '.join(ALLOWEDVERSIONS)))
    group.add_option('','--weightedppi',action='store_true',default=False,\
                         help='Run with weighted PPI.')
    group.add_option('','--onlynetpathwnt',action='store_true',default=False,\
                         help='Only run NetPath Wnt.')
    group.add_option('','--netpath',action='store_true',default=False,\
                         help='Run with NetPath inputs.')
    group.add_option('','--kegg',action='store_true',default=False,\
                         help='Run with KEGG inputs.')
    group.add_option('','--wntforexperiments',action='store_true',default=False,\
                     help='Run special wnt that includes FZD4/FZD6 receptors, for analyzing via networks.')
    group.add_option('','--dbcompare',action='store_true',default=False,\
                     help='Run the 6 KEGG and NetPath pathways in common.')
    group.add_option('','--aggunion',action='store_true',default=False,\
                         help='Run aggregate union special runs')
    group.add_option('','--missingnpkegg',action='store_true',default=False,\
                         help='Run Linker for interactomes with missing evidence.')
    parser.add_option_group(group)

    group = OptionGroup(parser,'Algorithms')
    group.add_option('','--pathlinker',action='store_true',default=False,\
                          help='Run PathLinker (KSP) on input files.')
    group.add_option('','--shortestpaths',action='store_true',default=False,\
                     help='Compute shortest paths from each receptor to each TF, incl. ties (reviewer comment)')
    group.add_option('','--inducedsubgraph',action='store_true',default=False,\
                     help='Re-rank Linker nodes by taking induced subgraph.')
    group.add_option('','--rerank',action='store_true',default=False,\
                     help='Re-rank Linker runs by unique nodes and edges.')
    group.add_option('','--pagerank',action='store_true',default=False,\
                     help='Run PageRank on input files.')
    group.add_option('','--eqed',action='store_true',default=False,\
                     help='Run eQED on input files.')
    group.add_option('','--pcsf',action='store_true',default=False,\
                         help='Run prize collecting steiner forest.')
    group.add_option('','--anat',action='store_true',default=False,\
                         help='Run ANAT.')
    group.add_option('','--responsenet',action='store_true',default=False,\
                         help='Run ResponseNet')
    group.add_option('','--degree',action='store_true',default=False,\
                         help='Compute degree of PPI')
    group.add_option('','--ipa',action='store_true',default=False,\
                         help='Compute IPA\'s Network Generation Algorithm on the undirected PPI.')

    parser.add_option_group(group)

    group = OptionGroup(parser,'Visualizations')
    group.add_option('','--computeprecrec',action='store_true',default=False,\
                         help='Compute precision and recall curves; write to file.')
    group.add_option('','--precrecviz',action='store_true',default=False,\
                         help='Display precision recall curves.')
    group.add_option('','--ranktfs',action='store_true',default=False,\
                     help='Plot TF ranks for PageRank, PathLinker, and IPA.')
    group.add_option('','--bullseye',action='store_true',default=False,\
                         help='Make bullseye plots.')
    group.add_option('','--falsenegs',action='store_true',default=False,\
                         help='Make false negative plots.')
    group.add_option('','--dbcompareviz',action='store_true',default=False,\
                         help='Post Wnt KEGG vs. NetPath DB comparison to GraphSpace.')
    group.add_option('','--varyparamsviz',action='store_true',default=False,\
                         help='run precrec/precrecviz with varying parameters.')
    group.add_option('','--venn',action='store_true',default=False,\
                         help='plot venn diagrams. Only with dbpcompare.')
    group.add_option('','--graphspace',action='store_true',default=False,\
                         help='post images to graphspace. Only with wntforexperimentsviz and dbcompare.')
    group.add_option('','--paper',action='store_true',default=False,\
                         help='Make plots for paper.')
    group.add_option('','--weightviz',action='store_true',default=False,\
                         help='Post weighted vs. unweighted graph to GraphSpace')
    parser.add_option_group(group)

    group = OptionGroup(parser,'Optional Arguments')
    group.add_option('','--varyparams',action='store_true',default=False,\
                     help='Run the algorithms that vary parameters with different parameter values.  These algorithms are PageRank, ResponseNet, PCSF, ANAT, and IPA.')
    group.add_option('','--q',type='float',default='0.5',metavar='FLOAT',\
                          help='probability of teleporting in the random walk. Default is 0.5.')
    group.add_option('','--k',type='int',default=20000,metavar='INT',\
                          help='# of Shortest Paths; only useful when LINKER (-l) or VIZ (-v) is specified. Default is 20000.')
    group.add_option('','--topk',type='int',default=200,metavar='INT',\
                         help='# of paths to visualize (Default is 200).')
    group.add_option('','--inputcurrent',type='int',default=10000,metavar='INT',\
                     help='Input current for eQED. Default = 10000')
    group.add_option('','--gamma',type='int',default=20,metavar='INT',\
                         help='Gamma value for ResponseNet. Default is 20.')
    group.add_option('','--prize',type='float',default=5,metavar='FLOAT',\
                         help='Prize for PCSF. Default is 5.')
    group.add_option('','--omega',type='float',default='0.01',metavar='FLOAT',\
                     help='Cost of adding additional trees to PCSF forest. Default is 0.01.')
    group.add_option('','--alpha',type='float',default=0.0,metavar='FLOAT',\
                     help='Tradeoff between shortest paths and steiner trees in ANAT. Default is 0.00.')
    group.add_option('','--nmax',type='int',default=10,metavar='INT',\
                     help='Size for IPA\'s network generation algorithm. Default is 10')
    group.add_option('','--subsamplefps',type='int',default=50,metavar='INT',\
                         help='Subsample so the number of FPs is this many times the # of TPs for each pathway.  If this value is -1, then all negatives are used.')
    group.add_option('','--bullseyerecall',type='float',default=0.5,metavar='FLOAT',\
			 help='Recall for bullseye plots.  Default is 0.5.')
    parser.add_option_group(group)
    
    # parse the command line arguments
    (opts, args) = parser.parse_args()

    if opts.ppiversion not in ALLOWEDVERSIONS:
        sys.exit('ERROR: --ppiversion must be one of %s. You have %s. ' % (','.join(ALLOWEDVERSIONS),opts.ppiversion))
    PPIVERSION = opts.ppiversion

    ## BACKGROUND INTERACTOME
    ## set PPI, determined by whether the --weightedppi argument is specified.
    if opts.weightedppi:
        ## ppifile contains the probabilities from weighting scheme
        ppifile = '/data/annaritz/datasets/svn-data/interactomes/human/%s-weighted.txt' % (PPIVERSION)
        resultprefix = 'results/%s/weighted/' % (PPIVERSION)
        PPIDIR = '%s/%s/weighted/' % (PPIDIR,PPIVERSION)
    else:
        ## ppifile has weights of 1 for all edges.
        ppifile = '/data/annaritz/datasets/svn-data/interactomes/human/%s.txt' % (PPIVERSION)
        resultprefix = 'results/%s/unweighted/' % (PPIVERSION)
        PPIDIR = '%s/%s/unweighted/' % (PPIDIR,PPIVERSION)
    checkDir(resultprefix)
    checkDir(PPIDIR)

    ## PATHWAY-SPECIFIC INTERACTOMES
    generatePathwaySpecificInteractomes(ppifile)

    # DATASETS
    ## pathways is a set() of (pathwayname,resultdir,datadir) tuples.
    pathways = set()
    kegg2netpath = {} # will be populated if kegg pathways are specified.
    if opts.netpath or opts.onlynetpathwnt or opts.dbcompare:
        # Get NetPath pathway names.  
        pathwaynames = getNetPathPathways(opts.onlynetpathwnt,opts.dbcompare)
        # create result directory 
        resultdir = '%s/netpath/' %(resultprefix)
        checkDir(resultdir) 
        # data directory is netpath directory.
        datadir = NETPATHDIR
        # ppi file 
        ppidir = '%s/netpath/' % (PPIDIR)
        # add (pathwayname,resultdir,datadir) tuples
        pathways.update([(p,resultdir,datadir,ppidir) for p in pathwaynames])
        print 'Running %d NetPath pathways' % (len(pathwaynames))
    if opts.wntforexperiments:
        resultdir = '%s/wnt-all-receptors/' % (resultprefix)
        checkDir(resultdir)
        datadir = 'data/wnt-all-receptors/'
        ppidir = '%s/wnt-all-receptors/' % (PPIDIR)
        pathways.update([('Wnt',resultdir,datadir,ppidir)])
    if opts.kegg or opts.dbcompare:
        # get KEGG pathway names
        pathwaynames,kegg2netpath = getKEGGPathways() 
        # create result directory
        resultdir = '%s/kegg/' % (resultprefix)
        checkDir(resultdir) 
        # data directory is kegg directory
        datadir = KEGGDIR
        ppidir = '%s/kegg/'% (PPIDIR)
        pathways.update([(p,resultdir,datadir,ppidir) for p in pathwaynames])
        print 'Running %d KEGG pathways' % (len(pathwaynames))
    if len(pathways)==0:
        print 'WARNING: no datasets specified. This is OK when visualizing output or running special runs.'
    else:
        print '%d total pathways considered.' % (len(pathways))

    ## ALGORITHMS ##

    ## PathLinker ##
    if opts.pathlinker and not opts.missingnpkegg:
        print 'Running PathLinker:'
        for (pathway,resultdir,datadir,ppidir) in pathways:
            runPathLinker(pathway,resultdir,datadir,ppidir,opts.k,opts.forcealg,opts.printonly)
        print 'Done Running PathLinker\n'

    ## Shortest Paths ##
    if opts.shortestpaths:
        print 'Running Shortest Paths'
        for (pathway,resultdir,datadir,ppidir) in pathways:
            runShortestPaths(pathway,resultdir,datadir,ppidir,opts.forcealg,opts.printonly)
        print 'Done Running Shortest Paths\n'

    ## Induced Subgraph ##
    if opts.inducedsubgraph:
        print 'Getting induced subgraph from PathLinker'
        for (pathway,resultdir,datadir,ppidir) in pathways:
            runInducedSubgraph(pathway,resultdir,datadir,ppidir,opts.forcealg,opts.printonly)
        print 'Done getting induced subgraph from PathLinker.\n'

    ## Reranking PathLinker ##
    if opts.rerank:
        print 'Reranking PathLinker'
        for (pathway,resultdir,datadir) in pathways:
            rerankPathLinker(pathway,resultdir,datadir,opts.forcealg,opts.printonly)
        print 'Done reranking PathLinker.\n'

    ## PageRank ##    
    if opts.pagerank:
        print 'Running PageRank'
        if not opts.varyparams: # just run with opts.q value
            for (pathway,resultdir,datadir,ppidir) in pathways:
                runPageRank(pathway,resultdir,datadir,ppidir,opts.q,opts.forcealg,opts.printonly)
        else: # vary opts.q value
            for varyq in VARYPARAMS['q']:
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    runPageRank(pathway,resultdir,datadir,ppidir,varyq,opts.forcealg,opts.printonly)
        print 'Done runing PageRank.\n'

    ## EQED ##
    if opts.eqed:
        print 'Running eQED:'
        for (pathway,resultdir,datadir,ppidir) in pathways:
            runEQED(pathway,resultdir,datadir,ppidir,opts.inputcurrent,opts.forcealg,opts.printonly)
        print 'Done running eQED\n'

    ## ResponseNet ##
    if opts.responsenet:
        print 'Running ResponseNet:'
        if not opts.varyparams: # just run with opts.gamma value
            for (pathway,resultdir,datadir,ppidir) in pathways:
                runResponseNet(pathway,resultdir,datadir,ppidir,opts.gamma,opts.forcealg,opts.printonly)
        else: # vary opts.gamma value
            for varygamma in VARYPARAMS['gamma']:
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    runResponseNet(pathway,resultdir,datadir,ppidir,varygamma,opts.forcealg,opts.printonly)
        print 'Done running ResponseNet\n'

    ## PCSF ##
    if opts.pcsf:
        print 'Running PCSF'
        if not opts.varyparams: # just run with opts.prize and opts.omega values
            for (pathway,resultdir,datadir,ppidir) in pathways:
                runPCSF(pathway,resultdir,datadir,ppidir,opts.prize,opts.omega,opts.forcealg,opts.printonly)
        else: # vary prize and omega
            for varyprize in VARYPARAMS['prize']:
                for varyomega in VARYPARAMS['omega']:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        runPCSF(pathway,resultdir,datadir,ppidir,varyprize,varyomega,opts.forcealg,opts.printonly)
        print 'Done running PCSF\n'

    ## ANAT  ##
    if opts.anat:
        print 'Running ANAT:'
        if not opts.varyparams: # just run with opts.alpha value
            for (pathway,resultdir,datadir,ppidir) in pathways:
                runANAT(pathway,resultdir,datadir,ppidir,opts.alpha,opts.forcealg,opts.printonly)
        else: # vary alpha
            for varyalpha in VARYPARAMS['alpha']:
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    runANAT(pathway,resultdir,datadir,ppidir,varyalpha,opts.forcealg,opts.printonly)
        print 'Done Running ANAT\n'

    ## DEGREE ##
    if opts.degree:
        print 'NOT UPDATED TO HANDLE PATHWAY-SPECIFIC INTERACTOMES!! skipping.'
        # print 'Computing Weighted Degree of PPIFILE:'
        # ## resultprefix is either "results-weighted" or "results-unweighted"
        # outprefix = '%s/ppidegree/weighted-degree' % (resultprefix)
        # checkDir(outprefix)

        # if opts.forcealg or not os.path.isfile('%s-edges.txt' % (outprefix)):
        #     script = '/data/annaritz/signaling/2014-06-linker/src/compute-ppi-degrees.py'
        #     cmd = 'python %s %s %s %s' % (script,PPIFILE,outprefix,MAPPINGFILE)
        #     print cmd
        #     if not opts.printonly:
        #         subprocess.check_call(cmd.split())
        # else:
        #     print 'Skipping weighting PPI: %s already exists.' % ('%s-edges.txt' % (outprefix))
        # print 'Done computing Degree of PPIFILE\n'

    ## IPA ##
    if opts.ipa:
        print 'Running IPA'
        if not opts.varyparams: # just run with opts.nmax value
            for (pathway,resultdir,datadir,ppidir) in pathways:
                runIPA(pathway,resultdir,datadir,ppidir,opts.nmax,opts.forcealg,opts.printonly)
        else: # vary nmax
            for varynmax in VARYPARAMS['nmax']:
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    runIPA(pathway,resultdir,datadir,ppidir,varynmax,opts.forcealg,opts.printonly)
        print 'Done IPA\n'

    ## VIZ SCRIPTS ##

    ## write precision/recall
    if opts.computeprecrec:
        print 'Write Precision Recall'
        for negtype in ['none','adjacent']:

            ## get sample output prefix (all methods use the same sampled positives/negatives)
            sampledir = '%s/samples-exclude-%s' % (resultprefix,negtype)
            checkDir(sampledir)

            if opts.wntforexperiments:
                wntsampledir = '%s/wntforexperiments-samples-exclude-%s' % (resultprefix,negtype)
                checkDir(wntsampledir)
            
            ## PATHLINKER ##
            if opts.pathlinker:
                sortcol = 3 # ksp
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/pathlinker/%s-k_%s-ranked-edges.txt' % (resultdir,pathway,opts.k)
                    outdir = '%s/precision-recall/pathlinker/' % (resultdir)
                    if 'wnt-all-receptors' in resultdir:
                        sampleoutprefix = '%s/%s' % (wntsampledir,pathway)
                    else:
                        sampleoutprefix = '%s/%s' % (sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly)
                if opts.netpath:
                    ## compute aggregate for NetPath
                    inputdir = '%s/netpath/precision-recall/pathlinker/' % (resultprefix)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,opts.forceprecrec,opts.printonly)
                    
            ## Shortest Paths ##
            if opts.shortestpaths:
                sortcol = None # no sorting; entire file is shortest paths subgraph.
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/shortestpaths/%s-shortest-paths.txt' % (resultdir,pathway)
                    outdir = '%s/precision-recall/shortestpaths/' % (resultdir)
                    if 'wnt-all-receptors' in resultdir:
                        sampleoutprefix = '%s/%s' % (wntsampledir,pathway)
                    else:
                        sampleoutprefix = '%s/%s' % (sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly)
                if opts.netpath:
                    ## compute aggregate for NetPath
                    inputdir = '%s/netpath/precision-recall/shortestpaths/' % (resultprefix)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,opts.forceprecrec,opts.printonly)
            
            ## Induced Subgraph ##
            if opts.inducedsubgraph:
                sortcol = 3 # ksp
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/inducedsubgraph/%s-induced-subgraph.txt' % (resultdir,pathway)
                    outdir = '%s/precision-recall/inducedsubgraph/' % (resultdir)
                    if 'wnt-all-receptors' in resultdir:
                        sampleoutprefix = '%s/%s' % (wntsampledir,pathway)
                    else:
                        sampleoutprefix = '%s/%s' % (sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly)
                if opts.netpath:
                    ## compute aggregate for NetPath
                    inputdir = '%s/netpath/precision-recall/inducedsubgraph/' % (resultprefix)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,opts.forceprecrec,opts.printonly)

            ## PageRank ##
            if opts.pagerank:
                # both sort columns are in decreasing order.
                edgesortcol = 3 # edge_flux
                nodesortcol = 2 # visitation_probability
                if not opts.varyparams: # just run with opts.q value
                    params = ['q_%.2f' % (opts.q)]
                else: # vary opts.q value
                    params = ['q_%.2f' % p for p in VARYPARAMS['q']]
                for param in params:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        edgefile = '%s/pagerank/%s-%s-edge-fluxes.txt' % (resultdir,pathway,param)
                        nodefile = '%s/pagerank/%s-%s-node-pagerank.txt' % (resultdir,pathway,param)
                        outdir = '%s/precision-recall/pagerank/' % (resultdir)
                        if 'wnt-all-receptors' in resultdir:
                            sampleoutprefix = '%s/%s' % (wntsampledir,pathway)
                        else:
                            sampleoutprefix = '%s/%s' % (sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,nodefile=nodefile,nodesortcol=nodesortcol,\
                                               descending=True,param=param)
                    if opts.netpath:
                        ## compute aggregate for NetPath
                        inputdir = '%s/netpath/precision-recall/pagerank/' % (resultprefix)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,descending=True,
                                                        param=param)

            ## eQED ##
            if opts.eqed:
                # both sort columns are in decreasing order.
                edgesortcol = 4 # abs(current)
                nodesortcol = 3 # positive input current
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/eqed/%s-eqed-edges.out' % (resultdir,pathway)
                    nodefile = '%s/eqed/%s-eqed-nodes.out' % (resultdir,pathway)
                    outdir = '%s/precision-recall/eqed/' % (resultdir)
                    if 'wnt-all-receptors' in resultdir:
                        sampleoutprefix = '%s/%s' % (wntsampledir,pathway)
                    else:
                        sampleoutprefix = '%s/%s' % (sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                           opts.printonly,nodefile=nodefile,nodesortcol=nodesortcol,descending=True)
                if opts.netpath:
                    ## compute aggregate for NetPath
                    inputdir = '%s/netpath/precision-recall/eqed/' % (resultprefix)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,opts.forceprecrec,opts.printonly,descending=True)

            ## ResponseNet ##
            if opts.responsenet:
                ## ANNA CHANGE: take union of edges/nodes with positive flow as a subgraph.
                # both sort columns are in decreasing order.
                edgesortcol = None 
                if not opts.varyparams: # just run with opts.gamma value
                    params = ['gamma_%d' % (opts.gamma)]
                else: # vary opts.gamma value
                    params = ['gamma_%d' % p for p in VARYPARAMS['gamma']]
                for param in params:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        edgefile = '%s/reponsenet/%s-%s_responsenet-edges.out' % (resultdir,pathway,param)
                        outdir = '%s/precision-recall/responsenet/' % (resultdir)
                        if 'wnt-all-receptors' in resultdir:
                            sampleoutprefix = '%s/%s' % (wntsampledir,pathway)
                        else:
                            sampleoutprefix = '%s/%s' % (sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,param=param)
                    if opts.netpath:
                        ## compute aggregate for NetPath
                        inputdir = '%s/netpath/precision-recall/responsenet/' % (resultprefix)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param)
            ## PCSF ##
            if opts.pcsf:
                sortcol = None # no sorting; entire file is shortest paths subgraph.
                if not opts.varyparams: # just run with opts.prize and opts.omega values
                    params = ['prize%d-omega%.2f' % (opts.prize,opts.omega)]
                else: # vary prize and omega
                    params = []
                    for varyprize in VARYPARAMS['prize']:
                        for varyomega in VARYPARAMS['omega']:
                            params.append('prize%d-omega%.2f' % (varyprize,varyomega))
                for param in params:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        edgefile = '%s/pcsf/%s-%s_PCSF-edges.out' % (resultdir,pathway,param)
                        outdir = '%s/precision-recall/pcsf/' % (resultdir)
                        if 'wnt-all-receptors' in resultdir:
                            sampleoutprefix = '%s/%s' % (wntsampledir,pathway)
                        else:
                            sampleoutprefix = '%s/%s' % (sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,param=param)
                    if opts.netpath:
                        ## compute aggregate for NetPath
                        inputdir = '%s/netpath/precision-recall/pcsf/' % (resultprefix)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param)
            ## ANAT ##
            if opts.anat:
                sortcol = None # no sorting; entire file is shortest paths subgraph.
                if not opts.varyparams: # just run with opts.alpha value
                    params = ['alpha%.2f' % (opts.alpha)]
                else: # vary opts.alpha value
                    params = ['alpha%.2f' % p for p in VARYPARAMS['alpha']]
                for param in params:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        edgefile = '%s/anat/%s-%s-edges.out' % (resultdir,pathway,param)
                        outdir = '%s/precision-recall/anat/' % (resultdir)
                        if 'wnt-all-receptors' in resultdir:
                            sampleoutprefix = '%s/%s' % (wntsampledir,pathway)
                        else:
                            sampleoutprefix = '%s/%s' % (sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,param=param)
                    if opts.netpath:
                        ## compute aggregate for NetPath
                        inputdir = '%s/netpath/precision-recall/anat/' % (resultprefix)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param)
            ## DEGREE ##

            ## IPA ##
            if opts.ipa:
                sortcol = None # no sort column; take entire file
                if not opts.varyparams: # just run with opts.nmax value
                    params = ['nmax%d' % (opts.nmax)]
                else: # vary opts.nmax value
                    params = ['nmax%d' % p for p in VARYPARAMS['nmax']]
                for param in params:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        edgefile = '%s/ipa/%s-%s.out' % (resultdir,pathway,param)
                        outdir = '%s/precision-recall/ipa/' % (resultdir)
                        if 'wnt-all-receptors' in resultdir:
                            sampleoutprefix = '%s/%s' % (wntsampledir,pathway)
                        else:
                            sampleoutprefix = '%s/%s' % (sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,descending=True,param=param)
                    if opts.netpath:
                        ## compute aggregate for NetPath
                        inputdir = '%s/netpath/precision-recall/ipa/' % (resultprefix)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        descending=True,param=param)
    if opts.precrecviz:
        print 'Plot Precision Recall'
        algcmd = ''
        if opts.pathlinker:
            algcmd += ' --alg pathlinker'
        if opts.pagerank:
            algcmd += ' --alg pagerank'
        if opts.shortestpaths:
            algcmd += ' --alg shortestpaths'
        if opts.inducedsubgraph:
            algcmd += ' --alg inducedsubgraph'
        if opts.eqed:
            algcmd += ' --alg eqed'
        if opts.responsenet:
            algcmd += ' --alg responsenet'
        if opts.pcsf:
            algcmd += ' --alg pcsf'
        if opts.anat: 
            algcmd += ' --alg anat'
        if opts.ipa:
            algcmd += ' --alg ipa'
        indir = '%s/netpath/precision-recall/' % (resultprefix)
        for (pathway,resultdir,datadir,ppidir) in pathways:
            if opts.wntforexperiments and 'wnt-all-receptors' in resultdir:
                outprefix = 'viz/precision-recall/%s-all-receptors' % (pathway)
            else:
                outprefix = 'viz/precision-recall/%s' % (pathway)
            cmd = 'python src/plot-precision-recall.py --indir %s --outprefix %s --pathway %s %s' % (indir,outprefix,pathway,algcmd)
            if opts.varyparams:
                if pathway != 'Wnt':
                    print 'WARNING: not plotting non-Wnt pathway %s with varying parameters' % (pathway)
                else:
                    cmd += '  --varyparams'
                    print cmd
                    if not opts.printonly:
                        subprocess.check_call(cmd.split())
            else:
                print cmd
                if not opts.printonly:
                    subprocess.check_call(cmd.split())
                
        if opts.netpath: # plot aggregate
            outprefix = 'viz/precision-recall/aggregate'
            cmd = 'python src/plot-precision-recall.py --indir %s --outprefix %s --pathway aggregate %s' % \
                  (indir,outprefix,algcmd)
            if opts.varyparams:
                cmd += '  --varyparams'

            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())

    if opts.graphspace:
        print 'Posting to GraphSpace...'
        
        if opts.pathlinker:
            ## pathlinkerfile is the ranked edges.
            infile = '%s/wnt-all-receptors/pathlinker/Wnt-k_%d-ranked-edges.txt' % (resultprefix,opts.k)
            gsid = 'Wnt-pathlinker-top%dpaths' % (opts.topk)
            postWntReconstructionsToGraphSpace(infile,opts.topk,gsid,opts.printonly,increase=True)

        if opts.pagerank:
            ## threshold is set to 150 edges and 0.1429 recall.
            thres = 0.001677
            infile = '%s/wnt-all-receptors/pagerank/Wnt-q_0.50-edge-fluxes.txt' % (resultprefix)
            gsid = 'Wnt-pagerank-thres0.001677'
            postWntReconstructionsToGraphSpace(infile,thres,gsid,opts.printonly,decrease=True)

        if opts.anat:
            infile = '%s/wnt-all-receptors/anat/Wnt-alpha0.00-edges.out' % (resultprefix)
            gsid = 'Wnt-anat-alpha0.00'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly)

        if opts.shortestpaths:
            infile = '%s/wnt-all-receptors/shortestpaths/Wnt-shortest-paths.txt' % (resultprefix)
            gsid = 'Wnt-shortest-paths'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly)

        if opts.responsenet:
            infile = '%s/wnt-all-receptors/reponsenet/Wnt-gamma_20_responsenet-edges.out' % (resultprefix)
            gsid = 'Wnt-responsenet-gamma20'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly)

        if opts.pcsf:
            infile = '%s/wnt-all-receptors/pcsf/Wnt-prize5-omega0.01_PCSF-edges.out' % (resultprefix)
            gsid = 'Wnt-pcsf-prize5-omega0.01'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly)

        if opts.ipa:
            infile = '%s/wnt-all-receptors/ipa/Wnt-nmax10.out' % (resultprefix)
            gsid = 'Wnt-ipa-nmax10'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly,undirected=True)

    if opts.ranktfs:
        outprefix = 'viz/misc/tr-ranks'
        indir = '%s/netpath/' % (resultprefix)
        nodefile = '/data/annaritz/datasets/svn-data/interactions/netpath/pathways/Wnt-nodes.txt' 
        cmd = 'python src/plot-tr-ranks.py --nodefile %s --o %s --indir %s' % (nodefile,outprefix,indir)
        print cmd
        if not opts.printonly:
            subprocess.check_call(cmd.split())

    if opts.bullseye:
        print 'Bullseye Plots'# at a recall of %.2f' % (opts.bullseyerecall)
        
        ## Do this ONCE: python ../2014-06-linker/src/shortest-paths-for-bullseye.py
        print 'Shortest Paths from nodes to pathways are already computed. See comment in master-script.py.'

        #cmd = 'python ../2014-06-linker/src/plot-bullseye.py -o viz/bullseye/aggregate- --prefix viz/precrecfiles-sample-once-per-pathway/precrec-exclude_none --pathway aggregate --pdf --edges --nodes -r %.2f --alg PRflux+KSP --alg PRflux --alg KSP' % (opts.bullseyerecall)
        #print cmd
        #if not opts.printonly:
        #    subprocess.check_call(cmd.split())

        #cmd = 'python ../2014-06-linker/src/plot-bullseye.py -o viz/bullseye/Wnt- --prefix viz/precrecfiles-sample-once-per-pathway/precrec-exclude_none --pathway Wnt --pdf --edges --nodes -r %.2f --alg PRflux+KSP --alg PRflux --alg KSP' % (opts.bullseyerecall)
        #print cmd
        #if not opts.printonly:
        #    subprocess.check_call(cmd.split())


        cmd = 'python ../2014-06-linker/src/bullseye-to-bar.py -o viz/bullseye-to-bar/aggregate- --prefix viz/precrecfiles-sample-once-per-pathway/precrec-exclude_none --pathway aggregate --pdf --edges --nodes -r 0.30 -r 0.60 --alg PRflux+KSP --alg PRflux --alg KSP --alg NG'
        print cmd
        if not opts.printonly:
            subprocess.check_call(cmd.split())

        cmd = 'python ../2014-06-linker/src/bullseye-to-bar.py -o viz/bullseye-to-bar/Wnt- --prefix viz/precrecfiles-sample-once-per-pathway/precrec-exclude_none --pathway Wnt --pdf --edges --nodes -r 0.30 -r 0.60 --alg PRflux+KSP --alg PRflux --alg KSP --alg NG'
        print cmd
        if not opts.printonly:
            subprocess.check_call(cmd.split())


    if opts.falsenegs:
        print 'False Negative Plots'
        
        ## Do this ONCE: python ../2014-06-linker/src/shortest-paths-for-false-negatives.py
        print 'Shortest Paths from/to nodes are already computed. See comment in master-script.py.'

        #cmd = 'python ../2014-06-linker/src/plot-false-negatives.py -o viz/false-negatives/false-negatives-aggregate- --prefix viz/precrecfiles-sample-once-per-pathway/precrec-exclude_none --pathway aggregate --pdf --nodes --edges --alg PRflux+KSP --alg KSP --alg PRflux -r 0.1 -r .2 -r .3 -r .4 -r .5 -r .6 -r .7 -r .8 -r .9 -r 1.0'
        cmd = 'python ../2014-06-linker/src/plot-false-negatives.py -o viz/false-negatives/false-negatives-aggregate- --prefix viz/precrecfiles-sample-once-per-pathway/precrec-exclude_none --pathway aggregate --pdf --edges --alg PRflux+KSP --alg KSP --alg PRflux -r 0.3 -r 0.6'
        print cmd
        if not opts.printonly:
            subprocess.check_call(cmd.split())

        #cmd = 'python ../2014-06-linker/src/plot-false-negatives.py -o viz/false-negatives/false-negatives-Wnt- --prefix viz/precrecfiles-sample-once-per-pathway/precrec-exclude_none --pathway Wnt --pdf --nodes --edges -r 0.5 --alg PRflux+KSP --alg KSP --alg PRflux -r 0.1 -r .2 -r .3 -r .4 -r .5 -r .6 -r .7 -r .8 -r .9 -r 1.0'
        cmd = 'python ../2014-06-linker/src/plot-false-negatives.py -o viz/false-negatives/false-negatives-Wnt- --prefix viz/precrecfiles-sample-once-per-pathway/precrec-exclude_none --pathway Wnt --pdf --edges --alg PRflux+KSP --alg KSP --alg PRflux -r 0.3 -r 0.6'
        print cmd
        if not opts.printonly:
            subprocess.check_call(cmd.split())

    if opts.weightviz:
        ## plot KSP with and without weighted network.
        topk = 200
        for (pathway,resultdir,datadir) in pathways:
            nodefile = '%s/%s-nodes.txt' % (datadir,pathway)
            edgefile = '%s/%s-edges.txt' % (datadir,pathway)

            unweightedfile = '%s/ksp/%s-q_0.5-none_weighted-ksp_20000_paths.txt' % (resultdir,pathway)
            unweightedgsid = '%s-KSP-unweighted-top%d' % (pathway,topk)
            cmd = 'python ../2014-06-linker/src/post-pathlinker-to-graphspace.py --nodefile %s --edgefile %s --kspfile %s --gsid %s --topk %d' % \
                (nodefile,edgefile,unweightedfile,unweightedgsid,topk)
            print cmd
            if not opts.printonly:
                os.system(cmd)

            weightedfile = '%s/ksp/%s-weighted-q_0.5-none_weighted-ksp_20000_paths.txt' % (resultdir,pathway)
            weightedgsid = '%s-KSP-weighted-top%d' % (pathway,topk)
            cmd = 'python ../2014-06-linker/src/post-weighted-linker-to-graphspace.py --nodefile %s --edgefile %s --kspfile %s --gsid %s --topk %d --weighted' % \
                (nodefile,edgefile,weightedfile,weightedgsid,topk)
            print cmd
            if not opts.printonly:
                os.system(cmd)
          

        ## compute weighted average of positive/negativesets
        cmd = 'python ../2014-06-linker/src/computeAvgEdgeWeights.py'
        print cmd
        if not opts.printonly:
            os.system(cmd)


    ## make venn diagrams if opts.dbcompare is specified.
    if opts.dbcompare and opts.venn:
        cmd = 'python ../2014-06-linker/src/make-venn.py'
        print cmd
        if not opts.printonly:
            subprocess.check_call(cmd.split())

    print 'DONE'
    return

############################################################
def generatePathwaySpecificInteractomes(ppifile):
    edges = [] # list of (u,v,line) tuples
    header = ''
    with open(ppifile) as fin:
        for line in fin:
            if line[0]=='#':
                header = line
                continue
            row = line.split('\t')
            edges.append((row[0],row[1],line))

    checkDir(PPIDIR+'/netpath/')
    pathways = getNetPathPathways(False,False)
    for p in pathways:
        interactomefile = '%s/netpath/%s-interactome.txt' % (PPIDIR,p)
        if not os.path.isfile(interactomefile):
            print 'Making NetPath %s Interactome' % (p)
            nodefile = '%s/%s-nodes.txt' % (NETPATHDIR,p)
            generatePPI(edges,nodefile,interactomefile,header)

    checkDir(PPIDIR+'/wnt-all-receptors/')
    interactomefile = '%s/wnt-all-receptors/Wnt-interactome.txt' % (PPIDIR)
    if not os.path.isfile(interactomefile):
        print 'Making Wnt All receptors Interactome'
        nodefile = 'data/wnt-all-receptors/Wnt-nodes.txt'
        generatePPI(edges,nodefile,interactomefile,header)

    checkDir(PPIDIR+'/kegg/')
    pathways,keggids = getKEGGPathways()
    for p in pathways:
        interactomefile = '%s/kegg/%s-interactome.txt' % (PPIDIR,p)
        if not os.path.isfile(interactomefile):
            print 'Making KEGG %s (%s) Interactome' % (p,keggids[p])
            nodefile = '%s/%s-nodes.txt' % (KEGGDIR,p)
            generatePPI(edges,nodefile,interactomefile,header)
    return

############################################################
def generatePPI(edges,nodefile,interactomefile,header):
    ## get receptors and tfs
    nodes = readColumns(nodefile,1,2)
    receptors = set([n for n,t in nodes if t == 'receptor'])
    tfs = set([n for n,t in nodes if t == 'tf'])    
    out = open(interactomefile,'w')
    out.write(header)
    numskipped = 0
    for u,v,line in edges:
        if u in tfs or v in receptors:
            numskipped+=1
            continue
        out.write(line)
    out.close()    
    
    # print receptors, tfs, num edges removed, percent removed for org file.
    print '| %d | %d | %d | %.2e |' % (len(receptors),len(tfs),numskipped,numskipped/float(len(edges)))
    #print '%d receptors and %d tfs.' % (len(receptors),len(tfs))
    #print 'Removed %d edges (%.2e)' % (numskipped,numskipped/float(len(edges)))
    return

############################################################
def checkDir(dirname):
    if not os.path.isdir(dirname):
        print 'WARNING: %s does not exist. Creating...' % (dirname)
        os.makedirs(dirname)
    return

############################################################
def getNetPathPathways(onlynetpathwnt,dbcompare):
    if onlynetpathwnt:
        return ['Wnt']
    if dbcompare: # only return the 6 pathways in common with KEGG.
        analyzedpathwayfile = 'data/netpath-dbcompare-pathways.txt'
    else:
        analyzedpathwayfile = 'data/netpath-analyzed-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,1)]
    return pathways

############################################################
def getKEGGPathways():
    analyzedpathwayfile = 'data/kegg-analyzed-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,2)]
    # dictionary of keggnames to netpath names.
    kegg2netpath = readDict(analyzedpathwayfile,2,1)
    return pathways,kegg2netpath

############################################################
def runPathLinker(pathway,resultdir,datadir,ppidir,k,forcealg,printonly):
    print '-'*25 + pathway + '-'*25
    
    # node file contains node annotated with 'tf' or 'receptor' or 'none'
    nodefile = '%s/%s-nodes.txt' % (datadir,pathway)

    # create output directory, make sure it exists, and
    # append pathway name for output prefix
    outdir = '%s/pathlinker/' % (resultdir)
    checkDir(outdir)
    outprefix = '%s/%s-' % (outdir,pathway)

    ## pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

    # pathlinker command
    if forcealg or not os.path.isfile('%sk_%d-paths.txt' % (outprefix,k)):
        ## Old script
        ##script = '/home/annaritz/src/python/PathLinker/PathLinker-1.0/hack-scripts-pre-refactoring/PathLinker-NoPR.py'
        ## Needed to pass PPIcosts (-log10) instead of PPIFILE for old script.

        script = '/home/annaritz/src/python/PathLinker/PathLinker-1.0/PathLinker.py'
        cmd = 'python %s -k %d --write-paths --output %s %s %s' % (script,k,outprefix,ppifile,nodefile) 
        print cmd              
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s exists. Use --forcealg to override.' % (pathway,'%sk_%d-paths.txt' % (outprefix,k))
    return

############################################################
def runShortestPaths(pathway,resultdir,datadir,ppidir,forcealg,printonly):
    print '-'*25 + pathway + '-'*25

    # node file contains node annotated with 'tf' or 'receptor' or 'none'
    nodefile = '%s/%s-nodes.txt' % (datadir,pathway)

    # create output directory, make sure it exists, and
    # append pathway name for output filename
    outdir = '%s/shortestpaths/' % (resultdir)
    checkDir(outdir)
    outfile = '%s/%s-shortest-paths.txt' % (outdir,pathway)
               
    ## pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

    if forcealg or not os.path.isfile(outfile):
        script = '/home/annaritz/src/python/CellCycle/shortest_paths.py'
        cmd = 'python %s --network %s --annotations %s --out %s --include-ties' % (script,ppifile,nodefile,outfile)
        print cmd
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s exists. Use --forcealg to override' % (pathway,outfile)
    return

############################################################
def runInducedSubgraph(pathway,resultdir,datadir,ppidir,forcealg,printonly):
    print '-'*25 + pathway + '-'*25

    ## paths file is pathlinker run for k=20,000.
    pathsfile = '%s/pathlinker/%s-k_20000-paths.txt' % (resultdir,pathway)
    if not os.path.isfile(pathsfile):
        sys.exit('ERROR: %s must exist. Run master-script.py with the --pathlinker option.' % (pathsfile))

    # create output directory, make sure it exists, and
    # append pathway name for output filename
    outdir = '%s/inducedsubgraph/' % (resultdir)
    checkDir(outdir)
    outfile = '%s/%s-induced-subgraph.txt' % (outdir,pathway)

    ## pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

    if forcealg or not os.path.isfile(outfile):
        script = '/data/annaritz/signaling/2014-06-linker/src/order-by-induced-subgraph.py'
        cmd = 'python %s --pathsfile %s --outfile %s --ppi %s' % (script,pathsfile,outfile,ppifile)
        print cmd
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s exists. Use --forcealg to override.' % (pathway,outfile)
    return

############################################################
def rerankPathLinker(pathway,resultdir,datadir,ppidir,forcealg,printonly):
    print '-'*25 + pathway + '-'*25

    ## paths file is pathlinker run for k=20,000.
    pathsfile = '%s/pathlinker/%s-k_20000-paths.txt' % (resultdir,pathway)
    if not os.path.isfile(pathsfile):
        sys.exit('ERROR: %s must exist. Run master-script.py with the --pathlinker option.' % (pathsfile))
        
    # create output directory, make sure it exists, and
    # append pathway name for output filename
    outdir = '%s/rerankedpathlinker/' % (resultdir)
    checkDir(outdir)
    outprefix = '%s/%s-reranked-pathlinker' % (outdir,pathway)

    if forcealg or not os.path.isfile('%s-unique-edges_paths.txt' % (outprefix)):
        script = '/data/annaritz/signaling/2014-06-linker/src/recount-ksp.py'
        cmd = 'python %s --pathsfile %s --outputprefix %s' % (script,pathsfile,outprefix)
        print cmd
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s exists. Use --forcealg to override.' % (pathway,'%s-unique-edges_paths.txt' % (outprefix))
    return

############################################################
def runPageRank(pathway,resultdir,datadir,ppidir,q,forcealg,printonly):
    print '-'*25 + pathway + '-'*25

    # node file contains node annotated with 'tf' or 'receptor' or 'none'
    nodefile = '%s/%s-nodes.txt' % (datadir,pathway)

    # create output directory, make sure it exists, and
    # append pathway name for output filename
    outdir = '%s/pagerank/' % (resultdir)
    checkDir(outdir)
    outprefix = '%s/%s-q_%.2f' % (outdir,pathway,q)

    ## pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

    if forcealg or not os.path.isfile('%s-node-pagerank.txt' % (outprefix)):
        ## old script
        ##script = '/home/annaritz/src/python/PathLinker/PathLinker-1.0/hack-scripts-pre-refactoring/PathLinker-PRonly.py'
    
        script = '/home/annaritz/src/python/PathLinker/PathLinker-1.0/PathLinker.py'
        cmd = 'python %s --PageRank -q %s -k %d --output %s %s %s' % (script,q,1,outprefix,ppifile,nodefile)     
        print cmd
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s Exists. Use --forcealg to override.' % (pathway,'%s-pagerank.txt' % (outprefix))
    return

############################################################
def runEQED(pathway,resultdir,datadir,ppidir,inputcurrent,forcealg,printonly):
    print '-'*25 + pathway + '-'*25

    # node file contains node annotated with 'tf' or 'receptor' or 'none'
    nodefile = '%s/%s-nodes.txt' % (datadir,pathway)

    # create output directory, make sure it exists, and
    # append pathway name for output filename
    outdir = '%s/eqed/' % (resultdir)
    checkDir(outdir)
    outprefix = '%s/%s' % (outdir,pathway)

    ## pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

    if forcealg or not os.path.isfile('%s-eqed-edges.out' % (outprefix)):
        script = '/data/annaritz/sig-path-other-methods/src/eQED.py'
        cmd = 'python %s -e %s -n %s -i %d -o %s -l' % (script,ppifile,nodefile,inputcurrent,outprefix)
        print cmd
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s Exists. Use --forcealg to override.' % (pathway,'%s-eqed-edges.out' % (outprefix))
    return


############################################################
def runResponseNet(pathway,resultdir,datadir,ppidir,gamma,forcealg,printonly):
    print '-'*25 + pathway + '-'*25
                   
    # node file contains node annotated with 'tf' or 'receptor' or 'none'
    nodefile = '%s/%s-nodes.txt' % (datadir,pathway)

    # create output directory, make sure it exists, and
    # append pathway name for output filename
    outdir = '%s/reponsenet/' % (resultdir)
    checkDir(outdir)
    outprefix = '%s/%s-gamma_%d' % (outdir,pathway,gamma)

    ## pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

    if forcealg or not os.path.isfile('%s_responsenet-edges.out' % (outprefix)):
        script = '/data/annaritz/sig-path-other-methods/src/ResponseNet.py'
        cmd = 'python %s -e %s -n %s -o %s -g %d' % (script,ppifile,nodefile,outprefix,gamma)
        print cmd
        if not printonly:
            try:
                subprocess.check_call(cmd.split())
            except subprocess.CalledProcessError:
                'Error with gamma=%d: skipping.' % (gamma)
    else:
        print 'Skipping %s w/ gamma %d: %s exists. Use --forcealg to override.' % (pathway,gamma,'%s_responsenet-edges.out' % (outprefix))
    return

############################################################
def runPCSF(pathway,resultdir,datadir,ppidir,prize,omega,forcealg,printonly):
    print '-'*25 + pathway + '-'*25
                   
    # node file contains node annotated with 'tf' or 'receptor' or 'none'
    nodefile = '%s/%s-nodes.txt' % (datadir,pathway)

    # create output directory, make sure it exists, and
    # append pathway name for output filename
    outdir = '%s/pcsf/' % (resultdir)
    checkDir(outdir)
    outprefix = '%s/%s-prize%d-omega%.2f' % (outdir,pathway,prize,omega)

    ## pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)
    
    # run PCSF
    if forcealg or not os.path.isfile('%s_PCSF-edges.out' % (outprefix)):
        script = '/data/annaritz/sig-path-other-methods/src/PCSF_weighted.py'
        cmd = 'python %s -e %s -n %s -o %s -p %d --omega %.2f' % (script,ppifile,nodefile,outprefix,prize,omega)
        print cmd
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s exists. Not running PCSF.' % (pathway,'%s_PCSF-edges.out' % (outprefix))
    return

############################################################
def runANAT(pathway,resultdir,datadir,ppidir,alpha,forcealg,printonly):
    print '-'*25 + pathway + '-'*25
                   
    # node file contains node annotated with 'tf' or 'receptor' or 'none'
    nodefile = '%s/%s-nodes.txt' % (datadir,pathway)

    # create output directory, make sure it exists, and
    # append pathway name for output filename
    outdir = '%s/anat/' % (resultdir)
    checkDir(outdir)
    outprefix = '%s/%s-alpha%.2f' % (outdir,pathway,alpha)

    ## pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

    script = '/data/annaritz/signaling/2014-06-linker/src/run-anat-weighted.py'
    if forcealg or not os.path.isfile('%s-edges.out' % (outprefix)):
        cmd = 'python %s -n %s -a %.2f -o %s --ppi %s'  % (script,nodefile,alpha,outprefix,ppifile)
        #if forcealg: # this is now handled at this function rather than passing to next.
        #    cmd+= ' --force'
        print cmd
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s exists. Not running ANAT.' % (pathway,'%s-edges.out' % (outprefix))
    return

############################################################
## TODO make sure ipa works on an undirected graph.
def runIPA(pathway,resultdir,datadir,ppidir,nmax,forcealg,printonly):
    print '-'*25 + pathway + '-'*25
                   
    # node file contains node annotated with 'tf' or 'receptor' or 'none'
    nodefile = '%s/%s-nodes.txt' % (datadir,pathway)

    # create output directory, make sure it exists, and
    # append pathway name for output filename
    outdir = '%s/ipa/' % (resultdir)
    checkDir(outdir)

    ## pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

    outfile = '%s/%s-nmax%d.out' % (outdir,pathway,nmax)
    if forcealg or not os.path.isfile(outfile):
        script = '/data/annaritz/signaling/2014-06-linker/src/ipa-network-generation.py'
        cmd = 'python %s --ppi %s --nodes %s --nmax %d --outfile %s'  % \
                (script,ppifile,nodefile,nmax,outfile)
        print cmd
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s exists. Use --forcealg to override.' % (pathway,outfile)
    return

############################################################
def postWntReconstructionsToGraphSpace(infile,thres,gsid,printonly,increase=False,decrease=False,undirected=False):
    ppifile = 'data/pathway-specific-interactomes/pathlinker-signaling-children-reg/weighted/wnt-all-receptors/Wnt-interactome.txt'
    cmd = 'python src/post-to-graphspace.py --infile %s --ppi %s --version %s --datadir %s --gsid %s --netpath Wnt --kegg Wnt --addfzd' % \
          (infile,ppifile,PPIVERSION,DATADIR,gsid)
    if increase: # ranked list - pass the threshold
        cmd += ' --increase --thres %f' % (thres)
    if decrease:
        cmd += ' --decrease --thres %f' % (thres)
    if undirected:
        cmd += ' --undirected'
    print cmd
    if not printonly:
        subprocess.check_call(cmd.split())

    ##print unlabeled
    cmd = 'python src/post-to-graphspace.py --infile %s --ppi %s --version %s --datadir %s --gsid %s --netpath Wnt --nolabels --addfzd' % \
          (infile,ppifile,PPIVERSION,DATADIR,gsid+'-nolabels')
    if increase: # ranked list - pass the threshold
        cmd += ' --increase --thres %f' % (thres)
    if decrease:
        cmd += ' --decrease --thres %f' % (thres)
    if undirected:
        cmd += ' --undirected'
    print cmd
    if not printonly:
        subprocess.check_call(cmd.split())

    return

############################################################
def computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,sampleoutprefix,\
                           subsamplefps,forceprecrec,printonly,nodefile=None,nodesortcol=None,descending=False,param=None):
    ## Get true edge file and true node file from the data directory.
    trueedgefile = '%s/%s-edges.txt' % (datadir,pathway)
    truenodefile = '%s/%s-nodes.txt' % (datadir,pathway)
        
    ## make output directory.
    checkDir(outdir)
    if param == None:
        outprefix = '%s%s' % (outdir,pathway)
    else: # e.g., 'q_0.50'
        outprefix = '%s%s-%s' % (outdir,pathway,param)

    ## check to see if the file exists:
    finalfile = '%s-exclude_%s-sample_%dX-edge-precision-recall.txt' % (outprefix,negtype,subsamplefps)
    if not forceprecrec and os.path.isfile(finalfile):
        print 'Skipping %s exclude %s: file exists. Use --forceprecrec to override.' % (pathway,negtype)
        return

    ## pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

    ## compute-precision-recally.py simply uses ppifile to get
    ## all possible edges.
    cmd = 'python src/compute-precision-recall.py --outprefix %s --edgefile %s --trueedgefile %s --truenodefile %s --sampledoutprefix %s --ppi %s --negtype %s --neg-factor %d' % \
          (outprefix,edgefile,trueedgefile,truenodefile,sampleoutprefix,ppifile,negtype,subsamplefps)
    if edgesortcol != None:
        cmd += ' --edgecol %d' % (edgesortcol)
    if nodefile != None:
        cmd += ' --nodefile %s' % (nodefile)
    if nodesortcol != None:
        cmd += ' --nodecol %d' % (nodesortcol)
    if descending:
        cmd += ' --descending'
    #if opts.forceprecrec: # this will cause a resampling of the negatives.  USed for debugging.
    #    cmd += ' --force'
    print cmd
    if not printonly:
        subprocess.check_call(cmd.split())
    return

############################################################
def computeAggregatePrecisionRecall(inputdir,negtype,subsamplefps,forceprecrec,printonly,descending=True,param=None):

    ## check to see if the file exists:
    if param == None:
        finalfile = '%saggregate-exclude_%s-sample_%dX-node-precision-recall.txt' % (inputdir,negtype,subsamplefps)
    else:
       finalfile = '%saggregate-%s-exclude_%s-sample_%dX-node-precision-recall.txt' % (inputdir,param,negtype,subsamplefps)
    if not forceprecrec and os.path.isfile(finalfile):
        print 'Skipping aggregate for NetPath, file exists. Use --forceprecrec to override.'
        return
    cmd = 'python src/compute-aggregate-precision-recall.py --inputdir %s --netpath --negtype %s --neg-factor %d' % (inputdir,negtype,subsamplefps)
    if param != None:
        cmd += ' --param %s' % (param)
    if descending:
        cmd += ' --descending'
    print cmd
    if not printonly:
        subprocess.check_call(cmd.split())
    return

if __name__=='__main__':
    main(sys.argv)
