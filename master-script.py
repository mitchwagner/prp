#!/usr/bin/python

## Anna Ritz, March 2015
## Master Script for running PathLinker for Nat. Methods submisssion.

from optparse import OptionParser,OptionGroup
from utilsPoirel import *
import sys
import os
import os.path
import subprocess

##############################################################
## GLOBAL VARIABLES

## PPIVERION is going to be populated in main() with one of the
## ALLOWEDVERSIONS.  We have settled on using "pathlinker-signalingg-children-reg",
## which is the 2015pathlinker interactome weighted using the 8 children of the 
## "signal transduction" GO term as well as "regulation of signal transduction."
## All other versions are commented out.
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

## DATADIR is the path to the data/ directory checked into SVN.  
## The interactomes, KEGG, and NetPath edge files are all checked in.
DATADIR = '/data/annaritz/datasets/svn-data/'

## PPIDIR is the directory that contains the pathway-specific interactomes.
## For each NetPath pathway, we take the PPIVERSION network and remove incoming 
## edges to receptors and outgoing edges from tfs.  Thus, there is one interactome
## for each NetPath pathway in the PPIDIR directory.  Note that 
## "pathway-specific" is a bit too strong of a phrase, and shouldn't be used
## in the paper. 
PPIDIR = '/data/annaritz/projects/2015-03-pathlinker/data/pathway-specific-interactomes/'

## MAPPING VARIABLES
## There is a human-gene-map.txt (originally from Chris) that is checked into
## SVN.  This is used to map between UniProtKB and the common name.
## TODO: Once CSBDB is bug-free, we should use CSBDB to do this mapping.  However
## that will require changing many of the auxilliary scripts.
MAPPINGFILE = '%s/namespace-mappers/human-gene-map.txt' % (DATADIR)
MAPFROMCOL = 6
MAPTOCOL = 1

## PATHWAY DIRECTORIES
## NetPath and KEGG edge files (*-edges.txt) are checked into SVN.
NETPATHDIR='%s/interactions/netpath/pathways/' % (DATADIR) 
KEGGDIR='%s/interactions/kegg/2015-03-23/hsa/edge-files/' % (DATADIR)

## VARYPARAMS
## To select parameters for competing methods, we vary each parameter.
## The VARYPARAMS dictionary provides the parameter list.
VARYPARAMS = {'q': [0.1, 0.25, 0.5, 0.75, 0.9],  ## PageRank teleportation probability
              'gamma':[10, 15, 20, 25, 30],      ## ResponseNet parameter
              'omega':[0, 0.01, 0.1],            ## PCSF penalty for adding trees
              'prize':[1, 3, 5, 7, 9],           ## PCSF prize
              'alpha':[0, 0.1, 0.25, 0.4, 0.5],  ## ANAT parameter
              'nmax':[5,10,15, 25, 35, 50, 75, 100, 200, 500], ## IPA parameter
          }

##############################################################
## The main method parses all parameters and runs all experiments.
def main(args):
    global PPIVERSION, PPIDIR

    ## parse arguments
    opts = parseArguments(args)

    ## set the PPIVERSION. 
    if opts.ppiversion not in ALLOWEDVERSIONS:
        sys.exit('ERROR: --ppiversion must be one of %s. You have %s. Exiting.' % (','.join(ALLOWEDVERSIONS),opts.ppiversion))
    PPIVERSION = opts.ppiversion

    ## BACKGROUND INTERACTOME
    ## set PPI, determined by whether the --weightedppi argument is specified.
    ## If --weightedppi, then the ppifile contains scores/probabilities 
    ## (*-weighted.txt).  If not --weightedppi, then ppifile contains
    ## scores of 1 for all the edges.
    ## We also set the result directory prefix, resultprefix, and the
    ## directory of pathway-specific interactomes (PPIDIR)
    if opts.weightedppi:
        ppifile = '/data/annaritz/datasets/svn-data/interactomes/human/%s-weighted.txt' % (PPIVERSION)
        resultprefix = 'results/%s/weighted/' % (PPIVERSION)
        PPIDIR = '%s/%s/weighted/' % (PPIDIR,PPIVERSION)
    else:
        ppifile = '/data/annaritz/datasets/svn-data/interactomes/human/%s.txt' % (PPIVERSION)
        resultprefix = 'results/%s/unweighted/' % (PPIVERSION)
        PPIDIR = '%s/%s/unweighted/' % (PPIDIR,PPIVERSION)
    ## Make sure the directories exist; if they don't, create them.
    checkDir(resultprefix)
    checkDir(PPIDIR)

    ## PATHWAY-SPECIFIC INTERACTOMES
    generatePathwaySpecificInteractomes(ppifile)

    # DATASETS
    ## pathways is a set() of (pathwayname,resultdir,datadir,ppidir) tuples.
    pathways = set()
    kegg2netpath = {} # will be populated if kegg pathways are specified.
    if opts.netpath or opts.onlynetpathwnt or opts.dbcompare:
        pathwaynames = getNetPathPathways(opts.onlynetpathwnt,opts.dbcompare)
        resultdir = '%s/netpath/' %(resultprefix)
        checkDir(resultdir) 
        datadir = NETPATHDIR
        ppidir = '%s/netpath/' % (PPIDIR)
        pathways.update([(p,resultdir,datadir,ppidir) for p in pathwaynames])
        print 'Running %d NetPath pathways' % (len(pathwaynames))
    if opts.wntforexperiments:
        resultdir = '%s/wnt-all-receptors/' % (resultprefix)
        checkDir(resultdir)
        datadir = 'data/wnt-all-receptors/'
        ppidir = '%s/wnt-all-receptors/' % (PPIDIR)
        pathways.update([('Wnt',resultdir,datadir,ppidir)])
        print 'Running 1 wnt-all-receptors pathway'
    if opts.kegg or opts.dbcompare:
        pathwaynames,kegg2netpath = getKEGGPathways() 
        resultdir = '%s/kegg/' % (resultprefix)
        checkDir(resultdir) 
        datadir = KEGGDIR
        ppidir = '%s/kegg/'% (PPIDIR)
        pathways.update([(p,resultdir,datadir,ppidir) for p in pathwaynames])
        print 'Running %d KEGG pathways' % (len(pathwaynames))
    if len(pathways)==0:
        print 'WARNING: no datasets specified. This is OK when visualizing output or running special runs.'
    else:
        print '%d total pathways considered.' % (len(pathways))

    ## ALGORITHMS ##
    ## For each algorithm, iterate through pathways and call the run() method.
    ## If --varyparams is specified, iterate through (pathway,param) combinations
    ## and cal the run() method.

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
        ## Now, always vary the nmax parameter; all param values are always plotted.
        #if not opts.varyparams: # just run with opts.nmax value
        #    for (pathway,resultdir,datadir,ppidir) in pathways:
        #        runIPA(pathway,resultdir,datadir,ppidir,opts.nmax,opts.forcealg,opts.printonly)
        #else: # vary nmax
        for varynmax in VARYPARAMS['nmax']:
            for (pathway,resultdir,datadir,ppidir) in pathways:
                runIPA(pathway,resultdir,datadir,ppidir,varynmax,opts.forcealg,opts.printonly)
        print 'Done IPA\n'

    ## VIZ SCRIPTS ##
    ## These are a "grab bag" of precision-recall computations, precision-recall
    ## plots, and other types of analyses.

    ## write precision/recall
    if opts.computeprecrec:
        print 'Write Precision Recall'

        ## Always run 'exclude-none' and 'exclude-adjacent'.
        ## If --ignorekeggpositives, also pass KEGG nodes and
        ## edges as ignored files along with the 'exclude-file'
        ## parameter.  Add 'file' to the negtypes.
        negtypes = ['none','adjacent']
        if opts.ignorekeggpositives:
            negtypes.append('file')

        ## Compute precision and recall for each negtype described above.
        for negtype in negtypes:

            ## All algorithms use the same set of positives and subsampled negatives 
            ## for the dataset specified.  The sampledir variable points to the directory,
            ## depending on the PPIVERSION and the negtype.  If --netpathkeggunion is 
            ## specified, this is a new set of positives (union of NetPath and KEGG).  Thus,
            ## the sample directory is prepended with "netpathkeggunion".
            if opts.netpathkeggunion:
                sampledir = '%s/netpathkeggunion-samples-exclude-%s' % (resultprefix,negtype)
            else:
                sampledir = '%s/samples-exclude-%s' % (resultprefix,negtype)
            checkDir(sampledir)

            ## If --wntforexperimets is specified, then there is one dataset in pathways variable
            ## that corresponds to wnt-all-receptors.  Since there is a different set of positives 
            ## than in NetPath Wnt pathway (e.g., FZD4/FZD6 are added), we need to subsample a 
            ## different set of negatives.  Thus, the wntsampledir directory is prepended with 
            ## "wntforexperiments"
            if opts.wntforexperiments:
                wntsampledir = '%s/wntforexperiments-samples-exclude-%s' % (resultprefix,negtype)
                checkDir(wntsampledir)
            
            ## For every algorithm, do the following:
            ## (1) determine the columns to sort on
            ## (2) for each pathway:
            ##   - get the file that contains the predicted edges (and possibly node file as well)
            ##   - determine the output directory: it is either:
            ##         <resultdir>/precision-recall/netpathkeggunion/<alg>/ if --netpathkeggunion is specified
            ##         <resultdir>/precision-recall/<alg>/ otherwise.
            ##   - determine the prefix of the subsampled negative files. It is either:
            ##         <wntsampledir>/<pathway> if --wntforexperiments is specified
            ##         <sampledir>/<pathway> otherwise
            ##   - compute precision and recall
            ## (3) if --netpath is specified, also compute aggregate precision and recall
            ##
            ## If --varyparams is specified, executes (2) and (3) for every parameter value.

            ## PATHLINKER ##
            if opts.pathlinker:
                sortcol = 3 # ksp
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/pathlinker/%s-k_%s-ranked-edges.txt' % (resultdir,pathway,opts.k)
                    outdir = getPRoutdir('pathlinker',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                           union=opts.netpathkeggunion)
                if opts.netpath:
                    ## compute aggregate for NetPath
                    inputdir = '%s/netpath/precision-recall/pathlinker/' % (resultprefix)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,opts.forceprecrec,opts.printonly)
                    
            ## Shortest Paths ##
            if opts.shortestpaths:
                sortcol = None # no sorting; entire file is shortest paths subgraph.
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/shortestpaths/%s-shortest-paths.txt' % (resultdir,pathway)
                    outdir = getPRoutdir('shortestpaths',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                           union=opts.netpathkeggunion)
                if opts.netpath:
                    ## compute aggregate for NetPath
                    inputdir = '%s/netpath/precision-recall/shortestpaths/' % (resultprefix)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,opts.forceprecrec,opts.printonly)
            
            ## Induced Subgraph ##
            if opts.inducedsubgraph:
                sortcol = 3 # ksp
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/inducedsubgraph/%s-induced-subgraph.txt' % (resultdir,pathway)
                    outdir = getPRoutdir('inducedsubgraph',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                           union=opts.netpathkeggunion)
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
                        outdir = getPRoutdir('pagerank',resultdir,opts.netpathkeggunion)
                        sampleoutprefix = getPRsubsampleprefix(wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,nodefile=nodefile,nodesortcol=nodesortcol,\
                                               descending=True,param=param,\
                                               union=opts.netpathkeggunion)
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
                    outdir = getPRoutdir('eqed',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                           opts.printonly,nodefile=nodefile,nodesortcol=nodesortcol,descending=True,\
                                           union=opts.netpathkeggunion)
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
                        outdir = getPRoutdir('responsent',resultdir,opts.netpathkeggunion)
                        sampleoutprefix = getPRsubsampleprefix(wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,param=param,\
                                               union=opts.netpathkeggunion)
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
                        outdir = getPRoutdir('pcsf',resultdir,opts.netpathkeggunion)
                        sampleoutprefix = getPRsubsampleprefix(wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,param=param,\
                                               union=opts.netpathkeggunion)
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
                        outdir = getPRoutdir('anat',resultdir,opts.netpathkeggunion)
                        sampleoutprefix = getPRsubsampleprefix(wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,param=param,\
                                               union=opts.netpathkeggunion)
                    if opts.netpath:
                        ## compute aggregate for NetPath
                        inputdir = '%s/netpath/precision-recall/anat/' % (resultprefix)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param)
            ## DEGREE ##
            if opts.degree:
                sys.exit('ERROR: DEGREE is not implemented in this version. Exiting.')

            ## IPA ##
            if opts.ipa:
                sortcol = None # no sort column; take entire file
                # always run all params, since these are always plotted.
                #if not opts.varyparams: # just run with opts.nmax value
                #    params = ['nmax%d' % (opts.nmax)]
                #else: # vary opts.nmax value
                params = ['nmax%d' % p for p in VARYPARAMS['nmax']]
                for param in params:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        edgefile = '%s/ipa/%s-%s.out' % (resultdir,pathway,param)
                        outdir = getPRoutdir('ipa',resultdir,opts.netpathkeggunion)
                        sampleoutprefix = getPRsubsampleprefix(wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,descending=True,param=param,\
                                               union=opts.netpathkeggunion)
                    if opts.netpath:
                        ## compute aggregate for NetPath
                        inputdir = '%s/netpath/precision-recall/ipa/' % (resultprefix)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        descending=True,param=param)

    ## Plot precision and recall, once the values have been computed.
    if opts.precrecviz:
        print 'Plot Precision Recall'
        ## The algcmd variable contains all of the algorithms that have been specified in the command
        ## line. These are the algorithms that will be plotted in the precision-recall curves.
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

        ## Speficy the input directory to read the precision-recall values from.
        ## if --netpathkeggunion is specified, add the netpathkeggunion directory
        ## so we look in the correct place.
        if opts.netpathkeggunion:
            indir = '%s/netpath/precision-recall/netpathkeggunion/' % (resultprefix)
        else:
            indir = '%s/netpath/precision-recall/' % (resultprefix)

        ## For each pathway, determine the output prefix, construct the call to 
        ## plot-precision-recall.py, and execute it.
        for (pathway,resultdir,datadir,ppidir) in pathways:
            ## if --wntforexperiments is specified, add "all-receptors" 
            ## to label. If --ignorekeggpositives is specified, add "ignorekeggpositives"
            ## to label.  if --netpathkeggunion is specified, add "netpathkeggunion" to label.
            ## Otherwise, label is simply the pathway name.
            if opts.wntforexperiments and 'wnt-all-receptors' in resultdir:
                outprefix = 'viz/precision-recall/%s-all-receptors' % (pathway)
            elif opts.ignorekeggpositives:
                outprefix = 'viz/precision-recall/%s-ignorekeggpositives' % (pathway)
            elif opts.netpathkeggunion:
                outprefix = 'viz/precision-recall/%s-netpathkeggunion' % (pathway)
            else:
                outprefix = 'viz/precision-recall/%s' % (pathway)

            ## Consruct the command.  
            cmd = 'python src/plot-precision-recall.py --indir %s --outprefix %s --pathway %s %s' % \
                  (indir,outprefix,pathway,algcmd)
            if opts.ignorekeggpositives:
                cmd += ' --ignorefromfile'

            ## Only plot the varying parameters for Wnt.  This reduces the amount of 
            ## clutter in the viz/ directory.
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
                
        ## If --netpath is specified, plot the aggregate precision-recall plots.
        if opts.netpath: 
            outprefix = 'viz/precision-recall/aggregate'
            cmd = 'python src/plot-precision-recall.py --indir %s --outprefix %s --pathway aggregate %s' % \
                  (indir,outprefix,algcmd)
            if opts.varyparams:
                cmd += '  --varyparams'
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())

    ## Post Wnt Pathways to Graphspace
    ## Post two different Wnt pathway runs.
    ## - Post NetPath Wnt pathways, used to compute precision and recall.  These are prepended with "pr"
    ## - Post wnt-all-receptors pathways, used to explore false positives for experimental followup
    if opts.graphspace:
        print 'Posting to GraphSpace...'
        
        if opts.pathlinker:
            infile = '%s/wnt-all-receptors/pathlinker/Wnt-k_%d-ranked-edges.txt' % (resultprefix,opts.k)
            gsid = 'Wnt-pathlinker-top%dpaths' % (opts.topk)
            postWntReconstructionsToGraphSpace(infile,opts.topk,gsid,opts.printonly,increase=True)

            infile = '%s/netpath/pathlinker/Wnt-k_%d-ranked-edges.txt' % (resultprefix,opts.k)
            gsid = 'pr-Wnt-pathlinker-top%dpaths' % (opts.topk)
            postWntReconstructionsToGraphSpace(infile,opts.topk,gsid,opts.printonly,increase=True,allreceptors=False)

        if opts.pagerank: # Manually-determined threshold
            ## threshold is set to 154, 270 edges (200,800 paths)
            for thres in [0.001677,0.0003247]:
                infile = '%s/wnt-all-receptors/pagerank/Wnt-q_0.50-edge-fluxes.txt' % (resultprefix)
                gsid = 'Wnt-pagerank-thres%f' % (thres)
                postWntReconstructionsToGraphSpace(infile,thres,gsid,opts.printonly,decrease=True)
                
            ## threshold is set to 153, 331 edges (200,800 paths)
            for thres in [0.001335,0.0002441]:
                infile = '%s/netpath/pagerank/Wnt-q_0.50-edge-fluxes.txt' % (resultprefix)
                gsid = 'pr-Wnt-pagerank-thres%f' % (thres)
                postWntReconstructionsToGraphSpace(infile,thres,gsid,opts.printonly,decrease=True,allreceptors=False)

        if opts.anat:
            infile = '%s/wnt-all-receptors/anat/Wnt-alpha0.00-edges.out' % (resultprefix)
            gsid = 'Wnt-anat-alpha0.00'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly)

            infile = '%s/netpath/anat/Wnt-alpha0.00-edges.out' % (resultprefix)
            gsid = 'pr-Wnt-anat-alpha0.00'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly,allreceptors=False)

        if opts.shortestpaths:
            infile = '%s/wnt-all-receptors/shortestpaths/Wnt-shortest-paths.txt' % (resultprefix)
            gsid = 'Wnt-shortest-paths'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly)

            infile = '%s/netpath/shortestpaths/Wnt-shortest-paths.txt' % (resultprefix)
            gsid = 'pr-Wnt-shortest-paths'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly,allreceptors=False)

        if opts.responsenet:
            infile = '%s/wnt-all-receptors/reponsenet/Wnt-gamma_20_responsenet-edges.out' % (resultprefix)
            gsid = 'Wnt-responsenet-gamma20'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly)

            infile = '%s/netpath/reponsenet/Wnt-gamma_20_responsenet-edges.out' % (resultprefix)
            gsid = 'pr-Wnt-responsenet-gamma20'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly,allreceptors=False)

        if opts.pcsf:
            infile = '%s/wnt-all-receptors/pcsf/Wnt-prize5-omega0.01_PCSF-edges.out' % (resultprefix)
            gsid = 'Wnt-pcsf-prize5-omega0.01'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly)

            infile = '%s/netpath/pcsf/Wnt-prize5-omega0.01_PCSF-edges.out' % (resultprefix)
            gsid = 'pr-Wnt-pcsf-prize5-omega0.01'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly,allreceptors=False)

        if opts.ipa:
            infile = '%s/wnt-all-receptors/ipa/Wnt-nmax10.out' % (resultprefix)
            gsid = 'Wnt-ipa-nmax10'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly,undirected=True)

            infile = '%s/netpath/ipa/Wnt-nmax10.out' % (resultprefix)
            gsid = 'pr-Wnt-ipa-nmax10'
            postWntReconstructionsToGraphSpace(infile,None,gsid,opts.printonly,undirected=True,allreceptors=False)

    ## RANK TFS
    ## Little script that ranks the TRs in PathLinker predictions vs. PageRank predictions
    ## TRs that have a path from a receptor in IPA are starred.
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
        if not os.path.isfile('data/shortest-paths-for-false-positives/Wnt-dist.txt'):
            cmd = 'python src/shortest-paths-for-false-positives.py'
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())
        else:
            print 'Shortest Paths from nodes to pathways are already computed.'

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
        
        if not os.path.isfile('data/shortest-paths-for-false-negatives/Wnt-dist.txt'):
            cmd = 'python src/shortest-paths-for-false-negatives.py'
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())
        else:
            print 'Shortest Paths from/to nodes are already computed.'

        cmd = 'python ../2014-06-linker/src/plot-false-negatives.py -o viz/false-negatives/false-negatives-aggregate- --prefix viz/precrecfiles-sample-once-per-pathway/precrec-exclude_none --pathway aggregate --pdf --edges --alg PRflux+KSP --alg KSP --alg PRflux -r 0.3 -r 0.6'
        print cmd
        if not opts.printonly:
            subprocess.check_call(cmd.split())

        cmd = 'python ../2014-06-linker/src/plot-false-negatives.py -o viz/false-negatives/false-negatives-Wnt- --prefix viz/precrecfiles-sample-once-per-pathway/precrec-exclude_none --pathway Wnt --pdf --edges --alg PRflux+KSP --alg KSP --alg PRflux -r 0.3 -r 0.6'
        print cmd
        if not opts.printonly:
            subprocess.check_call(cmd.split())

    # if opts.weightviz:
    #     ## plot KSP with and without weighted network.
    #     topk = 200
    #     for (pathway,resultdir,datadir) in pathways:
    #         nodefile = '%s/%s-nodes.txt' % (datadir,pathway)
    #         edgefile = '%s/%s-edges.txt' % (datadir,pathway)

    #         unweightedfile = '%s/ksp/%s-q_0.5-none_weighted-ksp_20000_paths.txt' % (resultdir,pathway)
    #         unweightedgsid = '%s-KSP-unweighted-top%d' % (pathway,topk)
    #         cmd = 'python ../2014-06-linker/src/post-pathlinker-to-graphspace.py --nodefile %s --edgefile %s --kspfile %s --gsid %s --topk %d' % \
    #             (nodefile,edgefile,unweightedfile,unweightedgsid,topk)
    #         print cmd
    #         if not opts.printonly:
    #             os.system(cmd)

    #         weightedfile = '%s/ksp/%s-weighted-q_0.5-none_weighted-ksp_20000_paths.txt' % (resultdir,pathway)
    #         weightedgsid = '%s-KSP-weighted-top%d' % (pathway,topk)
    #         cmd = 'python ../2014-06-linker/src/post-weighted-linker-to-graphspace.py --nodefile %s --edgefile %s --kspfile %s --gsid %s --topk %d --weighted' % \
    #             (nodefile,edgefile,weightedfile,weightedgsid,topk)
    #         print cmd
    #         if not opts.printonly:
    #             os.system(cmd)
          

    #     ## compute weighted average of positive/negativesets
    #     cmd = 'python ../2014-06-linker/src/computeAvgEdgeWeights.py'
    #     print cmd
    #     if not opts.printonly:
    #         os.system(cmd)


    ## make venn diagrams
    if opts.venn:
        cmd = 'python ../2014-06-linker/src/make-venn.py'
        print cmd
        if not opts.printonly:
            subprocess.check_call(cmd.split())

    print 'DONE'
    return

############################################################
## Parses (lots and lots) of options.  
def parseArguments(args):
    usage = 'master-script.py [options]\n'
    parser = OptionParser(usage=usage)

    ## Common options
    parser.add_option('','--forcealg',action='store_true',default=False,\
                          help='Run algorithms even if they will overwrite existing files.  Default is that algorithms are not run if output files exist.')
    parser.add_option('','--forceprecrec',action='store_true',default=False,\
                      help='Run computations of precision and recall, even if they will overwrite existing files. Default is that precision and recall are not run if output files exist.')
    parser.add_option('','--printonly',action='store_true',default=False,\
                          help='Print the commands to stdout, but do not execute.')

    ## Datasets.  
    group = OptionGroup(parser,'Datasets')
    group.add_option('','--ppiversion',type='string',default='pathlinker-signaling-children-reg',\
                     help='Version of the PPI to run.  Options are %s. Default is "pathlinker-signaling-children-reg."' % \
                     (', '.join(ALLOWEDVERSIONS)))
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

    ## Algorithms
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

    ## Precision Recall and Other Visualizations
    group = OptionGroup(parser,'Visualizations')
    group.add_option('','--computeprecrec',action='store_true',default=False,\
                         help='Compute precision and recall curves; write to file.')
    group.add_option('','--precrecviz',action='store_true',default=False,\
                         help='Display precision recall curves.')
    group.add_option('','--ignorekeggpositives',action='store_true',default=False,\
                     help='Ignore KEGG positives when computing & visualizing precision/recall. Will ignore pathways that do not appear in KEGG.')
    group.add_option('','--netpathkeggunion',action='store_true',default=False,\
                     help='Compute & visualize precision/recall using both NetPath and KEGG as positives.')
    group.add_option('','--ranktfs',action='store_true',default=False,\
                     help='Plot TF ranks for PageRank, PathLinker, and IPA.')
    group.add_option('','--bullseye',action='store_true',default=False,\
                         help='Make bullseye plots.')
    group.add_option('','--falsenegs',action='store_true',default=False,\
                         help='Make false negative plots.')
    group.add_option('','--venn',action='store_true',default=False,\
                         help='plot venn diagrams. Only with dbpcompare.')
    group.add_option('','--graphspace',action='store_true',default=False,\
                         help='post Wnt networks to graphspace.')
    group.add_option('','--paper',action='store_true',default=False,\
                         help='Make plots for paper.')
    group.add_option('','--weightviz',action='store_true',default=False,\
                         help='Post weighted vs. unweighted graph to GraphSpace')
    parser.add_option_group(group)

    ## Optional Arguments for Algorithms
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

    return opts

############################################################
## For each pathway in KEGG, NetPath, and the wnt-all-receptors
## datasets, remove the incoming edges to receptors and outgoing
## edges to TRs and save them in the PPIDIR.
## ppifile: original PPI file (e.g., pathlinker-signaling-children-reg-weighted.txt)
def generatePathwaySpecificInteractomes(ppifile):
    ## Read original PPI file 
    edges = [] # list of (u,v,line) tuples
    header = ''
    with open(ppifile) as fin:
        for line in fin:
            if line[0]=='#':
                header = line
                continue
            row = line.split('\t')
            edges.append((row[0],row[1],line))

    ## Make NetPath interactomes, if not already present
    checkDir(PPIDIR+'/netpath/')
    pathways = getNetPathPathways(False,False)
    for p in pathways:
        interactomefile = '%s/netpath/%s-interactome.txt' % (PPIDIR,p)
        if not os.path.isfile(interactomefile):
            print 'Making NetPath %s Interactome' % (p)
            nodefile = '%s/%s-nodes.txt' % (NETPATHDIR,p)
            generatePPI(edges,nodefile,interactomefile,header)

    ## Make wnt-all-receptors interactome, if not already present.
    checkDir(PPIDIR+'/wnt-all-receptors/')
    interactomefile = '%s/wnt-all-receptors/Wnt-interactome.txt' % (PPIDIR)
    if not os.path.isfile(interactomefile):
        print 'Making Wnt All receptors Interactome'
        nodefile = 'data/wnt-all-receptors/Wnt-nodes.txt'
        generatePPI(edges,nodefile,interactomefile,header)

    ## Make KEGG interactomes, if not already present.
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
## generatePPI makes the pathway-specific interactome.
## edges: original PPI edges
## nodefile: file of nodes (to get receptors and TRs)
## interactomefile: name of output file
## header: common header to put at the top of the output file
def generatePPI(edges,nodefile,interactomefile,header):
    ## get receptors and tfs
    nodes = readColumns(nodefile,1,2)
    receptors = set([n for n,t in nodes if t == 'receptor'])
    tfs = set([n for n,t in nodes if t == 'tf'])    
    
    ## Write edges to file if they are not incoming
    ## to receptors or outgoing from TFs.
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
## checkDir checks to see if a directory exists. If it doesn't,
## then it is created and a warning is output.
def checkDir(dirname):
    if not os.path.isdir(dirname):
        print 'WARNING: %s does not exist. Creating...' % (dirname)
        os.makedirs(dirname)
    return

############################################################
## getNetPathPathways reads the pathways for NetPath and returns
## them as a list.
## onlynetpathwnt: if True, only return Wnt
## dbcompare: if True, only return the 6 pathways in common with KEGG
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
## getKEGGPathways reads the pathways for KEGG and returns
## them as a list. It also returns a dictionary of {keggid:netpathname}
def getKEGGPathways():
    analyzedpathwayfile = 'data/kegg-analyzed-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,2)]
    # dictionary of keggnames to netpath names.
    kegg2netpath = readDict(analyzedpathwayfile,2,1)
    return pathways,kegg2netpath

############################################################
## Run PathLinker
## pathway: pathway to run (e.g., Wnt)
## resultdir: directory for results.
## datadir: directory to find positives for the pathway
## ppidir: pathway-specific PPI (e.g., Wnt-interactome.txt)
## k: # of paths to run
## forcealg: if True, will not skip over pre-written files.
## printonly: if True, will never execute command.
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
        script = '/home/annaritz/src/python/PathLinker/PathLinker-1.0/PathLinker.py'
        cmd = 'python %s -k %d --write-paths --output %s %s %s' % (script,k,outprefix,ppifile,nodefile) 
        print cmd              
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s exists. Use --forcealg to override.' % (pathway,'%sk_%d-paths.txt' % (outprefix,k))
    return

############################################################
## Run ShortestPaths
## pathway: pathway to run (e.g., Wnt)
## resultdir: directory for results.
## datadir: directory to find positives for the pathway
## ppidir: pathway-specific PPI (e.g., Wnt-interactome.txt)
## forcealg: if True, will not skip over pre-written files.
## printonly: if True, will never execute command.
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
        cmd = 'python %s --network %s --annotations %s --out %s --include-ties --weight --log-transform' % (script,ppifile,nodefile,outfile)
        print cmd
        if not printonly:
            subprocess.check_call(cmd.split())
    else:
        print 'Skipping %s: %s exists. Use --forcealg to override' % (pathway,outfile)
    return

############################################################
## Run Induced Subgraph
## pathway: pathway to run (e.g., Wnt)
## resultdir: directory for results.
## datadir: directory to find positives for the pathway
## ppidir: pathway-specific PPI (e.g., Wnt-interactome.txt)
## forcealg: if True, will not skip over pre-written files.
## printonly: if True, will never execute command.
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
## Rerank PathLinker predictions by paths with AT LEAST ONE new edge
## pathway: pathway to run (e.g., Wnt)
## resultdir: directory for results.
## datadir: directory to find positives for the pathway
## ppidir: pathway-specific PPI (e.g., Wnt-interactome.txt)
## forcealg: if True, will not skip over pre-written files.
## printonly: if True, will never execute command.
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
## Run PageRank
## pathway: pathway to run (e.g., Wnt)
## resultdir: directory for results.
## datadir: directory to find positives for the pathway
## ppidir: pathway-specific PPI (e.g., Wnt-interactome.txt)
## q: teleportation probability
## forcealg: if True, will not skip over pre-written files.
## printonly: if True, will never execute command.
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
## Runs EQED
## pathway: pathway to run (e.g., Wnt)
## resultdir: directory for results.
## datadir: directory to find positives for the pathway
## ppidir: pathway-specific PPI (e.g., Wnt-interactome.txt)
## inputcurrent: amount of current to "inject" into receptors
## forcealg: if True, will not skip over pre-written files.
## printonly: if True, will never execute command.
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
## Run ResponseNet
## pathway: pathway to run (e.g., Wnt)
## resultdir: directory for results.
## datadir: directory to find positives for the pathway
## ppidir: pathway-specific PPI (e.g., Wnt-interactome.txt)
## gamma: parameter that varies the penalty for additional edges with flow
## forcealg: if True, will not skip over pre-written files.
## printonly: if True, will never execute command.
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
## Run PCSF
## pathway: pathway to run (e.g., Wnt)
## resultdir: directory for results.
## datadir: directory to find positives for the pathway
## ppidir: pathway-specific PPI (e.g., Wnt-interactome.txt)
## prize: prize to place on TFs (terminal nodes)
## omega: penalty for adding aditional trees to the forest
## forcealg: if True, will not skip over pre-written files.
## printonly: if True, will never execute command.
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
## Run ANAT
## pathway: pathway to run (e.g., Wnt)
## resultdir: directory for results.
## datadir: directory to find positives for the pathway
## ppidir: pathway-specific PPI (e.g., Wnt-interactome.txt)
## alpha: parameter for the tradeoff between shortest-paths and steiner trees.
## forcealg: if True, will not skip over pre-written files.
## printonly: if True, will never execute command.
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
## Run IPA
## pathway: pathway to run (e.g., Wnt)
## resultdir: directory for results.
## datadir: directory to find positives for the pathway
## ppidir: pathway-specific PPI (e.g., Wnt-interactome.txt)
## nmax: maximum sub-network size
## forcealg: if True, will not skip over pre-written files.
## printonly: if True, will never execute command.
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
## Auxiliary function to get the Precision/Recall output 
## directory, depending on whether --netpathkeggunion is
## specified.
def getPRoutdir(alg,resultdir,netpathkeggunion):
    if netpathkeggunion:
        return '%s/precision-recall/netpathkeggunion/%s/' % (resultdir,alg)
    else:
        return '%s/precision-recall/%s/' % (resultdir,alg)

############################################################
## Auxiliary function to get the PRecision/Recall subsample
## file prefix, depending on whether this result directory
## is for --wntforexperiments (contains "wnt-all-receptors")
def getPRsubsampleprefix(wntsampledir,sampledir,pathway):
    if 'wnt-all-receptors' in resultdir:
        return '%s/%s' % (wntsampledir,pathway)
    else:
        return '%s/%s' % (sampledir,pathway)

############################################################
## Posts Wnt Reconstructions to GraphSpace
## Posts two different visualizations:
## - Post annotated networks with gene names and ranked value, if predictions are ranked.
## - Post un-annotated networks with labels removed.
##
## infile: file of edges, which may be ranked by a 3rd column
## thres: if increase==True or decrease==True, use this threshold
## gsid: GraphSpace ID
## printonly: if True, will never execute command.
## increase: if True, rank edges in inceasing order
## decrease: if True, rank edges in decreasing order
## (note: if increase==False and decrease==False, then the edges
## are considered an entire set and not ranked by thres)
## undirected: if True, checks both (u,v) and (v,u) for evidence sources
## allreceptors: if True, takes Wnt interactome from wnt-all-receptors/ instead of netpath/
def postWntReconstructionsToGraphSpace(infile,thres,gsid,printonly,increase=False,\
            decrease=False,undirected=False,allreceptors=True):
    if allreceptors:
        ppifile = 'data/pathway-specific-interactomes/pathlinker-signaling-children-reg/weighted/wnt-all-receptors/Wnt-interactome.txt'
    else:
        ppifile = 'data/pathway-specific-interactomes/pathlinker-signaling-children-reg/weighted/netpath/Wnt-interactome.txt'

    ## print annotated
    cmd = 'python src/post-to-graphspace.py --infile %s --ppi %s --version %s --datadir %s --gsid %s --netpath Wnt --kegg Wnt ' % \
          (infile,ppifile,PPIVERSION,DATADIR,gsid)
    if allreceptors: ## add FZD4/FZD6 as receptors
        cmd += ' --addfzd'
    if increase: # ranked list - pass the threshold
        cmd += ' --increase --thres %f' % (thres)
    if decrease: # ranked list - pass the threshold
        cmd += ' --decrease --thres %f' % (thres)
    if undirected: 
        cmd += ' --undirected'
    print cmd
    if not printonly:
        subprocess.check_call(cmd.split())

    ##print unlabeled
    cmd = 'python src/post-to-graphspace.py --infile %s --ppi %s --version %s --datadir %s --gsid %s --netpath Wnt --kegg Wnt --nolabels' % \
          (infile,ppifile,PPIVERSION,DATADIR,gsid+'-nolabels')
    if allreceptors:
        cmd += ' --addfzd'
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
## Computes precision and recall and writes values to file.
## pathway: pathway to compute precision and recall for (e.g., Wnt).
## datadir: directory for true edge file and true node file for pathway
## ppidir: directory for pathway-specific interactomes
## edgefile: predicted edges
## outdir: output directory
## edgesortcol: sort column for edges. If None, then take edges as a set.
## negtype: one of 'none','adjacent', or 'file'.
## sampleoutprefix: prefix of subsampled negatives and positives for the pathway
## subsamplefps: Number of negatives to sample (a factor of the size of the positives)
## forceprecrec: Continue even if files have been written
## printonly: if True, never execute commands
## nodefile: predicted nodes (if different from edgefile)
## nodesortcol: sort column for nodes
## descending: if True, rank in descending order
## param: additional string to append to outfile
## union: if True, use union of NetPath & KEGG as positives.
def computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,sampleoutprefix,\
                           subsamplefps,forceprecrec,printonly,nodefile=None,nodesortcol=None,\
                           descending=False,param=None,union=False):

    ## Get true edge file and true node file from the data directory.
    if union:
        ## TODO remove hard-coded links
        trueedgefile = '/data/annaritz/projects/2015-03-pathlinker/data/netpath-kegg-union/%s-edges.txt' % (pathway)
        truenodefile = '/data/annaritz/projects/2015-03-pathlinker/data/netpath-kegg-union/%s-nodes.txt' % (pathway)
        if not os.path.isfile(trueedgefile):
            print 'Skipping %s exclude %s for netpath-kegg-union evaultion: KEGG pathway doesn\'t exist. Skipping.' % (pathway,negtype)
            return
    else:
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

    ## if negtype == 'file', add kegg ignored nodes and edges
    if negtype == 'file':
        pathwaynames,kegg2netpath = getKEGGPathways() 
        keggid = [k for k in kegg2netpath if kegg2netpath[k]==pathway]
        if len(keggid)==0:
            print 'Skipping %s exclude %s: KEGG pathway does not exist.' % (pathway,negtype)
            return
        keggid = keggid[0]
        ignorededgefile = '%s/%s-edges.txt' % (KEGGDIR,keggid)
        ignorednodefile = '%s/%s-nodes.txt' % (KEGGDIR,keggid)

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
    if negtype == 'file':
        cmd += ' --ignorededgefile %s --ignorednodefile %s' % (ignorededgefile,ignorednodefile)
    #if opts.forceprecrec: # this will cause a resampling of the negatives.  USed for debugging.
    #    cmd += ' --force'
    print cmd
    if not printonly:
        subprocess.check_call(cmd.split())
    return

############################################################
## Computes aggregate precision and recall
## inputdir: directory of precision-recall files for individual pathways
## negtype: one of 'none' or 'adjacent'
## subsamplefps: Number of negatives to sample (a factor of the size of the positives)
## forceprecrec: Continue even if files have been written
## printonly: if True, never execute commands
## descending: if True, rank in descending order
## param: additional string to append to outfile
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


############################################################
if __name__=='__main__':
    main(sys.argv)
