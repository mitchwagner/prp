#!/usr/bin/python
"""
Master script for producing PathLinker comparative analyses.
"""

from optparse import OptionParser,OptionGroup
from utilsPoirel import *
import sys
import os
import os.path
import subprocess
import glob
import time
from distutils.dir_util import mkpath

##############################################################
# GLOBAL VARIABLES

# PPIVERION is going to be populated in main() with one of the ALLOWEDVERSIONS.
# We have settled on using "pathlinker-signalingg-children-reg", which is the
# 2015pathlinker interactome weighted using the 8 children of the "signal
# transduction" GO term as well as "regulation of signal transduction." All
# other versions are commented out.
PPIVERSION = ''
ALLOWEDVERSIONS = [
    #'2013linker', # Chris's version w/ vinayagam, etc.\
    #'2015pathlinker', # new version, reg weighting\
    #'pathlinker-old-kegg-buggy', # debug version\
    #'pathlinker-old-kegg-fixed', # debug version\
    #'pathlinker-old-kegg-new-parse', # debug version\
    #'pathlinker-old-kegg-buggy-oldannotations', # debug version\
    #'pathlinker-signaling',  # only weighted by sig-trans     \
    #'pathlinker-signaling-children',  # weighting\
    'pathlinker-signaling-children-reg',  # children and reg \
    'pathlinker-signaling-children-reg-no-netpath', # pathlinker-signaling-children-reg with no NetPath sources of evidence.
]

# DATADIR is the path to the data/ directory checked into SVN.  
# The interactomes, KEGG, and NetPath edge files are all checked in.
DATADIR = 'inputs/'

# Directory to write results to
OUTDIR = 'outputs/'

# ORIGINALPPI is the name of the original interactome, before making the
# pathway specific interactomes. THis is used when posting to GraphSpace so the
# directionality of the edges reflects the original interactome.  It is set in
# the main method.
ORIGINALPPI = ''

# PPIDIR is the directory that contains the pathway-specific interactomes.  For
# each NetPath pathway, we take the PPIVERSION network and remove incoming
# edges to receptors and outgoing edges from tfs.  Thus, there is one
# interactome for each NetPath pathway in the PPIDIR directory.  Note that
# "pathway-specific" is a bit too strong of a phrase, and shouldn't be used in
# the paper. 
PPIDIR = 'outputs/pathway-specific-interactomes/'

# MAPPING VARIABLES
# There is a human-gene-map.txt (originally from Chris) that is checked into
# SVN.  This is used to map between UniProtKB and the common name.  

# TODO: Once CSBDB is bug-free, we should use CSBDB to do this mapping.
# However that will require changing many of the auxilliary scripts.
MAPPINGFILE = '%s/namespace-mappers/human-gene-map.txt' % (DATADIR)
MAPFROMCOL = 6
MAPTOCOL = 1

# PATHWAY DIRECTORIES
# NetPath and KEGG edge files (*-edges.txt) are checked into SVN.
NETPATHDIR = \
    '%s/interactions/netpath/pathways/' % (DATADIR) 

KEGGDIR = \
    '%s/interactions/kegg/2015-03-23/hsa/edge-files-deadends/' % (DATADIR)

ORIGKEGGDIR = \
    '%s/interactions/kegg/2015-03-23/hsa/edge-files/' % (DATADIR)

PPIFILEDIR = \
    "%s/interactomes/human/" % (DATADIR)

# VARYPARAMS
# To select parameters for competing methods, we vary each parameter.
# The VARYPARAMS dictionary provides the parameter list.
VARYPARAMS = {'q': [0.1, 0.25, 0.5, 0.75, 0.9],  # PageRank teleportation probability
              'gamma':[10, 15, 20, 25, 30],      # ResponseNet parameter
              'omega':[0, 0.01, 0.1],            # PCSF penalty for adding trees
              'prize':[1, 3, 5, 7, 9],           # PCSF prize
              'alpha':[0, 0.1, 0.25, 0.4, 0.5],  # ANAT parameter
              'nmax':[5,10,15, 25, 35, 50, 75, 100, 200, 500], # IPA parameter
          }

def main(args):
    global PPIVERSION, PPIDIR, ORIGINALPPI

    opts = parseArguments(args)

    PPIVERSION = opts.ppiversion

    # BACKGROUND INTERACTOME
    # set PPI, determined by whether the --weightedppi argument is specified.
    # If --weightedppi, then the ppifile contains scores/probabilities 
    # (*-weighted.txt).  If not --weightedppi, then ppifile contains
    # scores of 1 for all the edges.
    # We also set the result directory prefix, resultprefix, and the
    # directory of pathway-specific interactomes (PPIDIR)
    if opts.weightedppi:
        ppifile = '%s/%s-weighted.txt' % (PPIFILEDIR, PPIVERSION)
        resultprefix = 'outputs/%s/weighted/' % (PPIVERSION)
        PPIDIR = '%s/%s/weighted/' % (PPIDIR,PPIVERSION)

    # TODO: This might not actually work...
    else:
        ppifile = '%s/%s.txt' % (PPIFILEDIR, PPIVERSION)
        resultprefix = 'outputs/%s/unweighted/' % (PPIVERSION)
        PPIDIR = '%s/%s/unweighted/' % (PPIDIR,PPIVERSION)
    ORIGINALPPI = ppifile

    # Make sure the directories exist; if they don't, create them.
    mkpath(resultprefix)
    mkpath(PPIDIR)

    # PATHWAY-SPECIFIC INTERACTOMES
    generatePathwaySpecificInteractomes(ppifile, PPIDIR, KEGGDIR, NETPATHDIR)

    # DATASETS
    # pathways is a set() of (pathwayname,resultdir,datadir,ppidir) tuples.
    pathways = set()
    kegg2netpath = {} # will be populated if kegg pathways are specified.
    if opts.netpath or opts.onlynetpathwnt or opts.allnetpath:
        pathwaynames = getNetPathPathways(
            opts.onlynetpathwnt, 
            opts.ignorekeggpositives or opts.netpathkeggunion,
            opts.allnetpath)

        resultdir = '%s/netpath/' %(resultprefix)
        mkpath(resultdir) 
        datadir = NETPATHDIR
        ppidir = '%s/netpath/' % (PPIDIR)
        pathways.update([(p,resultdir,datadir,ppidir) for p in pathwaynames])
        print 'Running %d NetPath pathways' % (len(pathwaynames))

    if opts.wntforexperiments:
        resultdir = '%s/wnt-all-receptors/' % (resultprefix)
        mkpath(resultdir)
        datadir = 'inputs/pathways/wnt-all-receptors/'
        ppidir = '%s/wnt-all-receptors/' % (PPIDIR)
        pathways.update([('Wnt', resultdir, datadir, ppidir)])
        print 'Running 1 wnt-all-receptors pathway'

    if opts.kegg:
        pathwaynames,kegg2netpath = getKEGGPathways(
            opts.ignorenetpathpositives or opts.netpathkeggunion) 

        resultdir = '%s/kegg/' % (resultprefix)
        mkpath(resultdir) 
        datadir = KEGGDIR
        ppidir = '%s/kegg/'% (PPIDIR)
        pathways.update([(p,resultdir,datadir,ppidir) for p in pathwaynames])
        print 'Running %d KEGG pathways' % (len(pathwaynames))
    if len(pathways)==0:
        print 'WARNING: no datasets specified. This is OK when visualizing output or running special runs.'
    else:
        print '%d total pathways considered.' % (len(pathways))


    # TODO: Flagging for deletion in the future. We should add proper
    # parallelization support, which should be trivial if we move to a better
    # structure
    # reverse the order of the pathways. Useful for running the same algorithm
    # on two different computers.
    if opts.rev_pathways: pathways = list(pathways)[::-1]

    # TODO: This text will need to change, if it is even included after I am
    # done with this pipeline
    # ALGORITHMS ##
    # For each algorithm, iterate through pathways and call the run() method.
    # If --varyparams is specified, iterate through (pathway,param) combinations
    # and cal the run() method.

    algorithms_to_run = get_algorithm_object_list_from_opts(opts)

    # Okay, so it looks like the old code was creating tuples that had the
    # "correct" output directory for each pathway, whether it was in
    # kegg or it was special wnt or it was in netpath.
    for (pathway, resultdir, datadir, ppidir) in pathways:
        file_location_context = FileLocationContext(
            pathway, datadir, ppidir, resultdir)

        for algorithm in algorithms_to_run:
            algorithm.run_if_forced(file_location_context, opts.forcealg)

    # VIZ SCRIPTS #
    # These are a "grab bag" of precision-recall computations, 
    # precision-recall plots, and other types of analyses.

    # write precision/recall
    if opts.computeprecrec:
        print 'Write Precision Recall'

        # Always run 'exclude-none' and 'exclude-adjacent'.
        # If --ignorekeggpositives, also pass KEGG nodes and
        # edges as ignored files along with the 'exclude-file'
        # parameter.  Add 'file' to the negtypes.
        negtypes = ['none','adjacent']
        if opts.ignorekeggpositives or opts.ignorenetpathpositives:
            negtypes.append('file')

        # Compute precision and recall for each negtype described above.
        for negtype in negtypes:

            # All algorithms use the same set of positives and subsampled
            # negatives for the dataset specified.  The sampledir variable
            # points to the directory, depending on the PPIVERSION and the
            # negtype.  If --netpathkeggunion is specified, this is a new set
            # of positives (union of NetPath and KEGG).  Thus, the sample
            # directory is prepended with "netpathkeggunion".
            if opts.netpathkeggunion:
                sampledir = '%s/netpathkeggunion-samples-exclude-%s' % (resultprefix,negtype)
            else:
                sampledir = '%s/samples-exclude-%s' % (resultprefix,negtype)
            mkpath(sampledir)

            # If --wntforexperimets is specified, then there is one dataset in
            # pathways variable that corresponds to wnt-all-receptors.  Since
            # there is a different set of positives than in NetPath Wnt pathway
            # (e.g., FZD4/FZD6 are added), we need to subsample a different set
            # of negatives.  Thus, the wntsampledir directory is prepended with
            # "wntforexperiments"
            if opts.wntforexperiments:
                wntsampledir = '%s/wntforexperiments-samples-exclude-%s' % (resultprefix,negtype)
                mkpath(wntsampledir)
            else:
                wntsampledir = ''
            
            # For every algorithm, do the following:
            # (1) determine the columns to sort on
            # (2) for each pathway:
            #   - get the file that contains the predicted edges (and possibly node file as well)
            #   - determine the output directory: it is either:
            #         <resultdir>/precision-recall/netpathkeggunion/<alg>/ if --netpathkeggunion is specified
            #         <resultdir>/precision-recall/<alg>/ otherwise.
            #   - determine the prefix of the subsampled negative files. It is either:
            #         <wntsampledir>/<pathway> if --wntforexperiments is specified
            #         <sampledir>/<pathway> otherwise
            #   - compute precision and recall
            # (3) if --netpath is specified, also compute aggregate precision and recall
            #
            # If --varyparams is specified, executes (2) and (3) for every
            # parameter value.

            # PATHLINKER #
            if opts.pathlinker:
                sortcol = 3 # ksp
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/pathlinker/%s-k_%s-ranked-edges.txt' % (resultdir,pathway,opts.k)
                    outdir = getPRoutdir('pathlinker',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                           union=opts.netpathkeggunion)
                if opts.netpath:
                    # compute aggregate for NetPath
                    inputdir = getPRoutdir('pathlinker',resultprefix+'/netpath/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion)
		if opts.allnetpath:
		    # comput aggregate for all of Netpath
		    inputdir = getPRoutdir('pathlinker',resultprefix+'/netpath/',opts.netpathkeggunion)
		    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion,allpathways=True)
                if opts.kegg:
                    # compute aggregate for KEGG
                    inputdir = getPRoutdir('pathlinker',resultprefix+'/kegg/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,netpath=False,\
                                                    union=opts.netpathkeggunion)
                    
            # PATHLINKER No Div #
            if opts.pathlinker_no_div:
                sortcol = 3 # ksp
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/pathlinker-no-div/%s-k_%s-ranked-edges.txt' % (resultdir,pathway,opts.k)
                    outdir = getPRoutdir('pathlinker-no-div',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                           union=opts.netpathkeggunion)
                if opts.netpath:
                    # compute aggregate for NetPath
                    inputdir = getPRoutdir('pathlinker-no-div',resultprefix+'/netpath/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion)
		if opts.allnetpath:
		    # comput aggregate for all of Netpath
		    inputdir = getPRoutdir('pathlinker-no-div',resultprefix+'/netpath/',opts.netpathkeggunion)
		    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion,allpathways=True)
                if opts.kegg:
                    # compute aggregate for KEGG
                    inputdir = getPRoutdir('pathlinker-no-div',resultprefix+'/kegg/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,netpath=False,\
                                                    union=opts.netpathkeggunion)
                    
            # PAGERANK PATHLINKER #
            if opts.pagerank_pathlinker:
                sortcol = 3 # ksp
                if not opts.varyparams: # just run with opts.q value
                    params = ['q_%.2f' % (opts.q)]
                else: # vary opts.q value
                    params = ['q_%.2f' % p for p in VARYPARAMS['q']]
                for param in params:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        edgefile = '%s/pagerank-pathlinker/%s-q_%.2fk_%s-ranked-edges.txt' % (resultdir,pathway,opts.q,opts.k)
                        outdir = getPRoutdir('pagerank-pathlinker',resultdir,opts.netpathkeggunion)
                        sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                               param=param,\
                                               union=opts.netpathkeggunion)
                    if opts.netpath:
                        # compute aggregate for NetPath
                        inputdir = getPRoutdir('pagerank-pathlinker',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,\
                                                        union=opts.netpathkeggunion)
                    if opts.allnetpath:
                        # comput aggregate for all of Netpath
                        inputdir = getPRoutdir('pagerank-pathlinker',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,\
                                                        union=opts.netpathkeggunion,allpathways=True)
                    if opts.kegg:
                        # compute aggregate for KEGG
                        inputdir = getPRoutdir('pagerank-pathlinker',resultprefix+'/kegg/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,netpath=False,\
                                                        param=param,\
                                                        union=opts.netpathkeggunion)
                    
            # CycLinker #
            if opts.cyclinker:
                sortcol = 3 # ksp
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/cyclinker/%s-ranked-edges.txt' % (resultdir,pathway)
                    outdir = getPRoutdir('cyclinker',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                           union=opts.netpathkeggunion)
                if opts.netpath:
                    # compute aggregate for NetPath
                    inputdir = getPRoutdir('cyclinker',resultprefix+'/netpath/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion)
		if opts.allnetpath:
		    # comput aggregate for all of Netpath
		    inputdir = getPRoutdir('cyclinker',resultprefix+'/netpath/',opts.netpathkeggunion)
		    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion,allpathways=True)
                if opts.kegg:
                    # compute aggregate for KEGG
                    inputdir = getPRoutdir('cyclinker',resultprefix+'/kegg/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,netpath=False,\
                                                    union=opts.netpathkeggunion)
            # PAGERANK CYCLINKER #
            if opts.pagerank_cyclinker:
                sortcol = 3 # ksp
                if not opts.varyparams: # just run with opts.q value
                    params = ['q_%.2f' % (opts.q)]
                else: # vary opts.q value
                    params = ['q_%.2f' % p for p in VARYPARAMS['q']]
                for param in params:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        edgefile = '%s/pagerank-cyclinker/%s-q_%.2f-ranked-edges.txt' % (resultdir,pathway,opts.q)
                        outdir = getPRoutdir('pagerank-cyclinker',resultdir,opts.netpathkeggunion)
                        sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                               param=param,\
                                               union=opts.netpathkeggunion)
                    if opts.netpath:
                        # compute aggregate for NetPath
                        inputdir = getPRoutdir('pagerank-cyclinker',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,\
                                                        union=opts.netpathkeggunion)
                    if opts.allnetpath:
                        # comput aggregate for all of Netpath
                        inputdir = getPRoutdir('pagerank-cyclinker',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,\
                                                        union=opts.netpathkeggunion,allpathways=True)
                    if opts.kegg:
                        # compute aggregate for KEGG
                        inputdir = getPRoutdir('pagerank-cyclinker',resultprefix+'/kegg/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,netpath=False,\
                                                        param=param,\
                                                        union=opts.netpathkeggunion)
                    
                    
            # Shortest Paths #
            if opts.shortestpaths:
                # no sorting; entire file is shortest paths subgraph.
                sortcol = None 
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/shortestpaths/%s-shortest-paths.txt' % (resultdir,pathway)
                    outdir = getPRoutdir('shortestpaths',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                           union=opts.netpathkeggunion)
                if opts.netpath:
                    # compute aggregate for NetPath
                    inputdir = getPRoutdir('shortestpaths',resultprefix+'/netpath/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion)
		if opts.allnetpath:
		    # comput aggregate for all of Netpath
		    inputdir = getPRoutdir('shortestpaths',resultprefix+'/netpath/',opts.netpathkeggunion)
		    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion,allpathways=True)
                if opts.kegg:
                    # compute aggregate for KEGG
                    inputdir = getPRoutdir('shortestpaths',resultprefix+'/kegg/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,netpath=False,\
                                                    union=opts.netpathkeggunion)
            # BowTieBUilder #                                                                                
            if opts.bowtiebuilder:

                # no sorting; entire file is the set of edges returned by
                # BowTieBuilder.
                sortcol = None 

                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/bowtiebuilder/%s-bowtiebuilder.txt' % (resultdir,pathway)
                    outdir = getPRoutdir('bowtiebuilder',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                           union=opts.netpathkeggunion)
                if opts.netpath:
                    # compute aggregate for NetPath                                                                                                  
                    inputdir = getPRoutdir('bowtiebuilder',resultprefix+'/netpath/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion)
		if opts.allnetpath:
		    # comput aggregate for all of Netpath
		    inputdir = getPRoutdir('bowtiebuilder',resultprefix+'/netpath/',opts.netpathkeggunion)
		    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion,allpathways=True)
                if opts.kegg:
                    # compute aggregate for KEGG                                                                                                  
                    inputdir = getPRoutdir('bowtiebuilder',resultprefix+'/kegg/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,netpath=False,\
                                                    union=opts.netpathkeggunion)
            # Induced Subgraph #
            if opts.inducedsubgraph:
                sortcol = 3 # ksp
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/inducedsubgraph/%s-induced-subgraph.txt' % (resultdir,pathway)
                    outdir = getPRoutdir('inducedsubgraph',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                           opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                           union=opts.netpathkeggunion)
                if opts.netpath:
                    # compute aggregate for NetPath
                    inputdir = getPRoutdir('inducedsubgraph',resultprefix+'/netpath/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion)
		if opts.allnetpath:
		    # comput aggregate for all of Netpath
		    inputdir = getPRoutdir('inducedsubgraph',resultprefix+'/netpath/',opts.netpathkeggunion)
		    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,\
                                                    union=opts.netpathkeggunion,allpathways=True)
                if opts.kegg:
                    # compute aggregate for KEGG
                    inputdir = getPRoutdir('inducedsubgraph',resultprefix+'/kegg/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,netpath=False,\
                                                    union=opts.netpathkeggunion)

            # PageRank #
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
                        sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,\
                                               opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,nodefile=nodefile,nodesortcol=nodesortcol,\
                                               descending=True,param=param,\
                                               union=opts.netpathkeggunion)
                    if opts.netpath:
                        # compute aggregate for NetPath
                        inputdir = getPRoutdir('pagerank',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,descending=True,\
                                                        param=param,\
                                                    union=opts.netpathkeggunion)

                    if opts.allnetpath:
                        # comput aggregate for all of Netpath
                        inputdir = getPRoutdir('pagerank',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,descending=True,\
                                                    param=param,\
                                                    union=opts.netpathkeggunion,allpathways=True)

                    if opts.kegg:
                        # compute aggregate for KEGG
                        inputdir = getPRoutdir('pagerank',resultprefix+'/kegg/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,\
                                                        opts.subsamplefps,opts.forceprecrec,opts.printonly,\
                                                        descending=True,param=param,netpath=False,\
                                                    union=opts.netpathkeggunion)

            # eQED #
            if opts.eqed:
                # both sort columns are in decreasing order.
                edgesortcol = 4 # abs(current)
                nodesortcol = 3 # positive input current
                for (pathway,resultdir,datadir,ppidir) in pathways:
                    edgefile = '%s/eqed/%s-eqed-edges.out' % (resultdir,pathway)
                    nodefile = '%s/eqed/%s-eqed-nodes.out' % (resultdir,pathway)
                    outdir = getPRoutdir('eqed',resultdir,opts.netpathkeggunion)
                    sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                    computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,\
                                           opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                           sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                           opts.printonly,nodefile=nodefile,nodesortcol=nodesortcol,descending=True,\
                                           union=opts.netpathkeggunion)
                if opts.netpath:
                    # compute aggregate for KEGG
                    inputdir = getPRoutdir('eqed',resultprefix+'/netpath/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,descending=True,\
                                                    union=opts.netpathkeggunion)

		if opts.allnetpath:
		    # comput aggregate for all of Netpath
		    inputdir = getPRoutdir('eqed',resultprefix+'/netpath/',opts.netpathkeggunion)
		    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,descending=True,\
                                                    union=opts.netpathkeggunion,allpathways=True)
                if opts.kegg:
                    # compute aggregate for KEGG
                    inputdir = getPRoutdir('eqed',resultprefix+'/kegg/',opts.netpathkeggunion)
                    computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                    opts.ignorenetpathpositives,opts.subsamplefps,\
                                                    opts.forceprecrec,opts.printonly,descending=True,netpath=False,\
                                                    union=opts.netpathkeggunion)

            # ResponseNet #
            if opts.responsenet:
                # Take union of edges/nodes with positive flow as
                # a subgraph.  both sort columns are in decreasing order.
                edgesortcol = None 
                if not opts.varyparams: # just run with opts.gamma value
                    params = ['gamma_%d' % (opts.gamma)]
                else: # vary opts.gamma value
                    params = ['gamma_%d' % p for p in VARYPARAMS['gamma']]
                for param in params:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        edgefile = '%s/reponsenet/%s-%s_responsenet-edges.out' % (resultdir,pathway,param)
                        outdir = getPRoutdir('responsenet',resultdir,opts.netpathkeggunion)
                        sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,edgesortcol,negtype,\
                                               opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,param=param,\
                                               union=opts.netpathkeggunion)
                    if opts.netpath:
                        # compute aggregate for NetPath
                        inputdir = getPRoutdir('responsenet',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,\
                                                        union=opts.netpathkeggunion)

                    if opts.allnetpath:
                        # comput aggregate for all of Netpath
                        inputdir = getPRoutdir('responsenet',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,\
                                                        union=opts.netpathkeggunion,allpathways=True)
                    if opts.kegg:
                        # compute aggregate for KEGG
                        inputdir = getPRoutdir('responsenet',resultprefix+'/kegg/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,netpath=False,\
                                                        union=opts.netpathkeggunion)
            # PCSF #
            if opts.pcsf:
                # no sorting; entire file is shortest paths subgraph.
                sortcol = None 

                # just run with opts.prize and opts.omega values
                if not opts.varyparams: 
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
                        sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,param=param,\
                                               union=opts.netpathkeggunion)
                    if opts.netpath:
                        # compute aggregate for NetPath
                        inputdir = getPRoutdir('pcsf',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,\
                                                    union=opts.netpathkeggunion)
                    if opts.allnetpath:
                        # comput aggregate for all of Netpath
                        inputdir = getPRoutdir('pcsf',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,\
                                                        union=opts.netpathkeggunion,allpathways=True)

                    if opts.kegg:
                        # compute aggregate for KEGG
                        inputdir = getPRoutdir('pcsf',resultprefix+'/kegg/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,netpath=False,\
                                                    union=opts.netpathkeggunion)
            # ANAT #
            if opts.anat:
                # no sorting; entire file is shortest paths subgraph.
                sortcol = None 

                # just run with opts.alpha value
                if not opts.varyparams: 
                    params = ['alpha%.2f' % (opts.alpha)]

                else: # vary opts.alpha value
                    params = ['alpha%.2f' % p for p in VARYPARAMS['alpha']]

                for param in params:
                    for (pathway,resultdir,datadir,ppidir) in pathways:
                        edgefile = '%s/anat/%s-%s-edges.out' % (resultdir,pathway,param)
                        outdir = getPRoutdir('anat',resultdir,opts.netpathkeggunion)
                        sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,param=param,\
                                               union=opts.netpathkeggunion)
                    if opts.netpath:
                        # compute aggregate for NetPath
                        inputdir = getPRoutdir('anat',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,\
                                                    union=opts.netpathkeggunion)
                    if opts.allnetpath:
                        # comput aggregate for all of Netpath
                        inputdir = getPRoutdir('anat',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,\
                                                        union=opts.netpathkeggunion,allpathways=True)
                    if opts.kegg:
                        # compute aggregate for KEGG
                        inputdir = getPRoutdir('anat',resultprefix+'/kegg/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        param=param,netpath=False,\
                                                    union=opts.netpathkeggunion)
            # DEGREE #
            if opts.degree:
                sys.exit('ERROR: DEGREE is not implemented in this version. Exiting.')

            # IPA #
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
                        sampleoutprefix = getPRsubsampleprefix(resultdir,wntsampledir,sampledir,pathway)
                        computePrecisionRecall(pathway,datadir,ppidir,edgefile,outdir,sortcol,negtype,\
                                               opts.ignorekeggpositives,opts.ignorenetpathpositives,\
                                               sampleoutprefix,opts.subsamplefps,opts.forceprecrec,\
                                               opts.printonly,descending=True,param=param,\
                                               union=opts.netpathkeggunion)
                    if opts.netpath:
                        # compute aggregate for NetPath
                        inputdir = getPRoutdir('ipa',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        descending=True,param=param,\
                                                    union=opts.netpathkeggunion)

                    if opts.allnetpath:
                        # compute aggregate for NetPath
                        inputdir = getPRoutdir('ipa',resultprefix+'/netpath/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        descending=True,param=param,\
                                                    union=opts.netpathkeggunion,allpathways=True)

                    if opts.kegg:
                        # compute aggregate for KEGG
                        inputdir = getPRoutdir('ipa',resultprefix+'/kegg/',opts.netpathkeggunion)
                        computeAggregatePrecisionRecall(inputdir,negtype,opts.ignorekeggpositives,\
                                                        opts.ignorenetpathpositives,opts.subsamplefps,\
                                                        opts.forceprecrec,opts.printonly,\
                                                        descending=True,param=param,netpath=False,\
                                                    union=opts.netpathkeggunion)

    # Plot precision and recall, once the values have been computed.
    if opts.precrecviz:
        print 'Plot Precision Recall'
        # The algcmd variable contains all of the algorithms that have been
        # specified in the command line. These are the algorithms that will be
        # plotted in the precision-recall curves.
        algcmd = ''
        if opts.pathlinker:
            algcmd += ' --alg pathlinker'
        if opts.pathlinker_no_div:
            algcmd += ' --alg pathlinker-no-div'
        if opts.pagerank_pathlinker:
            algcmd += ' --alg pagerank-pathlinker'
        if opts.cyclinker:
            algcmd += ' --alg cyclinker'
        if opts.pagerank_cyclinker:
            algcmd += ' --alg pagerank-cyclinker'
        if opts.pagerank:
            algcmd += ' --alg pagerank'
        if opts.shortestpaths:
            algcmd += ' --alg shortestpaths'
        if opts.bowtiebuilder:
            algcmd += ' --alg bowtiebuilder'
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

        # Speficy the input directory to read the precision-recall values from.
        # if --netpathkeggunion is specified, add the netpathkeggunion
        # directory so we look in the correct place.
        if opts.netpath or opts.onlynetpathwnt or opts.allnetpath:
            if opts.netpathkeggunion:
                indir = '%s/netpath/precision-recall/netpathkeggunion/' % (resultprefix)
            else:
                indir = '%s/netpath/precision-recall/' % (resultprefix)
        elif opts.kegg:
            if opts.netpathkeggunion:
                indir = '%s/kegg/precision-recall/netpathkeggunion/' % (resultprefix)
            else:
                indir = '%s/kegg/precision-recall/' % (resultprefix)

        # For each pathway, determine the output prefix, construct the call to 
        # plot-precision-recall.py, and execute it.
        for (pathway,resultdir,datadir,ppidir) in pathways:
            # if --wntforexperiments is specified, add "all-receptors" to
            # label. If --ignorekeggpositives is specified, add
            # "ignorekeggpositives" to label.  if --netpathkeggunion is
            # specified, add "netpathkeggunion" to label.  Otherwise, label is
            # simply the pathway name.
            if 'netpath' in datadir:
                outprefix = 'viz/precision-recall/netpath/%s' % (pathway)
            else:
                outprefix = 'viz/precision-recall/kegg/%s' % (pathway)
            mkpath(os.path.dirname(outprefix))
            if opts.wntforexperiments and 'wnt-all-receptors' in resultdir:
                outprefix += '-all-receptors'
            elif opts.ignorekeggpositives:
                outprefix += '-ignorekeggpositives' 
            elif opts.ignorenetpathpositives:
                outprefix += '-ignorenetpathpositives' 
            elif opts.netpathkeggunion:
                outprefix += '-netpathkeggunion'
	    if opts.varyparams:
                outprefix += '-varyparams'
                if opts.pcsf and not opts.pathlinker:
                    outprefix += '-%s'%('_'.join([str(o) for o in VARYPARAMS['omega']]))
                    outprefix += '-%s'%('_'.join([str(p) for p in VARYPARAMS['prize']]))
            if PPIVERSION == 'pathlinker-signaling-children-reg-no-netpath':
                outprefix += '-no_netpath'
            if opts.forceviz or not os.path.isfile('%s.pdf' % (outprefix)):
                # Consruct the command.  
                cmd = 'python src/plot-precision-recall.py --indir %s --outprefix %s --pathway %s %s --pdf' % \
                      (indir,outprefix,pathway,algcmd)
                if opts.ignorekeggpositives:
                    cmd += ' --ignorekegg'
                if opts.ignorenetpathpositives:
                    cmd += ' --ignorenetpath'

                # Only plot the varying parameters for Wnt.  This reduces the
                # amount of clutter in the viz/ directory.
                if opts.varyparams:
                    if pathway != 'Wnt':
                        print 'WARNING: not plotting non-Wnt pathway %s with varying parameters' % (pathway)
                    else:
                        cmd += '  --varyparams'
                        if opts.pcsf and not opts.pathlinker:
                            cmd += ' --varyomega %s'%('_'.join([str(o) for o in VARYPARAMS['omega']]))
                            cmd += ' --varyprize %s'%('_'.join([str(p) for p in VARYPARAMS['prize']]))
                        print cmd
                        if not opts.printonly:
                            subprocess.check_call(cmd.split())
                else:
                    print cmd
                    if not opts.printonly:
                        subprocess.check_call(cmd.split())
            else:
                print '%s.pdf exists; not overwriting. Use --forceviz to override.' % (outprefix)
                
        # If --netpath or --allnetpath is specified, plot the aggregate
        # precision-recall plots.
        if opts.netpath or opts.allnetpath: 
            outprefix = 'viz/precision-recall/netpath/aggregate'
            if opts.ignorekeggpositives:
                outprefix += '-ignorekeggpositives' 
            elif opts.netpathkeggunion:
                outprefix += '-netpathkeggunion'
	    elif opts.allnetpath:
	        outprefix += '-allnetpath'   
            if opts.varyparams:
                outprefix += '-varyparams'
                if opts.pcsf and not opts.pathlinker:
                    outprefix += '-%s'%('_'.join([str(o) for o in VARYPARAMS['omega']]))
                    outprefix += '-%s'%('_'.join([str(p) for p in VARYPARAMS['prize']]))
            if PPIVERSION == 'pathlinker-signaling-children-reg-no-netpath':
                outprefix += '-no_netpath'
            if opts.forceviz or not os.path.isfile('%s.pdf' % (outprefix)):
		if opts.allnetpath:
 		    cmd = 'python src/plot-precision-recall.py --indir %s --outprefix %s --pathway aggregate-allpathways %s --pdf' % \
                      (indir,outprefix,algcmd)
		else:
                    cmd = 'python src/plot-precision-recall.py --indir %s --outprefix %s --pathway aggregate %s --pdf' % \
                      (indir,outprefix,algcmd)
                if opts.varyparams:
                    cmd += '  --varyparams'
                    if opts.pcsf and not opts.pathlinker:
                        cmd += ' --varyomega %s'%('_'.join([str(o) for o in VARYPARAMS['omega']]))
                        cmd += ' --varyprize %s'%('_'.join([str(p) for p in VARYPARAMS['prize']]))
                if opts.ignorekeggpositives:
                    cmd += ' --ignorekegg'
                print cmd
                if not opts.printonly:
                    subprocess.check_call(cmd.split())
            else:
                print '%s.pdf exists; not overwriting. Use --forceviz to override.' % (outprefix)

        # If --kegg is specified, plot the aggregate precision-recall plots.
        if opts.kegg: 
            outprefix = 'viz/precision-recall/kegg/aggregate'
            if opts.ignorenetpathpositives:
                outprefix += '-ignorenetpathpositives' 
            elif opts.netpathkeggunion:
                outprefix += '-netpathkeggunion'        
            if opts.forceviz or not os.path.isfile('%s.pdf' % (outprefix)):
                cmd = 'python src/plot-precision-recall.py --indir %s --outprefix %s --pathway aggregate %s --pdf' % \
                      (indir,outprefix,algcmd)
                if opts.varyparams:
                    cmd += '  --varyparams'
                if opts.ignorenetpathpositives:
                    cmd += ' --ignorenetpath'
                print cmd
                if not opts.printonly:
                    subprocess.check_call(cmd.split())
            else:
                print '%s.pdf exists; not overwriting. Use --forceviz to override.' % (outprefix)


    # Post Wnt Pathways to Graphspace Post two different Wnt pathway runs.  -
    # Post NetPath Wnt pathways, used to compute precision and recall.  These
    # are prepended with "pr" - Post wnt-all-receptors pathways, used to
    # explore false positives for experimental followup New: only runs when
    # --nepath option is NOT specified. Otherwise all NetPath pathways are run
    # (see below)
    if (opts.graphspace or opts.oldgraphspace) and not opts.netpath:
        print 'Posting to GraphSpace...'
        outdir = 'viz/graphspace-json/'
        if opts.pathlinker:
            infile_allreceptors = '%s/wnt-all-receptors/pathlinker/Wnt-k_%d-ranked-edges.txt' % (resultprefix,opts.k)
            infile = '%s/netpath/pathlinker/Wnt-k_%d-ranked-edges.txt' % (resultprefix,opts.k)
            gsid = 'Wnt-pathlinker-top%dpaths' % (opts.topk)
            postReconstructionsToGraphSpace('Wnt',infile,infile_allreceptors,outdir,opts.topk,\
                                            gsid,opts.printonly,increase=True,oldgs=opts.oldgraphspace,posttag=opts.tag)

        if opts.cyclinker:
            infile_allreceptors = '%s/wnt-all-receptors/cyclinker/Wnt-ranked-edges.txt' % (resultprefix)
            infile = '%s/netpath/cyclinker/Wnt-ranked-edges.txt' % (resultprefix)
            gsid = 'Wnt-cyclinker-top%dpaths' % (opts.topk)
            postReconstructionsToGraphSpace('Wnt',infile,infile_allreceptors,outdir,opts.topk,\
                                            gsid,opts.printonly,increase=True,oldgs=opts.oldgraphspace,posttag=opts.tag)

        if opts.pagerank: # Manually-determined threshold
            # threshold is set to recall of 0.203
            for thres in [0.000825996]:
                infile_allreceptors = '%s/wnt-all-receptors/pagerank/Wnt-q_0.50-edge-fluxes.txt' % (resultprefix)
                infile = '%s/netpath/pagerank/Wnt-q_0.50-edge-fluxes.txt' % (resultprefix)
                gsid = 'Wnt-pagerank-thres%.5f' % (thres)
                postReconstructionsToGraphSpace('Wnt',infile,infile_allreceptors,outdir,thres,gsid,opts.printonly,decrease=True,oldgs=opts.oldgraphspace,posttag=opts.tag)
                
        if opts.anat:
            infile_allreceptors = '%s/wnt-all-receptors/anat/Wnt-alpha0.00-edges.out' % (resultprefix)
            infile = '%s/netpath/anat/Wnt-alpha0.00-edges.out' % (resultprefix)
            gsid = 'Wnt-anat-alpha0.00'
            postReconstructionsToGraphSpace('Wnt',infile,infile_allreceptors,outdir,None,gsid,opts.printonly,oldgs=opts.oldgraphspace,posttag=opts.tag)

        if opts.shortestpaths:
            infile_allreceptors = '%s/wnt-all-receptors/shortestpaths/Wnt-shortest-paths.txt' % (resultprefix)
            infile = '%s/netpath/shortestpaths/Wnt-shortest-paths.txt' % (resultprefix)
            gsid = 'Wnt-shortest-paths'
            postReconstructionsToGraphSpace('Wnt',infile,infile_allreceptors,outdir,None,gsid,opts.printonly,oldgs=opts.oldgraphspace,posttag=opts.tag)


        if opts.bowtiebuilder:
            infile_allreceptors = '%s/wnt-all-receptors/bowtiebuilder/Wnt-bowtiebuilder.txt' % (resultprefix)
            infile = '%s/netpath/bowtiebuilder/Wnt-bowtiebuilder.txt' % (resultprefix)
            gsid = 'Wnt-bowtiebuilder'
            postReconstructionsToGraphSpace('Wnt',infile,infile_allreceptors,outdir,None,gsid,opts.printonly,oldgs=opts.oldgraphspace,posttag=opts.tag)

        if opts.responsenet:
            infile_allreceptors = '%s/wnt-all-receptors/reponsenet/Wnt-gamma_20_responsenet-edges.out' % (resultprefix)
            infile = '%s/netpath/reponsenet/Wnt-gamma_20_responsenet-edges.out' % (resultprefix)
            gsid = 'Wnt-responsenet-gamma20'
            postReconstructionsToGraphSpace('Wnt',infile,infile_allreceptors,outdir,None,gsid,opts.printonly,oldgs=opts.oldgraphspace,posttag=opts.tag)

        if opts.pcsf:
            infile = '%s/wnt-all-receptors/pcsf/Wnt-prize5-omega0.01_PCSF-edges.out' % (resultprefix)
            infile_allreceptors = '%s/netpath/pcsf/Wnt-prize5-omega0.01_PCSF-edges.out' % (resultprefix)
            gsid = 'Wnt-pcsf-prize5-omega0.01'
            postReconstructionsToGraphSpace('Wnt',infile,infile_allreceptors,outdir,None,gsid,opts.printonly,oldgs=opts.oldgraphspace,posttag=opts.tag)

        if opts.ipa:
            infile_allreceptors = '%s/wnt-all-receptors/ipa/Wnt-nmax10.out' % (resultprefix)
            infile = '%s/netpath/ipa/Wnt-nmax10.out' % (resultprefix)
            gsid = 'Wnt-ipa-nmax10'
            postReconstructionsToGraphSpace('Wnt',infile,infile_allreceptors,outdir,None,gsid,opts.printonly,undirected=True,oldgs=opts.oldgraphspace,posttag=opts.tag)
    # Post NetPath Pathways:
    if (opts.graphspace or opts.oldgraphspace) and opts.netpath:
        print 'Posting to GraphSpace...'
        outdir = 'viz/graphspace-json/'
        if not opts.pathlinker:
            sys.exit("ERROR: when posting all NP pathways to GraphSpace, only PathLinker graphs are currently posted.")
        for (pathway,resultdir,datadir,ppidir) in pathways:
            infile = '%s/netpath/pathlinker/%s-k_%d-ranked-edges.txt' % (resultprefix,pathway,opts.k)
            gsid = 'NetPath-%s-pathlinker-top%dpaths' % (pathway,opts.topk)
            postNetPathReconstructionsToGraphSpace(pathway,infile,outdir,opts.topk,\
                                            gsid,opts.printonly,increase=True,oldgs=opts.oldgraphspace,posttag=opts.tag)
    # RANK TF
    # Little script that ranks the TRs in PathLinker predictions vs. PageRank
    # predictions Only plot Wnt and aggregate rankings.
    if opts.ranktfs:
        if opts.netpath:            
            indir = '%s/netpath/' % (resultprefix)
            nodesdir = '%s/interactions/netpath/pathways/' % (DATADIR)

            outprefix = 'viz/ranking-receptors-trs/netpath-Wnt'
            cmd = 'python src/plot-tr-ranks.py --datadir %s --o %s --indir %s --pathway Wnt' % (nodesdir,outprefix,indir)
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())

                outprefix = 'viz/ranking-receptors-trs/netpath-Wnt'
            cmd = 'python src/plot-tr-ranks.py --datadir %s --o %s --indir %s --pathway Wnt --truncate 400' % (nodesdir,outprefix,indir)
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())

            outprefix = 'viz/ranking-receptors-trs/netpath-aggregate'
            cmd = 'python src/plot-tr-ranks.py --datadir %s --o %s --indir %s --pathway aggregate' % (nodesdir,outprefix,indir)
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())

            outprefix = 'viz/ranking-receptors-trs/netpath-aggregate'
            cmd = 'python src/plot-tr-ranks.py --datadir %s --o %s --indir %s --pathway aggregate --truncate 1000' % (nodesdir,outprefix,indir)
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())

        if opts.kegg:            
            indir = '%s/kegg/' % (resultprefix)
            nodesdir = '%s/interactions/kegg/2015-03-23/hsa/edge-files/' % (DATADIR)

            outprefix = 'viz/ranking-receptors-trs/kegg-hsa04310'
            cmd = 'python src/plot-tr-ranks.py --datadir %s --o %s --indir %s --pathway hsa04310 --kegg' % (nodesdir,outprefix,indir)
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())

            outprefix = 'viz/ranking-receptors-trs/kegg-aggregate'
            cmd = 'python src/plot-tr-ranks.py --datadir %s --o %s --indir %s --pathway aggregate --kegg' % (nodesdir,outprefix,indir)
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())

            outprefix = 'viz/ranking-receptors-trs/kegg-aggregate'
            cmd = 'python src/plot-tr-ranks.py --datadir %s --o %s --indir %s --pathway aggregate --kegg --truncate 3000' % (nodesdir,outprefix,indir)
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())

    if opts.falsepos:
        print 'Plotting false positives'
        
        # compute distances
        for (pathway,resultdir,datadir,ppidir) in pathways:
            if opts.forceviz or not os.path.isfile('outputs/shortest-paths-for-false-positives/%s-edge-dists-undirected.txt' % (pathway)):
                cmd = 'python src/shortest-paths-for-false-positives.py --datadir %s --ppidir %s --pathway %s' % \
                      (datadir,ppidir,pathway)
                print cmd
                if not opts.printonly:
                    subprocess.check_call(cmd.split())
            else:
                print 'Shortest Paths from nodes to pathways for %s are already computed.' % (pathway)

        for (pathway,resultdir,datadir,ppidir) in pathways:
            if pathway != 'Wnt':
                print 'Warning: not plotting false positives for pathways other than Wnt.'
                continue

            indir = '%s/precision-recall/' % (resultdir)
            cmd = 'python src/plot-false-positive-distances.py --outprefix viz/false-positives/directed- --indir %s --pathway aggregate --pathway %s --pdf --edges -r 0.10 -r 0.20 -r 0.30 -r 0.40 -r 0.50 -r 0.60 -r 0.70  --alg pathlinker --alg pagerank' % (indir,pathway)
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())

            cmd = 'python src/plot-false-positive-distances.py --outprefix viz/false-positives/undirected- --indir %s --pathway aggregate --pathway %s --pdf --edges -r 0.10 -r 0.20 -r 0.30 -r 0.40 -r 0.50 -r 0.60 -r 0.70  --alg pathlinker --alg pagerank --undirected' % (indir,pathway)
            print cmd
            if not opts.printonly:
                subprocess.check_call(cmd.split())



    # make venn diagrams
    if opts.venn:
        cmd = 'python src/make-venn.py --datadir %s --inputdir %s' % (DATADIR,'%s/netpath/precision-recall/pathlinker/' % (resultprefix))
        print cmd
        if not opts.printonly:
            subprocess.check_call(cmd.split())

    # perform the subsampling analysis
    if opts.sample_sources_targets:
        performSubsampling(pathways, opts.k, opts.forcealg, opts.forceprecrec, opts.printonly, opts.batch_run, opts)

    print 'DONE'
    return

def parseArguments(args):
    """
    Parses command-line options
    """
    usage = 'master-script.py [options]\n'
    parser = OptionParser(usage=usage)

    # Common options
    parser.add_option('','--forcealg',action='store_true',default=False,\
                          help='Run algorithms even if they will overwrite existing files.  Default is that algorithms are not run if output files exist.')
    parser.add_option('','--forceprecrec',action='store_true',default=False,\
                      help='Run computations of precision and recall, even if they will overwrite existing files. Default is that precision and recall are not run if output files exist.')
    parser.add_option('','--forceviz',action='store_true',default=False,\
                      help='Run visualizations, even if they will overwrite existing files. Default is that some viz computations are not run if files exist..')
    parser.add_option('','--printonly',action='store_true',default=False,\
                          help='Print the commands to stdout, but do not execute.')

    # Datasets.  
    group = OptionGroup(parser,'Datasets')
    group.add_option('','--ppiversion',type='string',default='pathlinker-signaling-children-reg',\
                     help='Version of the PPI to run.  Options are %s. Default is "pathlinker-signaling-children-reg."' % (', '.join(ALLOWEDVERSIONS)))
    group.add_option('','--weightedppi',action='store_true',default=False,\
                         help='Run with weighted PPI.')
    group.add_option('','--onlynetpathwnt',action='store_true',default=False,\
                         help='Only run NetPath Wnt. Only one of --onlynetpathwnt,--netpath,--kegg,--wntforexperiments--allnetpath may be specified.')
    group.add_option('','--netpath',action='store_true',default=False,\
                         help='Run with NetPath inputs.  Only one of --onlynetpathwnt,--netpath,--kegg,--wntforexperiments,--allnetpath may be specified.')
    group.add_option('','--allnetpath',action='store_true',default=False,\
			help='Run with all NetPath pathways. Only one of --onlynetpathwnt, --netpath, --kegg, --wntforexpeirments, --allnetpath may be specified.')
    group.add_option('','--kegg',action='store_true',default=False,\
                         help='Run with KEGG inputs.  Only one of --onlynetpathwnt,--netpath,--kegg,--wntforexperiments,--allnetpath may be specified.')
    group.add_option('','--wntforexperiments',action='store_true',default=False,\
                     help='Run special wnt that includes FZD4/FZD6 receptors, for analyzing via networks.  Only one of --onlynetpathwnt,--netpath,--kegg, --wntforexperiments,--allnetpath may be specified')
    parser.add_option_group(group)

    # Algorithms
    group = OptionGroup(parser,'Algorithms')
    group.add_option('','--pathlinker',action='store_true',default=False,\
                          help='Run PathLinker (KSP) on input files.')
    group.add_option('','--pathlinker-no-div',action='store_true',default=False,\
                          help='Run PathLinker (KSP) without the large division of each edge weight on input files.')
    group.add_option('','--pagerank-pathlinker',action='store_true',default=False,\
                          help='Run PathLinker (KSP) by first re-weighting edges using PageRank and without the large division of each edge weight on input files.')
    group.add_option('','--cyclinker',action='store_true',default=False,\
                          help='Run CycLinker on input files.')
    group.add_option('','--pagerank-cyclinker',action='store_true',default=False,\
                          help='Run CycLinker by first re-weighting edges using PageRank on input files.')
    group.add_option('','--shortestpaths',action='store_true',default=False,\
                     help='Compute shortest paths from each receptor to each TF, incl. ties (reviewer comment)')
    group.add_option('','--bowtiebuilder',action='store_true',default=False,\
                     help='Run BowTieBuilder (reviewer comment)')
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

    # Robustness
    group = OptionGroup(parser,'Robustness')
    
    group.add_option('','--sample-sources-targets',action='store_true',default=False,\
                          help='Run a robustness analysis')
    group.add_option('','--only-generate-sample-sets',action='store_true',default=False,\
                          help='Generate the sample sets for a robsutness analysis.')
    group.add_option('','--batch-run',action='store_true',default=False,\
            help='Mark files so that many instances can be run in parallel. This means that force-killing the script can leave result files in a bad state. Quick instructions: To orchestrate a batch run from scratch, run once with --only-generate-sample-sets, then many times with --batch-run, then once with neither of these.')
    
    parser.add_option_group(group)

    # Precision Recall and Other Visualizations
    group = OptionGroup(parser,'Visualizations')
    group.add_option('','--computeprecrec',action='store_true',default=False,\
                         help='Compute precision and recall curves; write to file.')
    group.add_option('','--precrecviz',action='store_true',default=False,\
                         help='Display precision recall curves.')
    group.add_option('','--ignorekeggpositives',action='store_true',default=False,\
                     help='Ignore KEGG positives when computing & visualizing precision/recall (--netpath or --onlynetpathwnt options). Will ignore pathways that do not appear in KEGG.')
    group.add_option('','--ignorenetpathpositives',action='store_true',default=False,\
                     help='Ignore NetPath positives when computing & visualizing precision/recall (--kegg option). Will ignore pathways that do not appear in NetPath.')
    group.add_option('','--netpathkeggunion',action='store_true',default=False,\
                     help='Compute & visualize precision/recall using both NetPath and KEGG as positives.')
    group.add_option('','--ranktfs',action='store_true',default=False,\
                     help='Plot TF ranks for PageRank, PathLinker, and IPA.')
    group.add_option('','--falsepos',action='store_true',default=False,\
                         help='Visualize distances of false positives to true network.')
    group.add_option('','--falsenegs',action='store_true',default=False,\
                         help='Make false negative plots.')
    group.add_option('','--venn',action='store_true',default=False,\
                         help='plot venn diagrams. Only with dbpcompare.')
    group.add_option('','--graphspace',action='store_true',default=False,\
                         help='post Wnt networks to python (new) graphspace.')
    group.add_option('','--oldgraphspace',action='store_true',default=False,\
                         help='post Wnt networks to perl (old) graphspace.')
    group.add_option('','--paper',action='store_true',default=False,\
                         help='Make plots for paper.')
    group.add_option('','--weightviz',action='store_true',default=False,\
                         help='Post weighted vs. unweighted graph to GraphSpace')
    group.add_option('','--tag',action='store_true',default=False,\
                     help='Label graphspace graphs with public tag.')
    parser.add_option_group(group)

    # Optional Arguments for Algorithms
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
    group.add_option('','--prize',type='int',default=1,metavar='FLOAT',\
                         help='Prize for PCSF. Default is 1.')
    group.add_option('','--omega',type='float',default=0.01,metavar='FLOAT',\
                     help='Cost of adding additional trees to PCSF forest. Default is 0.01.')
    group.add_option('','--alpha',type='float',default=0.0,metavar='FLOAT',\
                     help='Tradeoff between shortest paths and steiner trees in ANAT. Default is 0.00.')
    group.add_option('','--nmax',type='int',default=10,metavar='INT',\
                     help='Size for IPA\'s network generation algorithm. Default is 10')
    group.add_option('','--subsamplefps',type='int',default=50,metavar='INT',\
                         help='Subsample so the number of FPs is this many times the # of TPs for each pathway.  If this value is -1, then all negatives are used.')
    group.add_option('','--bullseyerecall',type='float',default=0.5,metavar='FLOAT',\
			 help='Recall for bullseye plots.  Default is 0.5.')
    group.add_option('','--rev-pathways',action='store_true',default=False,\
			 help='Reverse the order of the pathways when running the algorithms. Useful for running the same algorithm on two different machines')
    parser.add_option_group(group)
    
    # parse the command line arguments
    (opts, args) = parser.parse_args()

    # check arguments
    if opts.ppiversion not in ALLOWEDVERSIONS:
        sys.exit('ERROR: --ppiversion must be one of %s. You have %s. Exiting.' % (','.join(ALLOWEDVERSIONS),opts.ppiversion))
    if opts.ignorekeggpositives and opts.ignorenetpathpositives: 
        sys.exit('ERROR: cannot ignore both KEGG positives and NetPath positives. Exiting.')

    if sum([x for x in [opts.netpath,opts.kegg,opts.onlynetpathwnt,opts.wntforexperiments,opts.allnetpath]]) > 1:
        sys.exit('ERROR: only one of --netpath, --onlynetpathwnt, --kegg, or --wntforexperiments, --allnetpath may be specified. Exiting.')

    return opts


def generatePathwaySpecificInteractomes(
        ppifile, ppidir, keggdir, netpathdir):
    """
    For each pathway in KEGG, NetPath, and the wnt-all-receptors 
    datasets, remove the incoming edges to receptors and outgoing
    edges to TRs and save them in the PPIDIR

    :param ppifile: original PPI file 
        (e.g., pathlinker-signaling-children-reg-weighted.txt)

    :param ppidir: Directory containing pathway-specific interactomes
    :param keggdir: Directory containing KEGG edge files
    :param netpathdir: Directory containing NETPATH edge files

    :return: None

    """

    # Read original PPI file 
    edges = [] # list of (u,v,line) tuples
    header = ''
    with open(ppifile) as fin:
        for line in fin:
            if line[0]=='#':
                header = line
                continue
            row = line.split('\t')
            edges.append((row[0],row[1],line))

    # Make NetPath interactomes, if not already present
    mkpath(ppidir+'/netpath/')
    pathways = getAllNetPathPathways()
    for p in pathways:
        interactomefile = '%s/netpath/%s-interactome.txt' % (ppidir, p)
        if not os.path.isfile(interactomefile):
            print 'Making NetPath %s Interactome' % (p)
            nodefile = '%s/%s-nodes.txt' % (netpathdir,p)
            generatePPI(edges, nodefile, interactomefile, header)

    callComputeMinCut(
        NETPATHDIR, PPIDIR + '/netpath/', 'outputs/min-cuts/netpath.txt')

    # Make wnt-all-receptors interactome, if not already present.
    mkpath(ppidir+'/wnt-all-receptors/')
    interactomefile = '%s/wnt-all-receptors/Wnt-interactome.txt' % (ppidir)
    if not os.path.isfile(interactomefile):
        print 'Making Wnt All receptors Interactome'
        nodefile = 'inputs/pathways/wnt-all-receptors/Wnt-nodes.txt'
        generatePPI(edges, nodefile, interactomefile, header)

    # Make KEGG interactomes, if not already present.
    mkpath(ppidir +'/kegg/')
    pathways = getAllKEGGPathways()
    for p in pathways:
        interactomefile = '%s/kegg/%s-interactome.txt' % (ppidir, p)
        if not os.path.isfile(interactomefile):
            print 'Making KEGG %s Interactome' % (p)
            nodefile = '%s/%s-nodes.txt' % (keggdir,p)
            generatePPI(edges,nodefile,interactomefile,header)

    callComputeMinCut(
        ORIGKEGGDIR, PPIDIR + 'kegg/', 'outputs/min-cuts/kegg.txt', 
        DATADIR + '/interactions/kegg/2015-03-23/hsa/HSA_PATHWAY_LIST_FORMATTED.txt')
        
    callComputeMinCut(
        KEGGDIR, PPIDIR + 'kegg/', 'outputs/min-cuts/kegg-deadends.txt', 
        DATADIR + '/interactions/kegg/2015-03-23/hsa/HSA_PATHWAY_LIST_FORMATTED.txt')

    return

def callComputeMinCut(datadir, ppidir, outfile, mapfile=None):
    mkpath(os.path.dirname(outfile))
    cmd = buildComputeMinCutCommand(datadir, ppidir, outfile, mapfile)
    print("Running: ", cmd)
    subprocess.call(cmd)
    
def buildComputeMinCutCommand(datadir, ppidir, outfile, mapfile=None): 
    cmd = ['python', 'src/compute-min-cut.py', 
           '--datadir', datadir, 
           '--ppidir', ppidir, 
           '--outfile', outfile]
    if mapfile:
        cmd = cmd + ['--mapfile', mapfile]
    
    return cmd

def generatePPI(edges, nodefile, interactomefile, header):
    """
    generatePPI makes the pathway-specific interactome.

    :param edges: Original PPI edges
    :param nodefile: File of nodes (to get receptors and TRs)
    :param interactomefile: Name of ouptut file
    :param header: Common header to put at the top of the output file

    :return: None

    """
    # get receptors and tfs
    nodes = readColumns(nodefile,1,2)
    receptors = set([n for n,t in nodes if t == 'receptor'])
    tfs = set([n for n,t in nodes if t == 'tf'])    
    
    # Write edges to file if they are not incoming
    # to receptors or outgoing from TFs.
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

    return


# TODO: This and the following function look to have a good
# deal of overlap
def getNetPathPathways(onlynetpathwnt, overlapwithkegg, allnetpath):
    """
    Reads the pathways for NetPath and returns them as a list.

    :param onlynetpathwnt: If True, only return Wnt

    :param dbcompare: If True, only return the 6 pathways in common
        with KEGG

    """
    if onlynetpathwnt:
        return ['Wnt']
    if overlapwithkegg: # only return the 6 pathways in common with KEGG.
        analyzedpathwayfile = 'inputs/pathways/netpath-dbcompare-pathways.txt'
    if allnetpath:
	analyzedpathwayfile = 'inputs/pathways/netpath-all-pathways.txt'
    else:
        analyzedpathwayfile = 'inputs/pathways/netpath-analyzed-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,1)]
    return pathways


def getAllNetPathPathways():
    """
    Reads the pathways in the NETPATHDIR directory and
    returns them
    """
    analyzedpathwayfile = 'inputs/pathways/netpath-all-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,1)]
    return pathways


def getKEGGPathways(overlapwithnetpath):
    """
    Reads the pathways for KEGG and returns them as a list.
    It also returns a dictionary of {keggid:netpathname}
    """
    if overlapwithnetpath:
        analyzedpathwayfile = 'inputs/pathways/kegg-dbcompare-pathways.txt'
    else:
        analyzedpathwayfile = 'inputs/pathways/kegg-analyzed-pathways.txt'
    pathways = [p for p in readItemSet(analyzedpathwayfile,2)]
    # dictionary of keggnames to netpath names.
    kegg2netpath = readDict(analyzedpathwayfile,2,1)
    return pathways,kegg2netpath


def getAllKEGGPathways():
    """
    Reads the pathways in the KEGGDIR directory and returns them.
    """
    pathways = glob.glob('%s/*-nodes.txt' % (KEGGDIR))
    pathways = [p.split('/')[-1] for p in pathways]
    pathways = [p.replace('-nodes.txt','') for p in pathways]
    return pathways


def getPRoutdir(alg, resultdir, netpathkeggunion):
    """
    Auxiliary function to get the Precision/Recall output
    directory, depending on whether --netpathkeggunion is
    specified.
    """
    if netpathkeggunion:
        return '%s/precision-recall/netpathkeggunion/%s/' % (resultdir,alg)
    else:
        return '%s/precision-recall/%s/' % (resultdir,alg)

def getPRsubsampleprefix(resultdir, wntsampledir, sampledir, pathway):
    """
    Auxiliary function to get the PRecision/Recall subsample
    file prefix, depending on whether this result directory
    is for --wntforexperiments (contains "wnt-all-receptors")
    """

    if 'wnt-all-receptors' in resultdir:
        return '%s/%s' % (wntsampledir,pathway)
    else:
        return '%s/%s' % (sampledir,pathway)

def postReconstructionsToGraphSpace(
        pathway, infile, infile_allreceptors, outdir, thres,gsid,
        printonly, increase=False, decrease=False, undirected=False,
        oldgs=False, posttag=False):
    """
    Posts Wnt Reconstructions to GraphSpace
    Posts two different visualizations:
        - Posts annotated networks with gene names and ranked value,
          if predictions are ranked

        - Posts un-annotated networks with labels removed

    :param infile: File of edges, which may be ranked by a 3rd column

    :param infile_allreceptors: File of edges from 'wnt-all-receptors' 
        experiments, ranked the same as infile

    :param thres: If increase==True or decrease==True, use this threshold

    :param gsid: GraphSpace ID

    :param printonly: If True, will never execute command.

    :param increase: If True, rank edges in inceasing order

    :param decrease: If True, rank edges in decreasing order

    :param undirected: If True, checks both (u,v) and (v,u) for 
        evidence sources

    :param allreceptors: If True, takes Wnt interactome from 
        wnt-all-receptors/ instead of netpath/

    :return: None

    """
    # (note: if increase==False and decrease==False, then the edges 
    # are considered an entire set and not ranked by thres)

    # PPI FILE is original interactome; this ensures that edges are 
    # directed as they were originally (not necessarily as they were 
    # after removing outgoing edges from TRs and incoming edges to 
    # receptors)

    # print annotated from infile_allreceptors
    labeledgsid = gsid+'-labeled'
    if oldgs:
        print 'ERROR: old GS functions deprecated.'
    else:
        cmd = 'python src/post-to-new-graphspace.py --infile %s --outdir %s --ppi %s --version %s --datadir %s --gsid %s --netpath %s --kegg %s --addfzd' % (infile_allreceptors,outdir,ORIGINALPPI,PPIVERSION,DATADIR,labeledgsid,pathway,pathway)
    if increase: # ranked list - pass the threshold
        cmd += ' --increase --thres %.8f' % (thres)
    if decrease: # ranked list - pass the threshold
        cmd += ' --decrease --thres %.8f' % (thres)
    if undirected: 
        cmd += ' --undirected'
    if posttag:
        cmd += ' --tag wnt-networks-all-methods'
    print cmd
    if not printonly:
        subprocess.check_call(cmd.split())

    if posttag: # only post labeled pathways.
        return

    #print unlabeled from infile
    unlabeledgsid=gsid+'-unlabeled'
    if oldgs:
        print 'ERRROR: old GS functions deprecated.'
    else:
        cmd = 'python src/post-to-new-graphspace.py --infile %s --outdir %s --ppi %s --version %s --datadir %s --gsid %s --netpath %s --kegg %s --nolabels' % (infile,outdir,ORIGINALPPI,PPIVERSION,DATADIR,unlabeledgsid,pathway,pathway)
    if increase: # ranked list - pass the threshold
        cmd += ' --increase --thres %.8f' % (thres)
    if decrease:
        cmd += ' --decrease --thres %.8f' % (thres)
    if undirected:
        cmd += ' --undirected'
    print cmd
    if not printonly:
        subprocess.check_call(cmd.split())

    return

def postNetPathReconstructionsToGraphSpace(
        pathway, infile, outdir, thres, gsid, printonly, increase=False,
        decrease=False, undirected=False, oldgs=False, posttag=False):
    """
    Post NetPath reconstructions (no extra FZD receptors, no KEGG
    coloring)
    """
    # PPI FILE is original interactome; this ensures that edges are directed as
    # they were originally  (not necessarily as they were after removing
    # outgoing edges from TRs and incoming edges to receptors)      

    # print annotated from infile_allreceptors                                                                   
    labeledgsid = gsid
    if oldgs:
        print 'ERROR: old GS functions deprecated.' 
    else:
        cmd = 'python src/post-to-new-graphspace.py --infile %s --outdir %s --ppi %s --version %s --datadir %s --gsid %s --netpath %s --nocrosstalk' % (infile,outdir,ORIGINALPPI,PPIVERSION,DATADIR,gsid,pathway)
    if increase: # ranked list - pass the threshold                                                         
        cmd += ' --increase --thres %.8f' % (thres)
    if decrease: # ranked list - pass the threshold                                                        
        cmd += ' --decrease --thres %.8f' % (thres)
    if undirected:
        cmd += ' --undirected'
    if posttag:
        cmd += ' --tag pathlinker-all-netpath-pathways'

    keggpathways,kegg2netpath = getKEGGPathways(True)
    if pathway in kegg2netpath.values():
        cmd += ' --kegg %s' % (pathway)
    print cmd
    if not printonly:
        subprocess.check_call(cmd.split())


def computePrecisionRecall(
        pathway, datadir, ppidir, edgefile, outdir, edgesortcol, negtype,
        ignorekeggpos, ignorenetpathpos, sampleoutprefix, subsamplefps,
        forceprecrec, printonly, nodefile=None, nodesortcol=None,
        descending=False, param=None, union=False):
    """
    Computes precision and recall and writes values to file.

    :param pathway: Pathway to compute precision and recall for 
        (e.g., Wnt)

    :param datadir: Directory for true edge file and true node file
        for pathway
        
    :param ppidir: Directory for pathway-specific interactomes

    :param edgefile: Predicted edges.

    :param outdir: Output directory

    :param edgesortcol: Sort column for edges. If None, then take
        edges as a set

    :param negtype: One of 'none', 'adjacent', or 'file'

    :param ignorekeggpos: ignore KEGG positives when sampling 
        negatives (negtype=='file')

    :param ignorenetpathpos: ignore NetPath positives when sampling 
        negatives (negtype='file')

    :param sampleoutprefix: prefix of subsampled negatives and 
        positives for the pathway

    :param subsamplefps: Number of negatives to sample (a factor of 
        the size of the positives

    :param forceprecrec: Continue even if files have been written

    :param printonly: If True, never execute commands

    :param nodefile: Predicted nodes (if different from edgefile)

    :param nodesortcol: Sort column for nodes

    :param descending: If True, rank in descending order

    :param param: Additional string to append to outfile

    :param union: If True, use union of NetPath and KEGG as positives

    :return: None
    """

    # Get true edge file and true node file from the data directory.
    if union:
        # TODO remove hard-coded links
        trueedgefile = 'inputs/pathways/netpath-kegg-union/%s-edges.txt' % (pathway)
        truenodefile = 'inputs/pathways/netpath-kegg-union/%s-nodes.txt' % (pathway)
    else:
        trueedgefile = '%s/%s-edges.txt' % (datadir,pathway)
        truenodefile = '%s/%s-nodes.txt' % (datadir,pathway)
        
    # make output directory.
    mkpath(outdir)
    if param == None:
        outprefix = '%s%s' % (outdir,pathway)
    else: # e.g., 'q_0.50'
        outprefix = '%s%s-%s' % (outdir,pathway,param)

    # check to see if the file exists:
    finalfile = '%s-exclude_%s-sample_%dX-node-precision-recall.txt' % (outprefix,negtype,subsamplefps)
    if not forceprecrec and os.path.isfile(finalfile):
        print finalfile
        print 'Skipping %s exclude %s: file exists. Use --forceprecrec to override.' % (pathway,negtype)
        return

    # pathway-specific interactome
    ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

    # if negtype == 'file', add kegg or netpath ignored nodes and edges
    if negtype == 'file':
        pathwaynames,kegg2netpath = getKEGGPathways(ignorenetpathpos or ignorekeggpos) 
        if ignorekeggpos:
            keggid = [k for k in kegg2netpath if kegg2netpath[k]==pathway]
            if len(keggid)==0:
                print kegg2netpath
                print pathway
                sys.exit('ERROR: Pathway %s does not exist in kegg. Exiting.' % (pathway))
            keggid = keggid[0]
            ignorededgefile = '%s/%s-edges.txt' % (KEGGDIR,keggid)
            ignorednodefile = '%s/%s-nodes.txt' % (KEGGDIR,keggid)
        elif ignorenetpathpos:
            ignorededgefile = '%s/%s-edges.txt' % (NETPATHDIR,kegg2netpath[pathway])
            ignorednodefile = '%s/%s-nodes.txt' % (NETPATHDIR,kegg2netpath[pathway])
        else:
            sys.exit('ERROR: if "file" is specified, either --ignorenetpathpositives or --ignorekeggpositives must be specified. Exiting.')

    # compute-precision-recall.py simply uses ppifile to get
    # all possible edges.
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


# ANNA POTENTIAL BUG CAUGHT SEPT 11, 2015!!  The default value for descending
# was True.  This was inconsistent with the computePrecisionRecall function
# above.  As a result, I suspect that PathLinker/ANAT/ResponseNet/APSP were
# ordered incorrectly when computing aggregate precision and recall.  I am
# confirming this now.
def computeAggregatePrecisionRecall(inputdir,negtype,ignorekeggpos,ignorenetpathpos,subsamplefps,forceprecrec,printonly,descending=False,param=None,netpath=True,union=False, allpathways=False):
    """
    Computes aggregate precision and recall

    :param inputdir: Directory of precision-recall files for individual
        pathways

    :param negtype: One of 'none' or 'adjacent'

    :param ignorekeggpos: Ignore KEGG positives when sampling negatives
        (negtype='file' and netpath = True)

    :param ignorenetpathpos: ignore NetPath positives when sampling negatives
        (negtype='file' and netpath=False)

    :param subsamplefps: Number of negatives to sample (a factor of the size
        of the positives)

    :param forceprecrec: Continue even if files have been written

    :param printonly: If True, never execute commands

    :param descending: If True, rank in descending order

    :param param: Additional string to append to outfile

    :param netpath: If True, run NetPath aggregate. Otherwise run KEGG
        aggregate.

    :param union: If True, use union of NetPath and KEGG as positives

    :return: None
    """

    if ignorekeggpos and not netpath:
        sys.exit('ERROR: cannot ignore kegg positives with kegg datasets. Exiting.')
    if ignorenetpathpos and netpath:
        sys.exit('ERROR: cannot ignore netpath positives with netpath datasets. Exiting.')        

    # check to see if the file exists:
    if param == None:
        if ignorekeggpos:
            finalfile = '%saggregate-pathways_shared_with_kegg-exclude_%s-sample_%dX-node-precision-recall.txt' % (inputdir,negtype,subsamplefps)
        elif ignorenetpathpos:
            finalfile = '%saggregate-pathways_shared_with_netpath-exclude_%s-sample_%dX-node-precision-recall.txt' % (inputdir,negtype,subsamplefps)
	elif allpathways:
            finalfile = '%saggregate-allpathways-exclude_%s-sample_%dX-node-precision-recall.txt' % (inputdir,negtype,subsamplefps)
        else:
            finalfile = '%saggregate-exclude_%s-sample_%dX-node-precision-recall.txt' % (inputdir,negtype,subsamplefps)
    else:
        if ignorekeggpos:
            finalfile = '%saggregate-pathways_shared_with_kegg-%s-exclude_%s-sample_%dX-node-precision-recall.txt' % (inputdir,param,negtype,subsamplefps)
        elif ignorenetpathpos:
             finalfile = '%saggregate-pathways_shared_with_netpath-%s-exclude_%s-sample_%dX-node-precision-recall.txt' % (inputdir,param,negtype,subsamplefps)
        elif allpathways:
            finalfile = '%saggregate-allpathways-%s-exclude_%s-sample_%dX-node-precision-recall.txt' % (inputdir,param,negtype,subsamplefps)
	else:
            finalfile = '%saggregate-%s-exclude_%s-sample_%dX-node-precision-recall.txt' % (inputdir,param,negtype,subsamplefps)
    if not forceprecrec and os.path.isfile(finalfile):
        print 'Skipping aggregate, file exists. Use --forceprecrec to override.'
        return
    #print '%s does not exist. Continuing.' % (finalfile)
    cmd = 'python src/compute-aggregate-precision-recall.py --inputdir %s --negtype %s --neg-factor %d' % (inputdir,negtype,subsamplefps)
    if netpath:
        cmd += ' --netpath'
    else:
        cmd += ' --kegg'
    if allpathways:
	cmd += ' --allpathways'
    if param != None:
        cmd += ' --param %s' % (param)
    if descending:
        cmd += ' --descending'
    if ignorekeggpos:
        cmd += ' --ignorekegg'
    if ignorenetpathpos:
        cmd += ' --ignorenetpath'
    if union:
        cmd += ' --union'
    print cmd
    if not printonly:
        subprocess.check_call(cmd.split())
    return

def performSubsampling(pathways, k, forceRecalc, forcePRRecalc, printonly, batchrun, opts):
    """
    Performs a subsampling analysis
    stuff: docs (what??)
    """

    # This analysis only makes sense to run on multiple pathways. Verify
    # that the opts obey this.
    if(opts.onlynetpathwnt or opts.wntforexperiments):
        print("Robustness study should be run on many pathways, not just wnt")
        exit(-1)

    # Validate args
    if(opts.only_generate_sample_sets and opts.batch_run):
        print("Options --only-generate-sample-sets and --batch-run should not be used together")
        exit(-1)

    # This is the directory in which all robustness related files will
    # be placed (both intermediate files and final results)
    weightedStr = "weighted" if opts.weightedppi else "unweighted"
    subsamplingDir = "outputs/%s/%s/tf-sampling/"%(opts.ppiversion, weightedStr)
    mkpath(subsamplingDir)

    # Sample sizes/increments are hardcoded constants for now
    sampleSizes = [50, 70, 90, 100, 110, 130, 150]
    nSamples = 25

    print("Using sample sizes = " + str(sampleSizes) + " nSamples = " + str(nSamples))

    print("Subsampling over pathways " + str([name for (name,_,_,_) in pathways]))

    # Generate upsampled and downsampled node sets
    if not batchrun:
        sampleNodeSets(
            pathways, subsamplingDir, forceRecalc, printonly, 
            sampleSizes, nSamples)
    
    if opts.only_generate_sample_sets:
        return

    # Run PathLinker on every sampled nodesets
    runPathLinkerSampledSets(
        pathways, k, subsamplingDir, forceRecalc, printonly, 
        batchrun, sampleSizes, nSamples)

    # Compute PR for each run, including aggregate
    if not batchrun:
        computeSampledPR(pathways, k, subsamplingDir, forceRecalc, 
        forcePRRecalc, printonly, sampleSizes, nSamples)
    
    # Plot the results
    if not batchrun:
        forcePlot = forceRecalc or forcePRRecalc
        plotRobustness(pathways, k, subsamplingDir, forcePlot, 
        printonly, sampleSizes, nSamples)

def sampleNodeSets(
        pathways, sampledSetDir, forceRecalc, printonly, sampleSizes, 
        nSamples):
    """
    Create sampled TR/Rec node sets for each pathway.
    """
    
    # Create the master list of all TFs and RECs in the network
    print("\n === Creating master receptor and TF sets ===")
    
    pathwayPathsList = " ".join([path + name + "-nodes.txt" for (name,_,path,_) in pathways])
    masterListsLoc = sampledSetDir + "master-"
    
    cmd = 'python src/create-master-TF-REC-lists.py -o %s %s'%(masterListsLoc, pathwayPathsList)
    if(forceRecalc):
        cmd += " --force-recalc"    

    print(cmd)
    if(not printonly):
        subprocess.check_call(cmd.split())
        
    
    # Generate the sampled nodeset for each pathway
    print("\n === Creating sampled nodesets ===")
    masterTFLoc = masterListsLoc + "TFs.txt"
    masterRECLoc = masterListsLoc + "RECs.txt"
    
    for (pathway,resultdir,datadir,ppidir) in pathways:
   
        # Get the pathway specific interactome
        ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)

        # Ensure the output directory exists
        pathwayDir = sampledSetDir + "/" + pathway
        if not printonly:
            mkpath(pathwayDir)

        for sampSize in sampleSizes:

            # Ensure the output directory exists
            sampleSizeDir = pathwayDir + "/" + str(sampSize)
            if not printonly:
                mkpath(sampleSizeDir)
            sampleDir = sampleSizeDir + "/sampled-nodesets/"
            if not printonly:
                mkpath(sampleDir)

            # Perform the sampling
            nodeFile = datadir + pathway + "-nodes.txt"
            cmd = 'python src/create-sample-sets.py %s -o %s -n %s -p %s --ppi %s --all-TFs %s --all-Recs %s'%(nodeFile, sampleDir, nSamples, sampSize, ppifile, masterTFLoc, masterRECLoc)

            if(forceRecalc):
                cmd += " --force-recalc"    

            print("\n" + cmd)
            if(not printonly):
                subprocess.check_call(cmd.split())
            
        
    print 'Done creating sampled nodesets\n'


def runPathLinkerSampledSets(
    pathways, k, sampledSetDir, forceRecalc, printonly, batchrun, 
    sampleSizes, nSamples):
    """
    Run PathLinker on each of the sampled nodesets
    """
        
    print("\n === Running PathLinker for all sampled sets ===")

    # Run PathLinker for each pathway x percent x sample
    for (pathway,resultdir,datadir,ppidir) in pathways:
   
        # Get the pathway specific interactome
        ppifile = '%s/%s-interactome.txt' % (ppidir,pathway)
        pathwayDir = sampledSetDir + "/" + pathway

        for sampSize in sampleSizes:

            sampleSizeDir = pathwayDir + "/" + str(sampSize)
            sampleDir = sampleSizeDir + "/sampled-nodesets/"

            outDir = sampleSizeDir + "/pathlinker-results/"
            if not printonly:
                mkpath(outDir)

            for iSamp in range(nSamples):

                sampledNodeFile = sampleDir + "sample_%s-%s-nodes.txt"%(sampSize, iSamp)
                outprefix = outDir + "sample-%s_%s-%s-"%(pathway, sampSize, iSamp)

                # pathlinker command
                testFile = '%sk_%d-ranked-edges.txt' % (outprefix,k)
                if forceRecalc or not os.path.isfile(testFile):

                    # Mark the file, to stop another instance of this program
                    # from overwriting it when run simultaneously. (Note: this
                    # is not a good solution to this problem, but it works well
                    # enough here)
                    if batchrun and not printonly:
                        testF = open(testFile, 'w')
                        testF.write("Working...\n")
                        testF.close()

                    script = '/home/annaritz/src/python/PathLinker/PathLinker-1.0/PathLinker.py'
                    cmd = 'python %s -k %d --write-paths --output %s %s %s' % (script,k,outprefix,ppifile,sampledNodeFile) 
                    print("\n" + cmd)
                    if not printonly:
                        subprocess.check_call(cmd.split())
                else:
                    print 'Skipping %s: %s exists. Use --forcealg to override.' % (pathway,'%sk_%d-paths.txt' % (outprefix,k))


def computeSampledPR(
        pathways, k, sampledSetDir, forceRecalc, forcePRRecalc, printonly, 
        sampleSizes, nSamples):
    """
    Compute interpolated precision-recall for each of the sampled runs
    """
        
    print("\n === Computing interpolated PR for all sampled sets ===")


    # TODO The way filenames/paths are managed here is somewhat sloppy,
    # could use a rewrite. It works like this because the PR file needs
    # to be able to find the ranks for each pathway to compute aggregate
    # PR.

    # TODO Right now the paths to to the subsampled negative node/edge
    # files are basically hardcoded. Need to work with Anna to make sure
    # it uses the proper ones.
    negTypes = ['none','adjacent']
    
    # Whether or not to recompute the PR if files already exist
    recomputePR = forceRecalc or forcePRRecalc
    
    for negType in negTypes:

        for percent in sampleSizes:

            for iTrial in range(nSamples):

                # This character is used to delimit paths when they are
                # passed as a list to do the aggregate calculation. It must
                # not appear in any path.
                delimChar = ','
                 

                # These lists acccumulate the lists of paths to be used
                # for the aggregate calculation
                aggEdgeRankLocs = []
                aggNegNodeLocs = []
                aggNegEdgeLocs = []
                aggNodeTypeLocs = []
                aggEdgeTypeLocs = []
                pathwayNames = []

                # Calculate PR for each percent x sample x pathway
                for (pathway,resultdir,datadir,ppidir) in pathways:
               
                    # The paths for the inputs in this PR calculation
                    edgeRankLoc = "%s%s/%d/pathlinker-results/sample-%s_%d-%d-k_%d-ranked-edges.txt"%(sampledSetDir,pathway,percent,pathway,percent,iTrial,k)
                    negNodeLoc = 'outputs/pathlinker-signaling-children-reg/weighted/samples-exclude-%s/%s-exclude_%s-50X-negative-nodes.txt'%(negType,pathway,negType)
                    negEdgeLoc = 'outputs/pathlinker-signaling-children-reg/weighted/samples-exclude-%s/%s-exclude_%s-50X-negative-edges.txt'%(negType,pathway,negType)
                    nodeTypeLoc = datadir + pathway + "-nodes.txt"
                    edgeTypeLoc = datadir + pathway + "-edges.txt"

                    # Add the paths to the aggregate lists
                    aggEdgeRankLocs.append(edgeRankLoc)
                    aggNegNodeLocs.append(negNodeLoc)
                    aggNegEdgeLocs.append(negEdgeLoc)
                    aggNodeTypeLocs.append(nodeTypeLoc)
                    aggEdgeTypeLocs.append(edgeTypeLoc)
                    pathwayNames.append(pathway)

                    # The destination for the output file
                    outLocDir = "%s/%s/%d/PR/"%(sampledSetDir, pathway, percent)
                    mkpath(outLocDir)
                    outLoc = "%s/exclude-%s_trial_%d"%(outLocDir, negType, iTrial)

                    # Assemble the command to the subprogram
                    cmd = "python src/calc-interpolated-PR.py %s --outname %s --node-type-file %s --edge-type-file %s --neg-edge-file %s --neg-node-file %s --pathway-name %s "%(edgeRankLoc, outLoc, nodeTypeLoc, edgeTypeLoc, negEdgeLoc, negNodeLoc, pathway) 
                    print("\n" + cmd)
  
                    # Run the subprogram
                    testFile = outLoc + "-interpedEdgePR.txt"
                    if recomputePR or not os.path.isfile(testFile):
                        if not printonly:
                            subprocess.check_call(cmd.split())
                    else:
                        print("Skipping computation of %s (and the corresponding node file)-- file exists."%(testFile))

                
                
                # Run the aggregate calculation for this percent x sample

                # The destination for the output file
                # TODO This directory structure isn't great
                outLocDir = "%s/aggregate-PR-exclude-%s"%(sampledSetDir, negType)
                mkpath(outLocDir)
                outLoc = "%s/percent_%d-trial_%d"%(outLocDir,percent,iTrial)


                cmd = "python src/calc-interpolated-PR.py --agg %s %s --outname %s --node-type-file %s --edge-type-file %s --neg-edge-file %s --neg-node-file %s  --pathway-name %s"%(delimChar, delimChar.join(aggEdgeRankLocs), outLoc, delimChar.join(aggNodeTypeLocs), delimChar.join(aggEdgeTypeLocs), delimChar.join(aggNegEdgeLocs), delimChar.join(aggNegNodeLocs), delimChar.join(pathwayNames)) 
                print("\n" + cmd)
  
                # Run the subprogram
                testFile = outLoc + "-interpedEdgePR.txt"
                print("Testing file " + testFile)
                if recomputePR or not os.path.isfile(testFile):
                    if not printonly:
                        subprocess.check_call(cmd.split())
                else:
                    print("Skipping computation of %s (and the corresponding node file)-- file exists."%(testFile))



def plotRobustness(pathways, k, sampledSetDir, forcePlot, printonly, sampleSizes, nSamples):
    """
    Visualize the results of the robustness study
    """
        
    # For now this only makes plots of the aggregate results, as
    # these are generally the ones we're interested in.

    print("\n === Plotting aggregate robustness results ===")

    # TODO Right now the paths to to the subsampled negative node/edge
    # files are basically hardcoded. Need to work with Anna to make sure
    # it uses the proper ones.
    negTypes = ['none','adjacent']
    
    outDir = sampledSetDir + "viz/"
    mkpath(outDir)

    for negType in negTypes:

        nodeFilePattern = sampledSetDir + "aggregate-PR-exclude-" + str(negType) + "/percent_%d-trial_%d-interpedNodePR.txt"
        edgeFilePattern = sampledSetDir + "aggregate-PR-exclude-" + str(negType) + "/percent_%d-trial_%d-interpedEdgePR.txt"
        outFile = outDir + "aggregate-exclude_%s-k_%d-viz-"%(negType,k)

        cmd = "python src/plot-robustness.py --outname %s --node-PR-file %s --edge-PR-file %s --n-samples %d"%(outFile, nodeFilePattern, edgeFilePattern, nSamples)
        print("\n" + cmd)
        
        testFile = outDir + "Edge-Medians.png"
        print("Testing file " + testFile)
        if forcePlot or not os.path.isfile(testFile):
            if not printonly:
                subprocess.check_call(cmd.split())
        else:
            print("Skipping plotting of exclude-%s because file exists: %s"%(negType, testFile))

# Okay, so it looks like some sort wrapper to run any algorithm on multiple
# pathways would be good. What would be BEST would be able to say:
# run_multiple_pathways(pathways, algorithms)
#     for pathway in pathways:
#         for algorithm in algorithms: 
#             run algorithm function 
#
###############################################################################
class RankingAlgorithm(object):
    def run(self, file_location_context):
        raise NotImplementedError()


    def run_if_forced(self, file_location_context, should_force):
        if (should_force):
            self.run(file_location_context)
        else:
            if not self.output_previously_written(file_location_context): 
                self.run(file_location_context)


    def get_name(self):
        raise NotImplementedError()


    def output_previously_written(self, file_location_context):
        raise NotImplementedError()


class FileLocationContext(object):
    def __init__(self, pathway_name, pathway_dir, ppi_file_dir, output_dir):
        self.pathway_name = pathway_name
        self.pathway_dir = pathway_dir
        self.ppi_file_dir = ppi_file_dir
        self.output_dir = output_dir


    def get_ppi_file_location(self):
        return '%s/%s-interactome.txt' % (self.ppi_file_dir, self.pathway_name)


    def get_sources_and_targets(self):
        return '%s/%s-nodes.txt' %(self.pathway_dir, self.pathway_name)


class OriginalPathLinker(RankingAlgorithm):
    # TODO: The intent here is to mimic the functionality of the original 
    # PathLinker algorithm, which had a heavy weight penalty that amounted to
    # dividing by the sum of all of the weights in the network. Not sure how
    # I want to go about handling this. Each pathway specific network would
    # technically have a different weight, so hardcoding it might be rough.
    # I could INSTEAD add it to PathLinker.
    None

# TODO: To mimic the old functionality of printing instead of actually running
# as well as by default not re-running if the algorithm has already run,
# we can define more methods here
class PathLinker(RankingAlgorithm):
    def __init__(self, k):
        self.k = k


    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        subprocess.call([
            "python",
            "src/external/pathlinker/PathLinker.py",
            "-k",
            str(self.k),
            "--write-paths",
            "--output",
            self.get_output_prefix_from_context(file_location_context),
            file_location_context.get_ppi_file_location(),
            file_location_context.get_sources_and_targets()
            ])


    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_output_prefix_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s-" % file_location_context.pathway_name
        return out_prefix


    def get_output_file_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s-" % file_location_context.pathway_name \
            + "k_%d-paths.txt" % self.k 
        return out_prefix


    def get_name(self):
        return "pathlinker"


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_prefix_from_context(file_location_context) \
            + "k_%d-paths.txt" % self.k 
        if os.path.exists(outfile):
            return True
        else:
            return False
        
# TODO: For this and other variations of PathLinker, if the only difference is
# a single option passed in, maybe its constructor should just take more
# parameters
class PageRankPathLinker(RankingAlgorithm):
    def __init__(self, q, k):
        self.q = q
        self.k = k


    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        subprocess.call([
            "python",
            "src/external/pathlinker/PathLinker.py",
            "--PageRank",
            "-q",
            str(self.q),
            "-k",
            str(self.k),
            "--write-paths",
            "--output",
            self.get_output_prefix_from_context(file_location_context),
            file_location_context.get_ppi_file_location(),
            file_location_context.get_sources_and_targets()
            ])


    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_output_prefix_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s-q_%.2f" % file_location_context.pathway_name, q

        return out_prefix


    def get_name(self):
        return "pagerank-pathlinker"


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_prefix_from_context(file_location_context) \
            + "k_%d-paths.txt" % self.k 
        if os.path.exists(outfile):
            return True
        else:
            return False
        

class CycLinker(RankingAlgorithm):
    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        cmd = 'java %s %s %s %s' % \
            ('CycLinker', 
             file_location_context.get_ppi_file_location(),
             file_location_context.get_sources_and_targets(),
             os.path.abspath(
                self.get_output_prefix_from_context(file_location_context)))
        print cmd
        curr_dir = os.getcwd()
        # TODO: Add https://github.com/jlaw9/Cyclinker as a submodule instead 
        os.chdir('/home/jeffl/git-workspace/CycLinker/src')
        subprocess.check_call(cmd.split())
        os.chdir(curr_dir)

    
    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_output_prefix_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name

        return out_prefix


    def get_name(self):
        return "cyclinker"


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_prefix_from_context(file_location_context) \
            + "-ranked-edges-.txt"
        if os.path.exists(outfile):
            return True
        else:
            return False


class PageRankCycLinker(RankingAlgorithm):
    def __init__(self, q):
        self.q = q


    def run(self, file_location_context):
        pagerank = PageRank(self.q)
        self.ensure_output_dir_exists(file_location_context)
        self.ensure_pagerank_output_exists(file_location_context)

        cmd = 'java %s %s %s %s' % \
            ('CycLinker', 
             pagerank.get_edge_flux_output_name_from_context(
                file_location_context),
             file_location_context.get_sources_and_targets(),
             os.path.abspath(
                self.get_output_prefix_from_context(file_location_context)))
        curr_dir = os.getcwd()
        # TODO: Add https://github.com/jlaw9/Cyclinker as a submodule instead 
        os.chdir('/home/jeffl/git-workspace/CycLinker/src')
        subprocess.check_call(cmd.split())
        os.chdir(curr_dir)


    def ensure_pagerank_output_exists(self, file_location_context):
        pagerank = PageRank(self.q)

        if not pagerank.output_previously_written(file_location_context):
            raise Exception("Error: This algorithm must be run alongside "
                "the --pathlinker option, as it utilizes PathLinker output "
                "as input. At the very least, the --pathlinker option must "
                "be used in a previous run to generate the output file for " 
                "the provided -k value (which has a default if no value is "
                "provided.")


    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_output_prefix_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name \
            + "-q_%.2f" % (self.q)

        return out_prefix


    def get_name(self):
        return "pagerank-cyclinker"


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_prefix_from_context(file_location_context) \
            + "-ranked-edges-.txt"
        if os.path.exists(outfile):
            return True
        else:
            return False
    

class ShortestPaths(RankingAlgorithm):
    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        script = '/home/annaritz/src/python/CellCycle/shortest_paths.py'
        cmd = ('python %s --network %s --annotations %s --out %s --include-ties'
               ' --weight --log-transform') % \
            (script,
             file_location_context.get_ppi_file_location(),
             file_location_context.get_sources_and_targets(),
             self.get_output_file_from_context(file_location_context))

        subprocess.check_call(cmd.split())   
    
    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_output_file_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name \
            + "-shortest-paths.txt" 

        return out_prefix


    def get_name(self):
        return "shortestpaths"


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_file_from_context(file_location_context)
        if os.path.exists(outfile):
            return True
        else:
            return False


class BowTieBuilder(RankingAlgorithm):
    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        script = '/home/annaritz/src/python/CellCycle/bowtiebuilder.py'         
        cmd = ('python %s --network %s --annotations %s --out %s '
               '--weight --log-transform') % \
            (script,
             file_location_context.get_ppi_file_location(),
             file_location_context.get_sources_and_targets(),
             self.get_output_file_from_context(file_location_context))

        subprocess.check_call(cmd.split())   
   

    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_output_file_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name \
            + "-bowtiebuilder.txt" 

        return out_prefix


    def get_name(self):
        return "bowtiebuilder"


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_file_from_context(file_location_context)
        if os.path.exists(outfile):
            return True
        else:
            return False


class InducedSubgraph(RankingAlgorithm):
    def __init__(self, k):
        self.k = k

    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        self.ensure_pathlinker_output_exists(file_location_context)
        script = '/data/annaritz/signaling/2014-06-linker/src/order-by-induced-subgraph.py'
        cmd = 'python %s --pathsfile %s --outfile %s --ppi %s' % \
            (script,
             self.get_corresponding_pathlinker_outfile(file_location_context),
             self.get_output_file_from_context(file_location_context),
             file_location_context.get_ppi_file_location())
        subprocess.check_call(cmd.split())   


    def ensure_pathlinker_output_exists(self, file_location_context):
        pathlinker_out = self.get_corresponding_pathlinker_outfile(
            self, file_location_context)

        if not os.path.exists(pathlinker_out):
            raise Exception("Error: This algorithm must be run alongside "
                "the --pathlinker option, as it utilizes PathLinker output "
                "as input. At the very least, the --pathlinker option must "
                "be used in a previous run to generate the output file for " 
                "the provided -k value (which has a default if no value is "
                "provided.")
    

    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_corresponding_pathlinker_outfile(self, file_location_context):
        pathlinker = PathLinker(self.k)
        return pathlinker.get_output_file_from_context(file_location_context)


    def get_output_file_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name \
            + "-induced-subgraph.txt" 

        return out_prefix


    def get_name(self):
        return "inducedsubgraph"


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_file_from_context(file_location_context)
        if os.path.exists(outfile):
            return True
        else:
            return False


class RerankPathLinker(RankingAlgorithm):
    def __init__(self, k):
        self.k = k

    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        self.ensure_pathlinker_output_exists(file_location_context)
        script = '/data/annaritz/signaling/2014-06-linker/src/recount-ksp.py'   
        cmd = 'python %s --pathsfile %s --outputprefix %s' % \
            (script,
             self.get_corresponding_pathlinker_outfile(file_location_context),
             self.get_output_prefix_from_context(file_location_context))

        subprocess.check_call(cmd.split())


    def ensure_pathlinker_output_exists(self, file_location_context):
        pathlinker_out = self.get_corresponding_pathlinker_outfile(
            self, file_location_context)

        if not os.path.exists(pathlinker_out):
            raise Exception("Error: This algorithm must be run alongside "
                "the --pathlinker option, as it utilizes PathLinker output "
                "as input. At the very least, the --pathlinker option must "
                "be used in a previous run to generate the output file for " 
                "the provided -k value (which has a default if no value is "
                "provided.")
    

    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_corresponding_pathlinker_outfile(self, file_location_context):
        pathlinker = PathLinker(self.k)
        return pathlinker.get_output_file_from_context(file_location_context)


    def get_output_prefix_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s-reranked-pathlinker" % file_location_context.pathway_name

        return out_prefix


    def get_name(self):
        return "rerankedpathlinker"


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_prefix_from_context(file_location_context) \
            + "-unique-edges_paths.txt"
        if os.path.exists(outfile):
            return True
        else:
            return False


class PageRank(RankingAlgorithm):
    def __init__(self, q):
        self.q = q

    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        script = '/home/annaritz/src/python/PathLinker/PathLinker-1.0/PathLinker.py'
        cmd = 'python %s --PageRank -q %s -k %d --output %s %s %s' % \
            (script,
             self.q,
             1,
             self.get_output_prefix_from_context(file_location_context),
             file_location_context.get_ppi_file_location(),
             file_location_context.get_sources_and_targets())

        subprocess.check_call(cmd.split())


    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    def get_output_prefix_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name \
            + "-q_%.2f" % self.q 

        return out_prefix


    def get_edge_flux_output_name_from_context(self, file_location_context):
        return self.get_output_prefix_from_context(file_location_context) + \
            "-edge-fluxes.txt"


    def get_name(self):
        return "pagerank"


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_prefix_from_context(file_location_context) \
            + "-node-pagerank.txt"
        if os.path.exists(outfile):
            return True
        else:
            return False


class EQED(RankingAlgorithm):
    def __init__(self, input_current):
        self.input_current = input_current 


    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        script = '/data/annaritz/sig-path-other-methods/src/eQED.py'            
        cmd = 'python %s -e %s -n %s -i %d -o %s -l' % \
            (script,
             file_location_context.get_ppi_file_location(),
             file_location_context.get_sources_and_targets(),
             self.input_current,
             self.get_output_prefix_from_context(file_location_context))

        subprocess.check_call(cmd.split())


    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_output_prefix_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name

        return out_prefix


    def get_name(self):
        return "eqed"


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_prefix_from_context(file_location_context) \
            + "-eqed-edges.out"
        if os.path.exists(outfile):
            return True
        else:
            return False


class ResponseNet(RankingAlgorithm):
    def __init__(self, gamma):
        self.gamma = gamma


    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        script = '/data/annaritz/sig-path-other-methods/src/ResponseNet.py'
        cmd = 'python %s -e %s -n %s -o %s -g %d' % \
            (script,
             file_location_context.get_ppi_file_location(),
             file_location_context.get_sources_and_targets(),
             self.get_output_prefix_from_context(file_location_context),
             self.gamma)

        subprocess.check_call(cmd.split()) 


    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_name(self):
        return "responsenet"


    def get_output_prefix_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name \
            + "-gamma_%d" % self.gamma   


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_prefix_from_context(file_location_context) \
            + "_responsenet-edges.out"
        if os.path.exists(outfile):
            return True
        else:
            return False
        

class PCSF(RankingAlgorithm):
    def __init__(self, prize, omega):
        self.prize = prize
        self.omega = omega

    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)
        script = '/home/jeffl/svnrepo/src/python/Algorithms/PCSF_weighted.py'   
        cmd = 'python %s -e %s -n %s -o %s -p %d --omega %.2f' % \
            (script,
             file_location_context.get_ppi_file_location(),
             file_location_context.get_sources_and_targets(),
             self.get_output_prefix_from_context(file_location_context),
             self.prize,
             self.omega)

        subprocess.check_call(cmd.split())


    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_name(self):
        return "pcsf"


    def get_output_prefix_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name \
            + "-prize%d-omega%.2f" % (self.prize, self.omega)


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_prefix_from_context(file_location_context) \
            + "_PCSF-edges.out"
        if os.path.exists(outfile):
            return True
        else:
            return False


class ANAT(RankingAlgorithm):
    def __init__(self, alpha):
        self.alpha = alpha


    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)

        script = '/data/annaritz/signaling/2014-06-linker/src/run-anat-weighted.py'
        cmd = 'python %s -n %s -a %.2f -o %s --ppi %s' % \
            (script,
             file_location_context.get_sources_and_targets(),
             self.alpha,
             self.get_output_prefix_from_context(file_location_context),
             file_location_context.get_ppi_file_location())
        subprocess.check_call(cmd.split())


    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_name(self):
        return "anat"


    def get_output_prefix_from_context(self, file_location_context):
        out_prefix = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name \
            + "-alpha%.2f" % self.alpha


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_prefix_from_context(file_location_context) \
            + "-edges.out"
        if os.path.exists(outfile):
            return True
        else:
            return False


class IPA(RankingAlgorithm):
    def __init__(self, nmax):
        self.nmax = nmax


    def run(self, file_location_context):
        self.ensure_output_dir_exists(file_location_context)

        outfile = '%s/%s-nmax%d.out' % (outdir,pathway,nmax)                        
        script = '/home/jeffl/svnrepo/src/python/Algorithms/ipa-network-generation.py'
        cmd = 'python %s --ppi %s --nodes %s --nmax %d --outfile %s' % \
            (script,
             file_location_context.get_ppi_file_location(),
             file_location_context.get_sources_and_targets(),
             self.nmax,
             outfile)                          
        subprocess.check_call(cmd.split())             


    def ensure_output_dir_exists(self, file_location_context):
        out_dir = file_location_context.output_dir \
            + "/%s/" % self.get_name()

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


    def get_name(self):
        return "ipa"


    def get_output_file_from_context(self, file_location_context):
        outfile = file_location_context.output_dir \
            + "/%s/" % self.get_name() \
            + "/%s" % file_location_context.pathway_name \
            + "-nmax%d.out" % self.nmax


    def output_previously_written(self, file_location_context):
        outfile = self.get_output_file_from_context(file_location_context)
        if os.path.exists(outfile):
            return True
        else:
            return False


def get_algorithm_object_list_from_opts(opts):
    algs = []
    if opts.pathlinker:
        algs.append(PathLinker(opts.k))

    if opts.pagerank_pathlinker:
        if not opts.varyparams:
            algs.append(PageRankPathLinker(opts.q, opts.k))
        else:
            for q in VARYPARAMS['q']:
                algs.append(PageRankPathLinker(q, opts.k))
            
    if opts.pagerank_cyclinker:
        if not opts.varyparams:
            algs.append(PageRankCycLinker(opts.q))
        else:
            for q in VARYPARAMS['q']:
                algs.append(PageRankCycLinker(q))

    if opts.cyclinker:
        algs.append(CycLinker())

    if opts.shortestpaths:
        algs.append(ShortestPaths())

    if opts.bowtiebuilder:
        algs.append(BowTieBuilder)

    if opts.inducedsubgraph:
        algs.append(InducedSubgraph(opts.k))

    if opts.rerank:
        algs.append(RerankPathLinker(opts.k))

    if opts.pagerank:
        if not opts.varyparams:
            algs.append(PageRank(opts.q))
        else:
            for q in VARYPARAMS['q']:
                algs.append(PageRank(q))

    if opts.eqed:
        algs.append(EQED(opts.inputcurrent))

    if opts.responsenet:
        if not opts.varyparams:
            algs.append(ResponseNet(opts.gamma))
        else:
            for gamma in VARYPARAMS['gamma']:
                algs.append(ResponseNet(gamma))

    if opts.pcsf:
        if not opts.varyparams:
            algs.append(PCSF(opts.prize, opts.omega))
        else:
            for prize in VARYPARAMS['prize']:
                for omega in VARYPARAMS['omega']:
                    algs.append(PCSF(prize, omega))
                    
    if opts.anat:
        if not opts.varyparams:
            algs.append(ANAT(opts.alpha))
        else:
            for alpha in VARYPARAMS['alpha']:
                algs.append(ANAT(alpha))

    if opts.ipa:
        if not opts.varyparams:
            algs.append(IPA(opts.nmax))
        else:
            for nmax in VARYPARAMS['nmax']:
                algs.append(IPA(nmax))

    return algs


if __name__=='__main__':
    main(sys.argv)
