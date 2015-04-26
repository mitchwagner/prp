# Takes a golden edge set and an edge ranking and produces interpolated
# precrec files

import sys
import random as rand

import networkx as nx

from optparse import OptionParser, OptionGroup

import numpy as np

#from digraph import DiGraph


def main(args):

    # Set up the parser and options

    usage ='''calc-interpolated-PR.py EDGE_RANK_FILE [options]\n\n
         EDGE_RANK_FILE should a ranked edge file from pathlinker.
         '''

    parser = OptionParser(usage=usage)

    parser.add_option('-o','--outname',default='',\
                                      help='Prepend to results. Default: ""')

    parser.add_option('','--node-type-file',default='',\
                                      help='The pathway node file. Used to find the positive nodes. default: ""')
    parser.add_option('','--edge-type-file',default='',\
                                      help='The pathway edge file. Used to find the positive edges. default: ""')

    parser.add_option('','--neg-edge-file',default='',\
                                      help='positive/negative edge set file. default: ""')
    parser.add_option('','--neg-node-file',default='',\
                                      help='Positive/negative node set file. Default: ""')
    parser.add_option('','--pathway-name',default='',\
                                      help='The name of the pathway. Used only for debugging.. Default: ""')
    parser.add_option('','--agg',default='',\
                                      help='If true, compute an aggregate precision recall over many pathways. If this option is used, the edge-rank-file, pos-neg-edge-file, and pos-neg-node-file should all be lists of file paths, rather than just paths. The argument of this option is the character to split this lists on.')

    (opts, args) = parser.parse_args()

    EDGE_RANK_FILE = args[0]


    # The master lists
    negEdges = set()
    posEdges = set()

    negNodes = set()
    posNodes = set()

    # Lists of the files to process as inputs. These lists will only
    # have > 1 element if the --agg option is specified
    rankedEdgeFiles = []
    nodeTypeFiles = []
    edgeTypeFiles = []
    negNodeFiles = []
    negEdgeFiles = []
    pathwayNames = []

    # Process a single pathway

    # Read in the pathways
    if opts.agg:
        splitChar = opts.agg
        rankedEdgeFiles = EDGE_RANK_FILE.split(splitChar)
        nodeTypeFiles = opts.node_type_file.split(splitChar)
        edgeTypeFiles = opts.edge_type_file.split(splitChar)
        negNodeFiles = opts.neg_node_file.split(splitChar)
        negEdgeFiles = opts.neg_edge_file.split(splitChar)
        pathwayNames = opts.pathway_name.split(splitChar)
        print("Running aggregate PR on %d pathways: %s"%(len(rankedEdgeFiles), str(pathwayNames)))
    else:
        rankedEdgeFiles.append(EDGE_RANK_FILE)
        nodeTypeFiles.append(opts.node_type_file)
        edgeTypeFiles.append(opts.edge_type_file)
        negNodeFiles.append(opts.neg_node_file)
        negEdgeFiles.append(opts.neg_edge_file)
        pathwayNames.append(opts.pathway_name)
        print("Running single PR on one pathway: " + str(pathwayNames))

    nPathways = len(pathwayNames)

    print("Reading nodes and edges from files")

    ## Read in positive and negative nodes/edges
    for iPath in range(nPathways):

        pathway = pathwayNames[iPath]
        print("\nProcessing pathway " + pathway)

        # Read positive nodes
        posNodeFile = nodeTypeFiles[iPath]
        print("Positive node file = " + posNodeFile)
        for line in open(posNodeFile, 'r').readlines():
            if line[0] == '#':
                continue
            items = line.strip().split('\t')

            # All nodes in this file are positive
            posNodes.add((items[0], pathway))

        # Read negative nodes
        negNodeFile = negNodeFiles[iPath]
        print("Negative node file = " + negNodeFile)
        for line in open(negNodeFile, 'r').readlines():
            if line[0] == '#':
                continue
            items = line.strip().split('\t')

            # All nodes in this file are negative
            negNodes.add((items[0], pathway))


        # Read positive edges
        posEdgeFile = edgeTypeFiles[iPath]
        print("Positive edge file = " + posEdgeFile)
        for line in open(posEdgeFile, 'r').readlines():
            if(line[0] == '#'):
                continue

            items = line.strip().split('\t')
            sortedItems = sorted((items[0],items[1]))
            sortedEdge = (sortedItems[0], sortedItems[1], pathway)

            # All edges in this file are positive
            posEdges.add(sortedEdge)

        # Read negative edges
        negEdgeFile = negEdgeFiles[iPath]
        print("Negative edge file = " + negEdgeFile)
        for line in open(negEdgeFile, 'r').readlines():
            if(line[0] == '#'):
                continue

            items = line.strip().split('\t')
            sortedItems = sorted((items[0],items[1]))
            sortedEdge = (sortedItems[0], sortedItems[1], pathway)

            # All edges in this file are negative
            negEdges.add(sortedEdge)


    # Validate nodes
    if len(posNodes) == 0:
        print("Positive node set is empty")
        exit(-1)
    if len(negNodes) == 0:
        print("Negative node set is empty")
        exit(-1)
    if len(posNodes.intersection(negNodes)) > 0:
        print("Invalid edge sets. Nonzero intersection: " + posNodes.intersection(negNodes))
        negNodes = negNodes - posNodes # Fix the bad sets until this is resolved
        #exit(-1) FIXME

    # Validate  edges
    if len(posEdges) == 0:
        print("Positive edge set is empty")
        exit(-1)
    if len(negEdges) == 0:
        print("Negative edge set is empty")
        exit(-1)
    if len(posEdges.intersection(negEdges)) > 0:
        print("Invalid edge sets. Nonzero intersection: " + str(posEdges.intersection(negEdges)))
        negEdges = negEdges - posEdges # Fix the bad sets until this is resolved
        #exit(-1)

    print("Test set has %d positive nodes and %d negative nodes"%(len(posNodes), len(negNodes)))
    print("Test set has %d positive edges and %d negative edges"%(len(posEdges), len(negEdges)))

    ## Process through the rankfile and find the node and edge precrec curve
    print("Computing precrec")

    nodePts = []
    edgePts = []
    posNodeCount = 0
    negNodeCount = 0
    posEdgeCount = 0
    negEdgeCount = 0
    seenNodes = set()

    rankedEdges = []

    ## Read in positive and negative nodes/edges
    for iPath in range(nPathways):

        pathway = pathwayNames[iPath]
        print("Computing PR for pathway = " + pathway)

        edgeRankFile = rankedEdgeFiles[iPath]
        print("Processing " + edgeRankFile)
        for line in open(edgeRankFile, 'r').readlines():
            if(line[0] == '#'):
                continue

            items = line.strip().split('\t')
            rankVal = float(items[2])
            sortedItems = sorted((items[0],items[1]))
            sortedEdge = (rankVal, sortedItems[0], sortedItems[1], pathway)
            rankedEdges.append(sortedEdge)

    rankedEdges.sort()

    for edge in rankedEdges:

        sortedEdge = (edge[1],edge[2],edge[3])

        # Edge precrec
        if sortedEdge in posEdges:
            posEdgeCount += 1
        if sortedEdge in negEdges:
            negEdgeCount += 1
        if sortedEdge in posEdges or sortedEdge in negEdges:
            rec = float(posEdgeCount) / len(posEdges)
            prec = float(posEdgeCount) / (posEdgeCount + negEdgeCount)
            edgePts.append((rec, prec))

        # Node precrec
        newNode = False
        n1 = (sortedEdge[0], sortedEdge[2])
        if n1 not in seenNodes:
            seenNodes.add(n1)
            if n1 in posNodes:
                posNodeCount += 1
                newNode = True
            if n1 in negNodes:
                negNodeCount += 1
                newNode = True

        n2 = (sortedEdge[1], sortedEdge[2])
        if n2 not in seenNodes:
            seenNodes.add(n2)
            if n2 in posNodes:
                posNodeCount += 1
                newNode = True
            if n2 in negNodes:
                negNodeCount += 1
                newNode = True

        if newNode:
            rec = float(posNodeCount) / len(posNodes)
            prec = float(posNodeCount) / (posNodeCount + negNodeCount)
            nodePts.append((rec, prec))

    # Make sure neither list is actually empty, so subsequent scripts
    # don't error out. (Happens if there were no positives)
    if(len(edgePts) == 0):
        edgePts.append((0,0))
    if(len(nodePts) == 0):
        nodePts.append((0,0))

    edgePts = np.array(edgePts)
    nodePts = np.array(nodePts)

    ## Interpolate precrec over nPts and write to file
    Npts = 1000
    recPoints = np.linspace(0, 1, Npts)

    # Edges
    edgeOutFileName = opts.outname + "-interpedEdgePR.txt"
    print("Writing edgeOutFile to " + edgeOutFileName)
    edgeOutFile = open(edgeOutFileName, 'w')

    interpVals = np.interp(recPoints, edgePts[:,0], edgePts[:,1], right=-1)
    for i in range(Npts):

        # Only write the curve for the points the that showed up in the
        # search
        if(interpVals[i] == -1):
            break

        edgeOutFile.write("%f\t%f\n"%(recPoints[i],interpVals[i]))

    edgeOutFile.close()

    # Nodes
    nodeOutFileName = opts.outname + "-interpedNodePR.txt"
    print("Writing nodeOutFile to " + nodeOutFileName)
    nodeOutFile = open(nodeOutFileName, 'w')

    interpVals = np.interp(recPoints, nodePts[:,0], nodePts[:,1], right=-1)
    for i in range(Npts):

        # Only write the curve for the points the that showed up in the
        # search
        if(interpVals[i] == -1):
            break

        nodeOutFile.write("%f\t%f\n"%(recPoints[i],interpVals[i]))

    nodeOutFile.close()



















if __name__ == '__main__':
    main(sys.argv)
