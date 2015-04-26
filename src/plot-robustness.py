# Generate plots and output statistics for the robustness (subsampling)
# study

# Nicholas Sharp - nsharp3@vt.edu

import os
import sys
import subprocess
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt

def main(args):

    usage = 'python plot-robustness.py  [options]\n'
    parser = OptionParser(usage=usage)

    parser.add_option('-o','--outname',default='sampled-set',\
                          help='Prepend to results. Default: "sampled-set"')

    parser.add_option('','--node-PR-file',default='',\
                          help='The file containing the interpolated node PR.')

    parser.add_option('','--edge-PR-file',default='',\
                          help='The file containing the interpolated edge PR.')

    parser.add_option('-n','--n-samples',type='int',metavar='INT',\
                          help='# number of samples')

    parser.add_option('','--plot-title',default='%s Recall (exclude-none) at $\\rho = $ %d',\
                          help='Prepend to plot titles.')
    parser.add_option('','--plot-median-title',default='Median %s Recall (exclude-none)',\
                          help='Prepend to plot titles')


    parser.add_option('', '--show-plots', action='store_true', default=False, help='Show plots. Default:False')


    # parse the command line arguments
    (opts, args) = parser.parse_args()

    percents = [50,70,90,100,110,130,150] # This is hardcoded.

    #### Read net precision-recall ####

    nodePRs = {}
    edgePRs = {}

    # Read the interpolated PR values for each percent-run (only plot
    # aggregates for now)
    for percent in percents:

        for i in range(opts.n_samples):

            print opts
            nodeFileName = opts.node_PR_file%(percent,i)
            edgeFileName = opts.edge_PR_file%(percent,i)
            #nodeFileName = opts.outname + "aggregate" + "-sampleSet-" + str(percent) + "-run-" + str(i) + "-ksp" + str(opts.k) + "-interpedNodePR.txt"
            #edgeFileName = opts.outname + "aggregate" + "-sampleSet-" + str(percent) + "-run-" + str(i) + "-ksp" + str(opts.k) + "-interpedEdgePR.txt"

            edgePR = []
            for line in open(edgeFileName, 'r').readlines():
                items = line.split('\t')
                edgePR.append((float(items[0]), float(items[1])))
            nodePR = []
            for line in open(nodeFileName, 'r').readlines():
                items = line.split('\t')
                nodePR.append((float(items[0]), float(items[1])))

            nodePRs[(percent, i)] = np.array(nodePR)
            edgePRs[(percent, i)] = np.array(edgePR)


    # Merge the precision-recalls to a median to prepare for plotting
    medianNodePRs = {}
    medianEdgePRs = {}

    meanNodePRs = {}
    meanEdgePRs = {}

    stdDevNodePRs = {}
    stdDevEdgePRs = {}

    MADNodePRs = {}
    MADEdgePRs = {}


    for percent in percents:

        print("Merging PR's for percent = " + str(percent))

        medianNodePR = []
        medianEdgePR = []

        meanNodePR = []
        meanEdgePR = []

        stdDevNodePR = []
        stdDevEdgePR = []

        MADNodePR = []
        MADEdgePR = []

        nPts = 1000 # 1000 is hardcoded in
        recVals = np.linspace(0,1,nPts)

        for j in range(nPts):

            nodeValsList = []
            edgeValsList = []

            for i in range(opts.n_samples):

                nodePR = nodePRs[(percent, i)]
                if j < len(nodePR):
                    nodeValsList.append(nodePR[j,1])

                edgePR = edgePRs[(percent, i)]
                if j < len(edgePR):
                    edgeValsList.append(edgePR[j,1])

            if len(nodeValsList) > 0:
                median = np.median(nodeValsList)
                medianNodePR.append((recVals[j], np.median(nodeValsList)))
                MADNodePR.append((recVals[j], np.median(np.abs(nodeValsList - median))))

                meanNodePR.append((recVals[j], np.mean(nodeValsList)))
                stdDevNodePR.append((recVals[j], np.std(nodeValsList)))

            if len(edgeValsList) > 0:
                median = np.median(edgeValsList)
                medianEdgePR.append((recVals[j], median))
                MADEdgePR.append((recVals[j], np.median(np.abs(edgeValsList - median))))

                meanEdgePR.append((recVals[j], np.mean(edgeValsList)))
                stdDevEdgePR.append((recVals[j], np.std(edgeValsList)))

        medianNodePRs[percent] = np.array(medianNodePR)
        medianEdgePRs[percent] = np.array(medianEdgePR)

        meanNodePRs[percent] = np.array(meanNodePR)
        meanEdgePRs[percent] = np.array(meanEdgePR)

        stdDevNodePRs[percent] = np.array(stdDevNodePR)
        stdDevEdgePRs[percent] = np.array(stdDevEdgePR)

        MADNodePRs[percent] = np.array(MADNodePR)
        MADEdgePRs[percent] = np.array(MADEdgePR)


    #### Plot! ####

    orderedPercents = sorted(percents, key=lambda p: -abs(100-p))

    print("Making big plots")

    # Define colors
    darkBlue = '#000099'
    lightBlue = '#3399FF'
    colorDict = {50:'#0000FF', 70:'#0066FF', 90:'#3399FF', 100:'#000000', 110:'#00CC00', 130:'#006600', 150:'#003300'}

    plotDir = opts.outname

    # Make a large plot with the median for each percentage
    # Nodes
    names = []
    for percent in orderedPercents:
        medianPR = medianNodePRs[percent]
        plt.plot(medianPR[:,0], medianPR[:,1], lw=3, c=colorDict[percent], label = (percent-100))
        names.append(str((percent-100)))
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.xlabel('Recall')
    plt.ylabel('Precision')

    # Make sure the entries appear in order in the legend
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend([handles[labels.index(str(p-100))] for p in percents], [str(p-100) for p in percents])

    plt.title(opts.plot_median_title%("Node"))

    plt.savefig(plotDir + "Node-Medians.png")
    plt.savefig(plotDir + "Node-Medians.pdf")
    if(opts.show_plots):
        plt.show()
    plt.close()

    # Edges
    names = []
    for percent in orderedPercents:
        medianPR = medianEdgePRs[percent]
        plt.plot(medianPR[:,0], medianPR[:,1], lw=3, c=colorDict[percent], label = (percent-100))
        names.append(str((percent-100)))
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(opts.plot_median_title%("Edge"))

    # Make sure the entries appear in order in the legend
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend([handles[labels.index(str(p-100))] for p in percents], [str(p-100) for p in percents])

    plt.savefig(plotDir + "Edge-Medians.png")
    plt.savefig(plotDir + "Edge-Medians.pdf")
    if(opts.show_plots):
        plt.show()
    plt.close()

    # Make small plots for each percent
    print("Making small plots")
    # Nodes
    print("Nodes:")
    for percent in percents:

        print("\tpercent = " + str(percent))

        # Plot all of the light ones
        for i in range(opts.n_samples):

            precrec = nodePRs[(percent, i)]
            plt.plot(precrec[:,0], precrec[:,1], lw=3, c=lightBlue)

        # Plot the darker median
        precrec = medianNodePRs[percent]
        plt.plot(precrec[:,0], precrec[:,1], lw=3, c=darkBlue)

        plt.xlim([0,1])
        plt.ylim([0,1])
        plt.xlabel('Recall')
        plt.ylabel('Precision')

        plt.title(opts.plot_title%("Node", percent - 100))

        plt.savefig(plotDir + "Node-All%d.png"%(percent))
        plt.savefig(plotDir + "Node-All%d.pdf"%(percent))
        if(opts.show_plots):
            plt.show()
        plt.close()

    # Edges
    print("Edges:")
    for percent in percents:

        print("\tpercent = " + str(percent))

        # Plot all of the light ones
        for i in range(opts.n_samples):

            precrec = edgePRs[(percent, i)]
            plt.plot(precrec[:,0], precrec[:,1], lw=3, c=lightBlue)

        # Plot the darker median
        precrec = medianEdgePRs[percent]
        plt.plot(precrec[:,0], precrec[:,1], lw=3, c=darkBlue)

        plt.xlim([0,1])
        plt.ylim([0,1])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title(opts.plot_title%("Edge", percent - 100))

        plt.savefig(plotDir + "Edge-All%d.png"%(percent))
        plt.savefig(plotDir + "Edge-All%d.pdf"%(percent))
        if(opts.show_plots):
            plt.show()
        plt.close()



    ## Write out some useful statistics in text form

    statsFileName = plotDir + "statistics.txt"
    print("Writing statistics to " + statsFileName)
    statsFile = open(statsFileName, 'w')

    statsFile.write("Statistics for EDGES: ")

    ind30 = np.abs(recVals-0.30).argmin()
    ind60 = np.abs(recVals-0.60).argmin()

    # For small k values, the results may not reach 30% or 60% recall
    # (which generates an exception in the code below).  If that
    # happens, just move on without printing anything.
    try:
        statsFile.write("At 30%: ")
        for percent in percents:
            statsFile.write("\tPercent = %d\t median@30 = %.2f MAD@30 = %.2f mean@30 = %.2f stdDev@30 = %.2f"%(percent, medianEdgePRs[percent][ind30][1], MADEdgePRs[percent][ind30][1], meanEdgePRs[percent][ind30][1], stdDevEdgePRs[percent][ind30][1]))
    except IndexError:
        print("PR did not reach 30%, nothing will be printed there.")
        pass

    statsFile.write("")
    statsFile.write("At 60%: ")
    try:
        for percent in percents:
            statsFile.write("\tPercent = %d\t median@60 = %.2f MAD@60 = %.2f mean@60 = %.2f stdDev@60 = %.2f"%(percent, medianEdgePRs[percent][ind60][1], MADEdgePRs[percent][ind60][1], meanEdgePRs[percent][ind60][1], stdDevEdgePRs[percent][ind60][1]))
    except IndexError:
        print("PR did not reach 60%, nothing will be printed there.")
        pass

    statsFile.write("")
    statsFile.write("As a LaTeX table:")

    try:
        for percent in percents:
            statsFile.write("\\hline %d & %.2f (%.2e) & %.2f (%.2e) & %.2f (%.2e) & %.2f (%.2e)\\\\"%((percent-100), medianEdgePRs[percent][ind30,1], MADEdgePRs[percent][ind30,1], meanEdgePRs[percent][ind30,1], stdDevEdgePRs[percent][ind30,1], medianEdgePRs[percent][ind60,1], MADEdgePRs[percent][ind60,1], meanEdgePRs[percent][ind60,1], stdDevEdgePRs[percent][ind60,1]))
    except IndexError:
        print("PR did not reach 60%, nothing will be printed there.")
        pass

    statsFile.close()

if __name__=='__main__':
    main(sys.argv)














