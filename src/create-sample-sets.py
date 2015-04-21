# Create up-sampled and down-sampled node sets

# Nicholas Sharp - nsharp3@vt.edu

import os
import sys
import subprocess
from optparse import OptionParser
from math import ceil
import random as rand

def main(args):

    #usage = 'SubsamplingPipeline.py PATHWAYS_LIST PERCENTS_LIST  [options]\n'
    usage = 'create-sample-sets.py NODE_FILE'
    parser = OptionParser(usage=usage)

    parser.add_option('-o','--outname',default='sampled-set',\
                          help='Prepend to results. Default: "sampled-set"')

    parser.add_option('-n', '', action='store', type='int', default=1,\
                help='The number of sets to generate (default=1)')

    parser.add_option('-p', '', action='store', type='int', default=100,\
                help='The percent (a positive integer around 100) at which to sample (<100 for subsampling, >100 for supersampling, 100 for no sampling) (default=1000)')

    parser.add_option('', '--force-recalc', action='store_true', default=False, help='Recompute everything, instead of following the usual policy of note recomputin existing files.')

    parser.add_option('','--ppi',default=None, help='Location of ppi file')
    parser.add_option('','--all-TFs',default=None, help='Location of global TF list')
    parser.add_option('','--all-Recs',default=None, help='Location of global receptor list')

    # parse the command line arguments
    (opts, args) = parser.parse_args()

    NODEFILE = args[0]

    print("Using nodefile: " + str(NODEFILE))
    print("Using percent: " + str(opts.p))

    ## Check and print all file locations
    # ppi
    if(opts.ppi == None):
        print("Error! No ppi file given")
    print("Using ppi file " + opts.ppi)

    # sets of TF and Recs
    if(opts.all_TFs == None):
        print("Error! No global TF file given")
    print("Using TF file " + opts.all_TFs)
    if(opts.all_Recs == None):
        print("Error! No global receptor file given")
    print("Using rec file " + opts.all_Recs)

    # directory to use for all results
    print("Using output directory " + opts.outname)

    # TODO implement "force-recalc" logic here

    realRecSet = set()
    realTFSet = set()
    realNoneSet = set()
    bgRecSet = set()
    bgTFSet = set()

    # Read from the real file
    realF = open(NODEFILE, 'r')
    for line in realF.readlines():
        if line[0] != '#':
            items = line.split('\t')
            if(items[1] == 'receptor'):
                realRecSet.add(items[0])
            elif(items[1] == 'tf'):
                realTFSet.add(items[0])
            else:
                realNoneSet.add(items[0])

    # Read from the bg rec file
    bgRecF = open(opts.all_Recs, 'r')
    for line in bgRecF.readlines():
        if line[0] != '#':
            bgRecSet.add(line.strip())
    # Read from the bg tf file
    bgTFF = open(opts.all_TFs, 'r')
    for line in bgTFF.readlines():
        if line[0] != '#':
            bgTFSet.add(line.strip())

    # Remove the real values from the bgSet to avoid sampling them
    bgRecSet.difference_update(realRecSet)
    bgTFSet.difference_update(realTFSet)

    print("After differencing:")
    print("    The real rec set contains " + str(len(realRecSet)) + " items")
    print("    The background rec set contains " + str(len(bgRecSet)) + " items")
    print("    The real TF set contains " + str(len(realTFSet)) + " items")
    print("    The background TF set contains " + str(len(bgTFSet)) + " items")

    baseFileName = opts.outname + "sample_" + str(opts.p) + "-"

    percent = opts.p
    if percent >= 100:
        sampleRecSize = len(realRecSet) * percent / 100
        sampleTFSize = len(realTFSet) * percent / 100
    else:
        sampleRecSize = int(ceil(len(realRecSet) * percent / 100.0))
        sampleTFSize = int(ceil(len(realTFSet) * percent / 100.0))

    print("NOTE because truncating arithmetic makes the sampling factor have seeminly different results for upsampling and downsampling, we adjust to make downsampling round up. Be sure the values below correspond to your understanding.")
    print("Sampling rec goes %d --> %d"%(len(realRecSet), sampleRecSize))
    print("Sampling tf goes %d --> %d"%(len(realTFSet), sampleTFSize))


    print("Writing " + str(opts.n) + " sampled sets to base name " + baseFileName)
    for i in range(opts.n):

        # Open the file
        outFileName = baseFileName + str(i) + "-nodes.txt"
        
        if(os.path.isfile(outFileName) and not opts.force_recalc):
            print("Skipping file " + outFileName + " because it already exists")
            continue

        outF = open(outFileName, 'w')

        # Sample the receptors and transcription factors
        if percent >= 100:
            outRecSet = realRecSet.copy()
            outRecSet.update(rand.sample(bgRecSet, sampleRecSize - len(realRecSet)))
            outTFSet = realTFSet.copy()
            outTFSet.update(rand.sample(bgTFSet, sampleTFSize - len(realTFSet)))
        else:
            outRecSet = rand.sample(realRecSet, sampleRecSize)
            outTFSet = rand.sample(realTFSet, sampleTFSize)

        # Write the result   
        # TODO carry the names around instead of printing N/A
        for x in outRecSet:
            outF.write(x + "\treceptor\tN/A\n")
        for x in outTFSet:
            outF.write(x + "\ttf\tN/A\n")
        for x in realNoneSet:
            outF.write(x + "\tnone\tN/A\n")

        outF.close()







if __name__=='__main__':
    main(sys.argv)
