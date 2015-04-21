# Create a master list of all TFs and RECs in the PPI by scanning many pathway lists

# Nicholas Sharp - nsharp3@vt.edu

import os
import sys
import subprocess
from optparse import OptionParser

def main(args):

    usage = 'create-master-TF-REC-lists.py PATHWAYS...'
    parser = OptionParser(usage=usage)

    parser.add_option('-o','--outname',default='master',\
                          help='Prepend to results. Default: "master"')

    parser.add_option('', '--force-recalc', action='store_true', default=False, help='Recompute everything, instead of following the usual policy of note recomputin existing files.')

    # parse the command line arguments
    (opts, args) = parser.parse_args()

    PATHWAYS = args
    print("Using pathways " + str(PATHWAYS))

    ## Read in pathways and percents
    allTFs = set()
    allRECs = set()

    # Read from the real file
    for pathway in PATHWAYS:
        pFile = open(pathway,'r')

        for line in pFile.readlines():
            if line[0] != '#':
                items = line.split('\t')
                if(items[1] == 'receptor'):
                    allRECs.add(items[0])
                elif(items[1] == 'tf'):
                    allTFs.add(items[0])

    print("Found %d TFs and %d RECs"%(len(allTFs), len(allRECs)))

    # Write the result

    # TFs
    outTFFile = opts.outname + "TFs.txt"
    print("Will write TFs to " + outTFFile)

    if(os.path.isfile(outTFFile) and not opts.force_recalc):
        print("Skipping writing of master TFs because file already exists")
    else:
        outF = open(outTFFile, 'w')
        for x in allTFs:
            outF.write(x+"\n")
        outF.close()

    # RECs
    outRECFile = opts.outname + "RECs.txt"
    print("Will write RECs to " + outRECFile)

    if(os.path.isfile(outRECFile) and not opts.force_recalc):
        print("Skipping writing of master RECs because file already exists")
    else:
        outF = open(outRECFile, 'w')
        for x in allRECs:
            outF.write(x+"\n")
        outF.close()



if __name__=='__main__':
    main(sys.argv)
