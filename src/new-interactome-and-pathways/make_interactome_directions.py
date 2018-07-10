'''
Use the evidence file to construct a map of edge to its directionality.
For example, (a, b) might map to True, and (c, b) might map to False
'''

import sys

edge_directedness = {}

evidence_file = sys.argv[1]
output_file = sys.argv[2]

# Read in the evidence file. 
# For every edge, get its directedness.
#    If directed, add it to the dictionary as such and delete the reverse
#    direction, if it exists
#
#   If undirected, add this line and its reverse to the dictionary, unless
#   this line or the reverse already exist as keys with the value "directed"
#
#with open('2018_01pathlinker-no-kegg-spike-extra-evidence.tsv') as f:
# with open('/home/mitchw94/Desktop/svnrepo/src/python/Interactomes/Human/evidence.txt') as f:
#with open('directionality-test.tsv', 'r') as f:
with open(evidence_file, 'r') as f:
    for line in f:
        if line.startswith("#"):
            continue
        else:
            toks = [x.strip() for x in line.split("\t")]
            edge = (toks[0], toks[1])
            rev_edge = (toks[1], toks[0])
            directedness = toks[2]

            if directedness == 'True':
                edge_directedness[edge] = True 

                # If the other direction is in the dictionary and it is 
                # undirected, remove it
                if rev_edge in edge_directedness:
                    if edge_directedness[rev_edge] == False:
                        edge_directedness.pop(rev_edge)
            
            else:
                if edge in edge_directedness or rev_edge in edge_directedness:
                    # Don't do anything.
                    # - If edge is True, we shouldn't do anything
                    # - If rev_edge is True, we shouldn't do anything
                    # - If they are both False (only other possibility), then
                    #   doing anything else is redundant
                    None

                else:
                    edge_directedness[edge] = False
                    edge_directedness[rev_edge] = False


# Write out the directionality file
# with open('interactome-directions.tsv', 'w') as f:
# with open('kegg-only-directions.tsv', 'w') as f:
with open(output_file, 'w') as f:
    for edge in edge_directedness:
        f.write(
            edge[0] + "\t" + 
            edge[1] + "\t" + 
            str(edge_directedness[edge]) + "\n")
