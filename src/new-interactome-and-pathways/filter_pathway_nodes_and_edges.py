'''
There are a number of nodes that are incorrect in the SBML-parsed nodes and 
edges files for pathways. The intent of this script is to filter those 
nodes and edges files by removing lines where the nodes are not in a given
evidence file.

We will subsequently use these nodes and edges files in the pathway, AFTER
S-T pruning.
'''

from pathlib import Path

pathway_dir = Path("..", "..", "inputs", "interactions",
    "netpath-s-t-pruned-updated-receptors", "pathways")


# First, read the evidence file to get a list of nodes in the interactome

evidence_file = Path("2018_01pathlinker-no-kegg-spike.tsv")

# Read in the evidence file and get a set of nodes

nodes = set()
counter = 0
with evidence_file.open('r') as f:
    for line in f:
        counter += 1
        # Ignore comments
        if line.startswith("#"):
            continue
        else:
            toks = [col.strip() for col in line.split("\t")]
            tail = toks[0]
            head = toks[1]
            nodes.add(tail)
            nodes.add(head)

pathway_names = \
    ["BDNF"
    ,"EGFR1"
    ,"IL1"
    ,"IL2"
    ,"IL3"
    ,"IL6"
    ,"IL-7"
    ,"KitReceptor"
    ,"Leptin"
    ,"Prolactin"
    ,"RANKL"
    ,"TCR"
    ,"TGF_beta_Receptor"
    ,"TNFalpha"
    ,"Wnt"
    ]

for name in pathway_names:
    print(name)
    edge_file = Path(pathway_dir, name + "-edges.txt")
    node_file = Path(pathway_dir, name + "-nodes.txt")

    filtered_edge_file = Path("filtered-pathways", name + "-edges.txt")
    filtered_node_file = Path("filtered-pathways", name + "-nodes.txt")

    filtered_edge_file.parent.mkdir(exist_ok=True, parents=True)
    
    edges_removed = 0
    nodes_removed = 0
    with edge_file.open('r') as f1, filtered_edge_file.open('w') as f2:
        for line in f1:
            # Ignore comments
            if line.startswith("#"):
                f2.write(line)
            else:
                toks = [x.strip() for x in line.split("\t")]
                tail = toks[0]
                head = toks[1]

                if tail in nodes and head in nodes:
                    f2.write(line)
                else:
                    print("    Skipping:", str((tail, head)))
                    edges_removed += 1

    print("    ", edges_removed, "edges removed")

    with node_file.open('r') as f1, filtered_node_file.open('w') as f2:
        for line in f1:
            # Ignore comments
            if line.startswith("#"):
                f2.write(line)
            else:
                toks = [x.strip() for x in line.split("\t")]
                node = toks[0]
                if node in nodes:
                    f2.write(line)
                else:
                    print("    Skipping:", str(node))
                    nodes_removed += 1
                

    print("    ", nodes_removed, "nodes removed")
