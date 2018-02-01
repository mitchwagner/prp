'''
Use abstraction of pathway and interactome to create an interactome that 
with all edges from the pathway and all edges from the interactome and all
nodes from all of the pathways, resulting in a giant franken-teractome
'''

from pathlib import Path

import external.utils.pathway.pathway as pw
import external.utils.interactome.interactome as inter


def main():
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

    pathway_dir = \
        Path(".."
            ,"inputs"
            ,"interactions"
            ,"netpath"
            ,"pathways"
            )

    interactome_file = \
        Path(".."
            ,"inputs"
            ,"interactomes"
            ,"human"
            ,"pathlinker-signaling-children-reg-weighted.txt"
            )

    interactome = None

    with interactome_file.open('r') as f: 
        interactome = inter.parse_csbdb_interactome_file(f)

    num_missing = {}
    for name in pathway_names:
        counter = 0
        counter2 = 0
        print(name)
        print("--------------")
        edge_file = Path(pathway_dir, name + "-edges.txt")
        node_file = Path(pathway_dir, name + "-nodes.txt")

        with edge_file.open('r') as f1, node_file.open('r') as f2:
            pathway = pw.parse_csbdb_pathway_file(f1, f2)

            for edge in pathway.get_edges(data=False):
                counter2 += 1
                if edge not in interactome.get_edges(data=False):
                    counter +=1
                    interactome.add_edge(edge[0], edge[1], {"weight":.75})
                    print(edge)

        num_missing[name] = (counter, counter2)

    for key in num_missing:
        print("pathway: %s num_missing: %d pathway_size: %d" % (key, num_missing[key][0], num_missing[key][1]))

    with Path("additional.txt").open('w') as f:
        inter.write_csbdb_interactome_file(interactome, f)


if __name__ == "__main__":
    main()
