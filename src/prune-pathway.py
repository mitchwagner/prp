'''
Given a pathway in the usual PathLinker pipeline format (consisting 
of a node and an edge file): 

1) Create a network from the nodes and edges file

2) Prune nodes and edges that are not on s-t paths 

3) For each line in the original nodes and edges file, check to see if the
   node or edges survived pruning. If so, writing to a new file, keep that
   line

'''
from pathlib import Path

import external.utils.graph_algorithms.prune as prune
import external.pathlinker.PathLinker as pl 
import external.pathlinker.parse as pl_parse

def write_pruned_pathway(
        old_nodes_file_handle, old_edges_file_handle,
        new_nodes_file_handle, new_edges_file_handle):

    net = create_network_from_pathway_files(
        old_nodes_file_handle, old_edges_file_handle)

    old_nodes_file_handle.seek(0)
    sources = pl_parse.get_source_set(old_nodes_file_handle)

    old_nodes_file_handle.seek(0)
    targets = pl_parse.get_target_set(old_nodes_file_handle)

    edges = set(net.edges())
    nodes = set(net.nodes())

    copy = net.copy()
    
    prune.remove_nodes_not_on_s_t_path(
        net, sources, targets, rm_ends=False, method="reachability")

    prune.remove_nodes_not_on_s_t_path(
        copy, sources, targets, rm_ends=True, method="reachability")

    new_edges = set(net.edges())
    new_nodes = set(net.nodes())
    
    new_edges2 = set(copy.edges())
    new_nodes2 = set(copy.nodes())

    #print(edges.difference(copy))
    #print(nodes.difference(copy))
    
    print("Extra edges deleted:")
    print(new_edges.difference(new_edges2))
    print("\nExtra nodes deleted:")
    a = new_nodes.difference(new_nodes2)

    for edge in a:
        if edge in sources:
            print("source", edge)
        else:
            print("target", edge)

    write_pruned_nodes_file(copy, old_nodes_file_handle, new_nodes_file_handle)

    write_pruned_edges_file(copy, old_edges_file_handle, new_edges_file_handle)


def create_network_from_pathway_files(nodes_file_handle, edges_file_handle):
    nodes_file_handle.seek(0)
    edges_file_handle.seek(0)

    net = None

    net = pl.readNetworkFile(edges_file_handle)

    nodes = set()
    for line in nodes_file_handle:
        if not line.lstrip().startswith("#"):
            toks = line.split("\t")
            nodes.add(toks[0])

    for node in nodes:
        net.add_node(node)

    return net


def write_pruned_nodes_file(net, old_file_handle, new_file_handle):
    old_file_handle.seek(0)
    
    nodes = set(net.nodes())
    
    for line in old_file_handle:
        if not line.lstrip().startswith("#"):
            toks = line.split("\t")
            node = toks[0]
            if node in nodes:
                new_file_handle.write(line)
        else:
            new_file_handle.write(line)


def write_pruned_edges_file(net, old_file_handle, new_file_handle):
    old_file_handle.seek(0)

    edges = set(net.edges())

    for line in old_file_handle:
        if not line.lstrip().startswith("#"):
            toks = line.split("\t")
            edge = (toks[0], toks[1])
            if edge in edges:
                new_file_handle.write(line)
        else:
            new_file_handle.write(line)
    

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

    new_pathway_dir = \
        Path(".."
            ,"inputs"
            ,"interactions"
            ,"netpath-s-t-pruned-including-ends"
            ,"pathways"
            )

    new_pathway_dir.mkdir(parents=True, exist_ok=True)

    for name in pathway_names:
        print(name)
        old_node_file = Path(pathway_dir, name + "-nodes.txt")
        old_edge_file = Path(pathway_dir, name + "-edges.txt")
        new_node_file = Path(new_pathway_dir, name + "-nodes.txt")
        new_edge_file = Path(new_pathway_dir, name + "-edges.txt")

        with old_node_file.open('r') as onf, old_edge_file.open('r') as oef, \
             new_node_file.open('w') as nnf, new_edge_file.open('w') as nef:
            
            write_pruned_pathway(onf, oef, nnf, nef)
        print("---------------------------------------------")

if __name__ == "__main__":
    main()
