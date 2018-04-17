from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np


from src.evaluators.Evaluator import Evaluator
import src.external.pathlinker.PathLinker as pl

class PathwayStats(Evaluator):
    '''
    Runnable analysis that outputs several statistics about the interactomes.
    '''

    def __init__(self, interactome, pathway_collection):
        self.interactome = interactome
        self.pathway_collection = pathway_collection


    def get_name(self):
        '''
        Return a name for the evaluation performed by this Evaluator
        '''
        return "Pathway Statistics Evaluation"


    def get_interactome_net(self):
        with self.interactome.path.open('r') as f:
            return pl.readNetworkFile(f)


    def get_output_dir(self):
        return Path("pathway-stats")


    def write_summary_statistics(self, output_dir):
        None


    def plot_edge_weight_distributions(self, output_dir):
        net = self.get_interactome_net()

        for pathway in self.pathway_collection.pathways:
            
            pathway_obj = pathway.get_pathway_obj()
            pathway_edges = set(pathway_obj.get_edges(data=False))

            weights = [edge[2]["weight"] for edge in net.edges(data=True) if 
                (edge[0], edge[1]) in pathway_edges]
            
            fig, ax = self.init_edge_weight_figure(pathway)
            
            ax.hist(weights, bins=np.arange(0, 1, .05), histtype='bar')
            self.save_edge_weight_distribution_fig(fig, pathway, output_dir)


    def init_edge_weight_figure(self, pathway):

        fig, ax = plt.subplots()
        ax.set_title("Edge Weight Distribution (%s)" % pathway.name)
        ax.set_xlabel("Edge weight (bin size = .05)")
        ax.set_ylabel("# of edges in bin")

        return fig, ax


    def save_edge_weight_distribution_fig(self, fig, pathway, output_dir):

        outfile_pdf = \
            Path(output_dir, self.get_output_dir(), self.interactome.name,
                pathway.name + "-edge-weight-distribution.pdf")

        outfile_png = \
            Path(output_dir, self.get_output_dir(), self.interactome.name,
                pathway.name + "-edge-weight-distribution.png")

        outfile_pdf.parent.mkdir(parents=True, exist_ok=True)

        fig.tight_layout()
        fig.savefig(str(outfile_pdf))
        fig.savefig(str(outfile_png))


    # TODO: purge_results is not implemented
    def run(self, output_dir=Path(), purge_results=False):
        '''
        Gathers statistics on the following:
        1) The number of nodes and edges in each pathway

        2) The number of nodes and edges from each pathway that exist in a
           given interactome

        3) The edge-weight distribution in each pathway, using the weights
           for a given interactome
        '''
        self.write_summary_statistics(output_dir)
        self.plot_edge_weight_distributions(output_dir)


"""

    def pathway_subset_analysis(self):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                results = []
                for pathway in pathway_collection.pathways:

                    #self.get_pathway_specific_interactome_file_path(
                    #    interactome, pathway)
                
                    # Read in the interactome from an edgelist
                    interactome_net = None
                    with specific_interactome.open('r') as f:
                        interactome_net = pl.readNetworkFile(f) 
                    
                    # Pathway node and edge files
                    node_file = \
                        pathway_collection.get_pathway_nodes_file(pathway)

                    edge_file = \
                        pathway_collection.get_pathway_edges_file(pathway)

                    # Read the pathway network from the edge file
                    net = None
                    with edge_file.open('r') as f:
                        net = pl.readNetworkFile(f) 
                    
                    # Read the pathway nodes
                    nodes = pathway_collection.get_nodes_from_pathway_nodes_file(
                        pathway)

                    # Add nodes to pathway if for whatever reason they are not
                    # on the edgelist
                    for node in nodes:
                        net.add_node(node)

                    sources = None 
                    with node_file.open('r') as f: 
                        sources = pl_parse.get_source_set(f)

                    targets = None
                    with node_file.open('r') as f:
                        targets = pl_parse.get_target_set(f)

                    interactome_nodes = set(interactome_net.nodes())
                    interactome_edges = set(interactome_net.edges())
                    
                    pathway_nodes = set(net.nodes())
                    pathway_edges = set(net.edges())

                    results.append((
                        str(interactome),
                        pathway_collection.name,
                        pathway,
                        len(pathway_nodes),
                        len(pathway_nodes.intersection(interactome_nodes)),
                        len(pathway_nodes.intersection(interactome_nodes)) / len(pathway_nodes),
                        len(sources),
                        len(sources.intersection(interactome_nodes)),
                        len(targets), 
                        len(targets.intersection(interactome_nodes)),
                        len(pathway_edges),
                        len(pathway_edges.intersection(interactome_edges)),
                        len(pathway_edges.intersection(interactome_edges)) / len(pathway_edges)
                        ))

                table_file = Path(
                    "outputs",
                    "other",
                    "subset-analysis",
                    pathway_collection.name,
                    "table.txt")

                table_file.parent.mkdir(parents=True, exist_ok=True)

                with table_file.open('w') as f:
                    f.write(
                        "Interactome\t"
                        "Pathway Collection\t"
                        "Pathway\t"
                        "Pathway Node #\t"
                        "Nodes in Interactome\t"
                        "Fraction in Interactome\t"
                        "# Sources\t"
                        "# Sources in Interactome\t"
                        "# Targets\t"
                        "# Targets in Interactome\t"
                        "Pathway Edge #\t"
                        "Edges in Interactome\t"
                        "Fraction in Interactome\n")

                    for result in results:
                        f.write("\t".join([str(elem) for elem in result]))
                        f.write("\n")
    '''
'''
def pathway_edge_weight_histograms(self):

    Create histograms, per pathway, of the weights of edges in that
    pathway.

    for interactome in self.input_settings.interactomes:
        # TODO: This is not a pathway specific interactome?
        specific_interactome = interactome.path
        for pathway_collection in self.input_settings.pathway_collections: 
            for pathway in pathway_collection.pathways:
                fig, ax = plt.subplots()

                ax.set_title(pathway.name)
                ax.set_xlabel("Edge weight (bin size = .05)")
                ax.set_ylabel("# of edges in bin")

                pathway_obj = pathway.get_pathway_obj()
                edges = pathway_obj.get_edges(data=False)

                interactome = None
                with specific_interactome.open('r') as f:
                    interactome = pl.readNetworkFile(f) 

                final_edges = set(edges).intersection(
                    set(interactome.edges()))
            
                weights = [interactome[x[0]][x[1]]["weight"] 
                    for x in final_edges] 


                ax.hist(weights, bins=np.arange(0,1,.05), fc=(0,0,1,.5))

                out = Path(
                    "outputs",
                    "other",
                    "edge-weight-distribution",
                    pathway.name + "-histogram.png")

                out.parent.mkdir(parents=True, exist_ok=True)

                fig.savefig(str(out))
"""
