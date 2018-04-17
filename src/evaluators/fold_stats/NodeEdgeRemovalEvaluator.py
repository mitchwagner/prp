from pathlib import Path

from src.evaluators.fold_stats.RemovalEvaluator import RemovalEvaluator
import src.fold_creators.NodeAndEdgeWithholdingFoldCreator as NodeEdge

class NodeEdgeRemovalEvaluator(RemovalEvaluator):
    '''
    Analyze the impact of removing a certain percentage of nodes/edges
    '''

    def get_fold_creator(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = NodeEdge.NodeAndEdgeWithholdingFoldCreator(
            self.interactome, pathway, self.options)

        return fc


    def analyze_fold_creation(self, pathway):

        fc = self.get_fold_creator(pathway)

        copies = fc.get_pathway_copies()

        stats = [[] for i in range(len(copies))]

        initial_node_count = len(copies[0].nodes())
        initial_edge_count = len(copies[0].edges())

        for row in stats:
            row.append(pathway.name)
            row.append(initial_node_count)
            row.append(initial_edge_count)

        fc.delete_pathway_node_percentage(copies)

        for i, copy in enumerate(copies):
            stats[i].append(len(copy.nodes()))
            stats[i].append(len(copy.edges()))

        fc.delete_pathway_edge_percentage(copies)

        for i, copy in enumerate(copies):
            stats[i].append(len(copy.edges()))


        # At the end here, I will append the header
        header = ["pathway", "# nodes initial", "# edges initial",
                  "# nodes after node deletion", "# edges after node deletion",
                  "# edges after edge deletion"]

        stats.insert(0, header)

        return stats


    def get_output_prefix(self):
        return Path("node-and-edge-deletion", "nodes-%f-edges-%f-iter-%d" % (
            self.options["percent_nodes_to_keep"],
            self.options["percent_edges_to_keep"],
            self.options["iterations"]))
