from pathlib import Path

from src.runners.reconstruction.AlgorithmRunner \
    import AlgorithmRunner

from src.runners.folds.NodeAndEdgeWithholdingFoldCreator \
    import NodeAndEdgeWithholdingFoldCreator


class NodeAndEdgeWithholdingRunner(AlgorithmRunner): 
    '''
    '''

    def get_name(self):
        return "node-and-edge-withholding evaluation"


    def get_details(self):
        return "Fraction Retained: Nodes = %f, Edges = %f" % (
            self.options["percent_nodes_to_keep"], 
            self.options["percent_edges_to_keep"])


    def get_fold_creator(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = NodeAndEdgeWithholdingFoldCreator(
            self.interactome, pathway, self.options)

        return fc


    def get_output_prefix(self):
        return Path("node-and-edge-deletion-v4", "nodes-%f-edges-%f-iter-%d" % (
            self.options["percent_nodes_to_keep"],
            self.options["percent_edges_to_keep"],
            self.options["iterations"]))
