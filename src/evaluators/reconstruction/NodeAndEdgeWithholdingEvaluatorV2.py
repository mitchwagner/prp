from pathlib import Path

from src.evaluators.reconstruction.AlgorithmEvaluatorV2 \
    import AlgorithmEvaluatorV2

from src.fold_creators.NodeAndEdgeWithholdingFoldCreatorV2 \
    import NodeAndEdgeWithholdingFoldCreatorV2

class NodeAndEdgeWithholdingEvaluatorV2(AlgorithmEvaluatorV2): 


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
        fc = NodeAndEdgeWithholdingFoldCreatorV2(
            self.interactome, pathway, self.options)

        return fc


    def get_output_prefix(self):
        return Path("node-and-edge-deletion", "nodes-%f-edges-%f-iter-%d" % (
            self.options["percent_nodes_to_keep"],
            self.options["percent_edges_to_keep"],
            self.options["iterations"]))
