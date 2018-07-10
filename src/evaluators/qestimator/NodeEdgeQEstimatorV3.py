from pathlib import Path

from src.evaluators.qestimator.QEstimatorV3 import QEstimatorV3

from src.fold_creators.NodeAndEdgeWithholdingFoldCreatorV3 \
    import NodeAndEdgeWithholdingFoldCreatorV3

class NodeEdgeQEstimatorV3(QEstimatorV3):
    '''
    Implement the factory methods necessary for a QEstimator
    '''
    
    def get_fold_creator(self, pathway):
        fc = NodeAndEdgeWithholdingFoldCreatorV3(
            self.interactome, pathway, self.options)

        return fc
   

    def get_output_prefix(self):
        return Path("node-and-edge-deletion", "nodes-%f-edges-%f-iter-%d" % (
            self.options["percent_nodes_to_keep"],
            self.options["percent_edges_to_keep"],
            self.options["iterations"]))
