from pathlib import Path

from src.evaluators.reconstruction.AlgorithmEvaluator import AlgorithmEvaluator

from src.fold_creators.EmpiricalEdgeSamplingFoldCreator \
    import EmpiricalEdgeSamplingFoldCreator

class EmpiricalEdgeSamplingEvaluator(AlgorithmEvaluator): 
    '''
    Note: we don't really deletet edge but the percentages here
    correspond to empirical values for node + edge deletion. Hence,
    the names here reflect that
    '''

    def get_name(self):
        return "empirical edge sampling evaluation"


    def get_details(self):
        return "Corresponding Fractions: Nodes = %f, Edges = %f" % (
            self.options["percent_nodes_to_keep"], 
            self.options["percent_edges_to_keep"])


    def get_fold_creator(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = EmpiricalEdgeSamplingFoldCreator(
            self.interactome, pathway, self.options)

        return fc


    def get_output_prefix(self):
        return Path("empirical-edge-sampling", "nodes-%f-edges-%f-iter-%d" % (
            self.options["percent_nodes_to_keep"],
            self.options["percent_edges_to_keep"],
            self.options["iterations"]))
