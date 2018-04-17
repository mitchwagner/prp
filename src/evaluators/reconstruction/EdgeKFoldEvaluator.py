from pathlib import Path

from src.evaluators.reconstruction.AlgorithmEvaluator import AlgorithmEvaluator

from src.fold_creators.EdgeKFoldFoldCreator import EdgeKFoldFoldCreator

class EdgeKFoldEvaluator(AlgorithmEvaluator): 


    def get_name(self):
        return "edge-withholding cross-validation"


    def get_details(self):
        return "Number folds: %d" % self.options["num_folds"]


    def get_fold_creator(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = EdgeKFoldFoldCreator(
            self.interactome, pathway, self.options)

        return fc


    def get_output_prefix(self):
        return Path("edge-cv", "%d-folds" % self.options["num_folds"])
