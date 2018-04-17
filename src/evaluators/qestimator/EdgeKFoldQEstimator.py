from pathlib import Path 
from src.evaluators.qestimator.QEstimator import QEstimator

from src.fold_creators.EdgeKFoldFoldCreator \
    import EdgeKFoldFoldCreator 

class EdgeKFoldQEstimator(QEstimator):
    '''
    Implement the factory methods necessary for a QEstimator
    '''
    
    def get_fold_creator(self, pathway):
        fc = EdgeKFoldFoldCreator(
            self.interactome, pathway, self.options)

        return fc
   

    def get_output_prefix(self):
        return Path("edge-kfold", "%d-folds" % self.options["num_folds"])
