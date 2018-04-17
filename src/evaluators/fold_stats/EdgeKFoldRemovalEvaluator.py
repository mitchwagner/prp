from pathlib import Path

from src.evaluators.fold_stats.RemovalEvaluator import RemovalEvaluator

import src.fold_creators.EdgeKFoldFoldCreator as EdgeKFold

class EdgeKFoldRemovalEvaluator(RemovalEvaluator):
    '''
    Analyze the impact of removing edges according to folds defined
    by cross validation.
    '''

    def get_fold_creator(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''

        fc = EdgeKFold.EdgeKFoldFoldCreator(
            self.interactome, pathway, self.options)

        return fc


    def analyze_fold_creation(self, pathway):
        fc = self.get_fold_creator(pathway)

        num_folds = self.options["num_folds"]

        stats = [[] for i in range(num_folds)]

        pathway_obj = pathway.get_pathway_obj()
        initial_node_count = len(pathway_obj.get_nodes(data=False))
        initial_edge_count = len(pathway_obj.get_edges(data=False))

        folds = fc.get_training_folds()

        for i, fold in enumerate(folds):
            edges_after_deletion = len(fold[0])
            
            stats[i].append(pathway.name)
            stats[i].append(initial_node_count)
            stats[i].append(initial_edge_count)
            stats[i].append(edges_after_deletion)

        # Append a header to the list 
        header = ["pathway", "# nodes initial", "# edges initial",
                  "# edges after edge deletion"]

        stats.insert(0, header)

        return stats


    def get_output_prefix(self):
        return Path("edge-kfold", "%d-folds" % self.options["num_folds"])
