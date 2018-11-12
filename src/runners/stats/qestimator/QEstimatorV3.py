from pathlib import Path

import numpy as np
import networkx as nx

from src.runners.Runner import Runner 

import src.runners.folds.FoldCreator as FoldCreator
from src.external.pathlinker import PathLinker as pl

class QEstimatorV3(Runner):
    '''
    Estimate q. For each pathway, for each fold, get all the shortest 
    paths using PathLinker. Count sequences of x's and find the average
    length of a sequence of x's in these paths.
    '''

    def __init__(
            self, interactome, pathway_collection, options={}):
        '''
        :param interactome: on-disk interactome object
        :param pathway_collection: PathwayCollection object
        :param algorithms: list of RankingAlgorithms
        :param options: map of options for the evaluator
        '''
        self.interactome = interactome
        self.pathway_collection = pathway_collection
        self.options = options


    def get_fold_creator(self, pathway):
        '''
        Factory method intended to supply a FoldCreator
        '''
        raise NotImplementedError()


    def get_output_prefix(self):
        '''
        Factory method intended to supply a name for writing output to
        '''
        raise NotImplementedError()
   

    def run(self, output_dir=Path(), purge_results=False):

        output_dir = Path(output_dir, "q-estimator") 

        # Read in the interactome

        interactome_net = None

        with self.interactome.path.open('r') as f:
            interactome_net = pl.readNetworkFile(f)

        number_of_xs_per_pathway = []

        for pathway in self.pathway_collection.pathways:
            print("Pathway: %s" % pathway.name)

            # Comment in the one we want!
            fc = self.get_fold_creator(pathway)

            # Make a pathway_obj for convenient access to sources/targets 
            pathway_obj = pathway.get_pathway_obj()

            # Make a network from the pathway 
            pathway_net = pathway_obj.get_net_from_pathway()

            targets = pathway_obj.get_tfs(data=False)
    
            '''
            This should no longer be necessary

            # Remove edges from the pathway net if they are not in the
            # interactome. Also, give them weights from the interactome
            filtered_edges = FoldCreator.get_filtered_pathway_edges(
                pathway_obj, self.interactome)

            for edge in pathway_net.edges():
                if edge not in filtered_edges:
                    pathway_net.remove_edge(edge[0], edge[1])
                else:
                    pathway_net[edge[0]][edge[1]]["weight"] = 1
                    pathway_net[edge[0]][edge[1]]["ksp_weight"] = 1
            '''

            '''
            This is the only component that should need replacing if I 
            do things correctly

            training_folds = fc.get_training_folds()
            test_folds = fc.get_test_folds()
            '''

            training_folds = []
            test_folds = []

            #folds = zip(fc.create_positive_folds(), fc.create_negative_folds())

            pos_folds, neg_folds = fc.create_folds_initial()

            folds = zip(pos_folds, neg_folds)

            for c, (pos_fold, neg_fold) in enumerate(folds):
                pos_dir_train, pos_dir_test, pos_undir_train, pos_undir_test = \
                    pos_fold

                neg_dir_train, neg_dir_test, neg_undir_train, neg_undir_test = \
                    neg_fold

                pos_merged_train = list(set(pos_dir_train).union(set(pos_undir_train)))
                pos_merged_test = list(set(pos_dir_test).union(set(pos_undir_test)))

                neg_merged_train = list(set(neg_dir_train).union(set(neg_undir_train)))
                neg_merged_test = list(set(neg_dir_test).union(set(neg_undir_test)))

                training_folds.append((pos_merged_train, neg_merged_train))
                test_folds.append((pos_merged_test, neg_merged_test))

            # TODO: Clean the above + variable names up. We are
            # interested in the length of the path, not in the 
            # number of x's.
            number_of_xs = []

            for i, _ in enumerate(training_folds):
                print("Fold %d" % i)

                # First, create a lookup map for quick determination of an
                # edge's label

                training_positives = training_folds[i][0]
                test_positives = test_folds[i][0]
                
                # List of sources == head nodes in training positives
                # So, we'll make this simultaneously with the labeling

                sources = []
                
                for edge in training_positives:
                    sources.append(edge[1])

                for source in sources:
                    for target in targets:
                        x_count = 0

                        try:
                            path = \
                                nx.shortest_path(pathway_net, source, target)
                            # Plus one because during RWER we first
                            # slide down from u to v?
                            number_of_xs.append(len(path) + 1)
                        except:
                            number_of_xs.append(99999999999)
                            continue
                            
            number_of_xs_per_pathway.append(number_of_xs) 

        self.write_stats(number_of_xs_per_pathway, output_dir)


    def write_stats(self, number_of_xs, output_dir):
        concatenated_list = []
        for i, pathway in enumerate(self.pathway_collection.pathways):
            ls = number_of_xs[i]

            name = pathway.name
            outfile = self.get_pathway_outfile(output_dir, pathway)
            outfile.parent.mkdir(parents=True, exist_ok=True)

            avg = sum(ls)/len(ls)
            median = np.median(ls)
            std = np.std(ls)
            q = self.estimate_q(median)

            with outfile.open('w') as f:
                f.write("Pathway: %s\n" % name)
                f.write("avg: %f\n" % avg)
                f.write("median: %f\n" % median)
                f.write("stddev: %f\n" % std)
                f.write(str(q))

            concatenated_list = concatenated_list + ls

        avg = sum(concatenated_list) / len(concatenated_list)
        median = np.median(concatenated_list)
        std = np.std(concatenated_list)
        q = self.estimate_q(median)

        outfile = self.get_overall_outfile(output_dir)
            
        with outfile.open('w') as f:
            f.write("Overall: %s\n" % name)
            f.write("avg: %f\n" % avg)
            f.write("median: %f\n" % median)
            f.write("stddev: %f\n" % std)
            f.write(str(q))


    def estimate_q(self, x):
        '''
        x = (1 - q) / q
        qx = 1 - q
        qx + q = 1
        q(x + 1) = 1
        q = (1 / (x + 1))
        '''
        q = (1 / (x + 1))

        return q


    def get_pathway_outfile(self, output_dir, pathway):
        return Path(
            output_dir, self.get_output_prefix(), pathway.name + "-stats.txt")


    def get_overall_outfile(self, output_dir):
        return Path(
            output_dir, self.get_output_prefix(), "overall-stats.txt")
