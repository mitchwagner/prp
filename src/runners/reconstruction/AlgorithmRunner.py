import time
import itertools
from pathlib import Path

import scipy as sp
import numpy as np
import networkx as nx

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors

import src.external.pathlinker.parse as pl_parse
import src.external.utils.precision_recall.precision_recall as precrec

import src.input_utilities as iu

import src.algorithms.RankingAlgorithm as RankingAlgorithm
from src.runners.Runner import Runner


def flatten_fold_aggregate(xs):
    '''
    [[a,b],[c]] -> [(a, 0), (b, 0), (c, 1)]
    
    Inner lists correspond to folds, and folds here corespond to int 
    labels:

    [[edgeA,edgeB],[edgeC]] -> [(edgeA, 0), (edgeB, 0), (edgeC, 1)]
    '''
    flat = [(y, i) for i, ys in enumerate(xs) for y in ys]
    return flat


def flatten_fold_predictions(xs):
    '''
    I need to decompose and re-group these predictions by weight

    1)
    [[{(edge, weight)}]] -> [((edge, weight), fold)]
    [[{a}]] -> [(a, 0)]

    2)
    Regrouping:
    [((edge, weight), fold)] -> [{((edge, weight),fold)}]
    
    3)
    Making the items match the positives/negatives:
    [{(edge, weight),fold}] -> [{(edge, fold)}]
    ''' 
    flat = [(z, i) for i, ys in enumerate(xs) for y in ys for z in y]

    weights = set([x[0][1] for x in flat])
    weights = list(weights)
    weights.sort(reverse=True)

    '''
    regrouped = []

    # Group by weights. This is slow; I can probably speed it up though
    for weight in weights:
        s = {x for x in flat if x[0][1] == weight}
        regrouped.append(s)
    '''
    groups = {}

    for weight in weights:
        groups[weight] = set()

    for x in flat:
        groups[x[0][1]].add(x)

    regrouped = [] 
    for k, v in sorted(groups.items(), reverse=True):
        regrouped.append(v)

    final = [{(x[0][0], x[1]) for x in xs} for xs in regrouped]

    return final


class AlgorithmRunner(Runner):
    '''
    This class serves as a basis for pathway reconstruction runners that:

    1) Creates training and test sets from the nodes or edges in an interactome
    2) Runs algorithms over the interactome using the training data
    3) Evaluates the results of running those algorithms
    4) Produces visualizations of the evalution results

    Directed and undirected edges must be split into folds
    independently. This is because undirected edges are represented
    with a pair of directed edges under our model, and we want
    both directions to either in the same bin (train or test).
    '''

    def __init__(self, *args, **kwargs) -> None:
        '''
        Sets instance variables
        '''
        super().__init__(*args, **kwargs)

        # Cache remaining keyword arguments
        self.options = kwargs
        
        print('initializing folds')
        self.training_folds = self.get_training_folds()
        self.test_folds = self.get_test_folds()
        print('done initializing folds')
    

    def get_name(self) -> str:
        raise NotImplementedError()


    def get_details(self) -> str:
        raise NotImplementedError()


    def get_output_prefix(self):
        raise NotImplementedError()


    def get_fold_creator(self, pathway):
        raise NotImplementedError()


    def get_fold_creators(self):
        # Check if fold creators exist in the cache
        if hasattr(self, "fcs"):
            return self.fcs

        fcs = []

        for pathway in self.collection.pathways:
            fcs.append(self.get_fold_creator(pathway))
        
        # Cache the fold creators because they cache things
        self.fcs = fcs

        return fcs


    def get_filtered_folds(self, folds):
        '''
        Filter the undirected edges of training/test folds
        so that only one direction is present.

        :param folds: a four-tuple:
            - dir_training edges
            - dir_test edges
            - undir_training edges
            - undir_test edges
        '''

        print("# of folds: ", len(folds))
        #dir_map = self.get_dir_map()

        real_folds = []

        for c, fold in enumerate(folds):
            # Filter each fold's undirected train/test edges
            print("filtering undir train") 

            # This is a little hacky but for this AlgorithmEvaluator I can
            # just hobble the filtering and achieve the desired results

            '''
            filtered_training = iu.filter_edge_direction(
                fold[2], dir_map)

            print("filtering undir test")
            filtered_test = iu.filter_edge_direction(
                fold[3], dir_map)
            '''

            filtered_training = fold[2]
            filtered_test = fold[3]

            # Now I need to create a list of true folds from these.
            training = list(set(filtered_training).union(
                set(fold[0])))

            test = list(set(filtered_test).union(set(fold[1])))

            real_folds.append((training, test))

        return real_folds


    def get_training_folds(self):
        '''
        Returns a list of list of tuples (positives, negatives)
        detailing the set of training positives/negatives for each
        fold
        '''
        fold_creators = self.get_fold_creators()
        #folds = []

        pos_list= []
        neg_list = []

        for fc in fold_creators:
            print("-------------------------------------------")
            print("Creating training folds for", fc.pathway.name)
            print("-------------------------------------------")
            pos_folds, neg_folds = fc.create_folds()
            pos_training = [x[0] for x in self.get_filtered_folds(pos_folds)]
            neg_training = [x[0] for x in self.get_filtered_folds(neg_folds)]

            pos_list.append(pos_training)
            neg_list.append(neg_training)

        # return folds
        intermediate = list(zip(pos_list, neg_list))

        return [list(zip(x, y)) for x, y in intermediate]


    def get_test_folds(self):
        '''
        Returns a list of tuples (positives, negatives) detailing the
        set of test positives/negatives for each fold
        '''
        fold_creators = self.get_fold_creators()

        pos_list= []
        neg_list = []

        for fc in fold_creators:
            print("-------------------------------------------")
            print("Creating test folds for", fc.pathway.name)
            print("-------------------------------------------")
            pos_folds, neg_folds = fc.create_folds()
            pos_test = [x[1] for x in self.get_filtered_folds(pos_folds)]
            neg_test = [x[1] for x in self.get_filtered_folds(neg_folds)]

            pos_list.append(pos_test)
            neg_list.append(neg_test)

        intermediate = list(zip(pos_list, neg_list))

        return [list(zip(x, y)) for x, y in intermediate]


    def run_reconstructions(self, output_dir=Path()):
        '''
        Run each algorithm over each fold of each pathway in the 
        pathway collection.
        '''

        # self.training_folds is a list of tuples of lists
        pairs = zip(self.collection.pathways, self.training_folds)

        # training_folds is a tuple of lists corresponding to the pathway
        for i, (pathway, training_fold) in enumerate(pairs):

            # OKAY, so the problem is that training fold is
            # a tuple of lists. It needs to be a list of tuples
            for j, fold in enumerate(training_fold):

                for algorithm in self.algorithms:

                    # First, write output directory
                    full_output_dir = Path(
                        output_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    alg_dir = algorithm.get_full_output_directory(
                        full_output_dir)
                    
                    # TODO: Is this step even necessary?
                    alg_dir.mkdir(parents=True, exist_ok=True)
                    
                    print("# of training positives: ", 
                        len(fold[0]))

                    print("# of training negatives: ", 
                        len(fold[1]))

                    alg_input = RankingAlgorithm.PathwayReconstructionInput(
                        self.interactome.path,
                        fold[0], # training positive edges
                        pathway.get_nodes_file(),
                        full_output_dir,
                        pathway.get_edges_file(),
                        fold[1]) # training negative edges

                    print("Running algorithm: ")
                    print("    " + self.interactome.name)
                    print("    " + self.collection.name)
                    print("    " + pathway.name)
                    print("    " + self.get_name())
                    print("    " + "fold-" + str(j))
                    self.run_alg(algorithm, alg_input)


    def run_alg(self, algorithm, alg_input):
        '''
        Run an algorithm, keeping track of and printing the time that it
        takes to run.
        '''
        print("    Running " + algorithm.get_descriptive_name())
        start = time.time()
        algorithm.run_wrapper(alg_input, should_force=False)
        end = time.time()
        print("    Time to run: " + str(end - start))
        print("-----------------------------------------------------")


    def evaluate_reconstructions(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Calculate performance metrics, like precision/recall scores.
        '''
        
        None
        # AUPRC, Average Precision, and Max F1-measure plots 
        #self.calculate_metrics(reconstruction_dir, evaluation_dir)
        #self.calculate_metrics_nodes(reconstruction_dir, evaluation_dir)
        #self.calculate_wilcoxon_scores(evaluation_dir)
        #self.calculate_wilcoxon_scores_nodes(evaluation_dir)

        ###### Precision/Recall plots
        #self.aggregate_pr_over_folds(reconstruction_dir, evaluation_dir)
        #self.aggregate_pr_over_pathways(evaluation_dir)

        #####self.successively_aggregate_pr_by_auprc(evaluation_dir)
        #####self.successively_aggregate_pr_by_auprc_nodes(evaluation_dir)
        self.successively_aggregate_pr_by_pathway_auprc(evaluation_dir)
        self.successively_aggregate_pr_by_pathway_auprc_nodes(evaluation_dir)
        #self.calculate_wilcoxon_scores_pathway(evaluation_dir)
        #self.calculate_wilcoxon_scores_pathway_nodes(evaluation_dir)

        #### Precision per rank 
        self.calculate_and_plot_precision_per_rank(
            reconstruction_dir, evaluation_dir)


    def calculate_metrics_nodes(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Calculate the precision recall curve and average precison
        for each fold independently, writing it to disk.
        '''

        pairs = zip(self.collection.pathways, self.test_folds)

        # TODO These varible names are AWFUL LIES.
        # dir_neg and undirected return ALL edges in the interactome,
        # including pathway edges.

        # Create an interactome to figure out node negatives
        # Get directed negatives
        dir_neg = iu.get_directed_interactome_edges(
            self.interactome.path,
            self.interactome.direction_file)

        # Get undirected negatives
        undir_neg = iu.get_undirected_interactome_edges(
            self.interactome.path,
            self.interactome.direction_file)

        interactome_a = nx.DiGraph()
        for edge in dir_neg:
            interactome_a.add_edge(edge[0], edge[1])

        for edge in undir_neg:
            interactome_a.add_edge(edge[0], edge[1])
            interactome_a.add_edge(edge[1], edge[0])

        for pathway, test_fold in pairs:

            # Create a pathway object to figure out node positives
            pathway_obj_a = pathway.get_pathway_obj().get_net_from_pathway()

            for algorithm in self.algorithms:
                for j, fold in enumerate(test_fold):
                    
                    # TODO: get a directory from a function instead

                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        reconstruction_file.touch()
        
                    # Where we will write precision/recall results
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))


                    fold_predictions = None
                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
               

                    fold_predictions = [
                        set([tup[0] for tup in s])
                        for s in fold_predictions]




                    '''
                    NEW
                    '''

                    interactome = interactome_a.copy()
                    pathway_obj = pathway_obj_a.copy()

                    # Remove pathway edges from the interactome here!
                    for edge in pathway_obj.edges():
                        interactome.remove_edge(edge[0], edge[1])

                    # We can't delete pathway nodes because it removes edges
                    # from the interactome, and we won't be able to determine
                    # which nodes are isolated because of pathway node deletion
                    # and which are isolated because we chose that node for
                    # deletion.

                    #for node in pathway_obj.nodes():
                    #    interactome.remove_node(node)

                    def expand_predictions(preds):
                        expansion = []
                        for rank in preds:
                            items = set()
                            for item in rank:
                                items.add(item[0])
                                items.add(item[0])

                            expansion.append(items)

                        return expansion
            
                    node_pred = expand_predictions(fold_predictions)
    
                    nodes = interactome.nodes()
                    nodes_with_no_neighbors = []

                    # This should not happen, but it never hurts
                    # to be safe. It cannot happen in an S-T pruned pathway.
                    # Changed my mind. This can totally happen if I delete
                    # pathway edges beforehand.

                    #for node in nodes:
                    #    if len(interactome.neighbors(node)) == 0 and len(interactome.predecessors(node)) == 0:
                    #        print("Hey, I just met you, and this is craaaazy")
                    #        print("but I have no friends; so delete me baby!")
                    #        print(node)
                    #        nodes_with_no_neighbors.append(node)
                    #        interactome.remove_node(node)

                    # Delete the edges (as obtained from the fold)

                    positives = set(fold[0])
                    negatives = set(fold[1])
                    print("positives", len(positives))
                    print(positives)
                    print("negatives", len(negatives))

                    node_test_pos = set()
                    node_test_neg = set()
                    for p in positives:
                        #print("Deleting edge", p, "in fold", j)
                        pathway_obj.remove_edge(p[0], p[1])

                    for n in negatives:
                        interactome.remove_edge(n[0], n[1])
                    lonely_pathway_nodes = []
                    lonely_interactome_nodes = []

                    pathway_nodes = pathway_obj.nodes()
                    interactome_nodes = interactome.nodes()

                    for node in pathway_nodes:
                        if len(pathway_obj.neighbors(node)) == 0 and len(pathway_obj.predecessors(node)) == 0:
                            lonely_pathway_nodes.append(node)

                    for node in interactome_nodes:
                        if len(interactome.neighbors(node)) == 0 and len(interactome.predecessors(node)) == 0:
                            lonely_interactome_nodes.append(node) 

                    lonely_interactome_nodes = set(lonely_interactome_nodes)

                    # Nodes without any edges that are not pathway edges are
                    # the nodes that we marked for deletion, as no node in
                    # the interactome is "lonely" to start
                    lonely_interactome_nodes = \
                        lonely_interactome_nodes - \
                            set(pathway_obj_a.nodes()) 

                    for node in lonely_pathway_nodes:
                        node_test_pos.add(node)

                    for node in lonely_interactome_nodes:
                        node_test_neg.add(node)

                    print(pathway.name, j)
                    print("computing p/r for nodes")
                    print("len node test pos", len(node_test_pos))
                    print("len node test neg", len(node_test_neg))

                    # Call existing precrec functions passing these things above
                    points = \
                        precrec.compute_precision_recall_curve_negatives_decimals(
                            node_pred, node_test_pos, node_test_neg)

                    print("done")

                    '''
                    END NEW
                    '''







                    # Calculate all of the measures and write them out
                    avg_prec = precrec.compute_average_precision(points)

                    f1_max = precrec.compute_f_max(points)

                    # TODO: Rename utils function to compute_AUPRC
                    auprc = precrec.compute_AUPRC(points)

                    avg_prec_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "average_precision-nodes.txt") 

                    f1_max_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "f1_max-nodes.txt") 

                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc-nodes.txt") 

                    avg_prec_file.parent.mkdir(parents=True, exist_ok=True)

                    with avg_prec_file.open('w') as f: 
                        f.write(str(avg_prec))

                    with f1_max_file.open('w') as f:
                        f.write(str(f1_max))

                    with auprc_file.open('w') as f:
                        f.write(str(auprc))


    def calculate_metrics(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Calculate the precision recall curve and average precison
        for each fold independently, writing it to disk.
        '''

        pairs = zip(self.collection.pathways, self.test_folds)

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                for j, fold in enumerate(test_fold):
                    
                    # TODO: get a directory from a function instead

                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        reconstruction_file.touch()
        
                    # Where we will write precision/recall results
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    positives = set(fold[0])
                    negatives = set(fold[1])

                    fold_predictions = None
                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
               
                    fold_predictions = [
                        set([tup[0] for tup in s])
                        for s in fold_predictions]

                    print("Calculating precision/recall points")
                    t1 = time.time()
                    points = \
                        precrec.compute_precision_recall_curve_negatives_decimals(fold_predictions, positives, negatives)
                    t2 = time.time()
                    print("Done! That took: ", t2 - t1)

                    # Calculate all of the measures and write them out
                    avg_prec = precrec.compute_average_precision(points)

                    f1_max = precrec.compute_f_max(points)

                    # TODO: Rename utils function to compute_AUPRC
                    auprc = precrec.compute_AUPRC(points)

                    avg_prec_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "average_precision.txt") 

                    f1_max_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "f1_max.txt") 

                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc.txt") 

                    avg_prec_file.parent.mkdir(parents=True, exist_ok=True)

                    with avg_prec_file.open('w') as f: 
                        f.write(str(avg_prec))

                    with f1_max_file.open('w') as f:
                        f.write(str(f1_max))

                    with auprc_file.open('w') as f:
                        f.write(str(auprc))



    def calculate_wilcoxon_scores_pathway(self, evaluation_dir=Path()):
   
        '''
        Write out a TSV relating Wilcoxon scores between algorithms


        Uses the 15 points from each pathway as opposed to the 150 points
        from each pathway fold.
        '''

        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            auprc_algorithm_map[name] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:

                eval_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pr_file = Path(
                    eval_dir,
                    algorithm.get_output_directory(),
                    "precision-recall.txt")

                with pr_file.open('r') as f:
                    pr_points = precrec.read_precision_recall_fractions(f)

                recall = [float(p[0]) / float(p[2]) for p in pr_points]
                precision = [float(p[0]) / float(p[1]) if p[1] != 0 else 1 for p in pr_points]

                decimals = list(zip(precision, recall))

                # Then get the AUPRC
                auprc = precrec.compute_AUPRC(decimals)

                name = algorithm.get_descriptive_name()

                auprc_algorithm_map[name].append(auprc)

        # Second, calculate the Wilcoxon stat for each measure
        alpha = .05

        correction = sp.special.comb(len(self.algorithms), 2)

        # This is the Bonferroni-corrected alpha value
        corrected_alpha = alpha / correction

        wilcoxon_auprc = []

        for i, alg_a in enumerate(self.algorithms):
            wilcoxon_auprc.append([])

            for j, alg_b in enumerate(self.algorithms):
                name_a = alg_a.get_descriptive_name()
                name_b = alg_b.get_descriptive_name()

                if name_a == name_b:
                    wilcoxon_auprc[i].append((-1, -1, -1))
                    continue

                auprc_list_a = auprc_algorithm_map[name_a]

                auprc_list_b = auprc_algorithm_map[name_b]

                print(alg_a.get_descriptive_name())
                print(alg_b.get_descriptive_name())

                print(len(auprc_list_a))
                print(len(auprc_list_b))


                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html

                auprc_stat, auprc_p_val = \
                    sp.stats.wilcoxon(auprc_list_a, auprc_list_b)
                print("wilcox3")

                auprc_a_median = np.median(auprc_list_a)
                auprc_b_median = np.median(auprc_list_b)

                wilcoxon_auprc[i].append(
                    (auprc_a_median, auprc_b_median, auprc_p_val))


        # Third, write these scores out to disk
        wilcoxon_auprc_file = Path(
            evaluation_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "wilcoxon-auprc.tsv")

        wilcoxon_auprc_file.parent.mkdir(parents=True, exist_ok=True)

        with wilcoxon_auprc_file.open('w') as f: 
            # Write header
            f.write("Alg1\tAlg2\tMedian1\tMedian2\tp-val (Wilcoxon)\n")

            # Write out matrix2 to the tsv, row by row
            for i, algorithm in enumerate(self.algorithms):

                for j, algorithm2 in enumerate(self.algorithms):
                    name1 = algorithm.get_descriptive_name()
                    name2 = algorithm2.get_descriptive_name()
                
                    median1 = wilcoxon_auprc[i][j][0]
                    median2 = wilcoxon_auprc[i][j][1]
                    p_val = wilcoxon_auprc[i][j][2]

                    f.write("%s\t%s\t%.10f\t%.10f\t%.10f\n" % (
                        name1,
                        name2,
                        median1,
                        median2,
                        p_val))




    def calculate_wilcoxon_scores_pathway_nodes(self, evaluation_dir=Path()):
   
        '''
        Write out a TSV relating Wilcoxon scores between algorithms


        Uses the 15 points from each pathway as opposed to the 150 points
        from each pathway fold.
        '''

        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            auprc_algorithm_map[name] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:

                eval_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pr_file = Path(
                    eval_dir,
                    algorithm.get_output_directory(),
                    "precision-recall-nodes.txt")

                with pr_file.open('r') as f:
                    pr_points = precrec.read_precision_recall_fractions(f)

                recall = [float(p[0]) / float(p[2]) for p in pr_points]
                precision = [float(p[0]) / float(p[1]) if p[1] != 0 else 1 for p in pr_points]

                decimals = list(zip(precision, recall))

                # Then get the AUPRC
                auprc = precrec.compute_AUPRC(decimals)

                name = algorithm.get_descriptive_name()

                auprc_algorithm_map[name].append(auprc)

        # Second, calculate the Wilcoxon stat for each measure
        alpha = .05

        correction = sp.special.comb(len(self.algorithms), 2)

        # This is the Bonferroni-corrected alpha value
        corrected_alpha = alpha / correction

        wilcoxon_auprc = []

        for i, alg_a in enumerate(self.algorithms):
            wilcoxon_auprc.append([])

            for j, alg_b in enumerate(self.algorithms):
                name_a = alg_a.get_descriptive_name()
                name_b = alg_b.get_descriptive_name()

                if name_a == name_b:
                    wilcoxon_auprc[i].append((-1, -1, -1))
                    continue

                auprc_list_a = auprc_algorithm_map[name_a]

                auprc_list_b = auprc_algorithm_map[name_b]

                print(alg_a.get_descriptive_name())
                print(alg_b.get_descriptive_name())

                print(len(auprc_list_a))
                print(len(auprc_list_b))


                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html

                auprc_stat, auprc_p_val = \
                    sp.stats.wilcoxon(auprc_list_a, auprc_list_b)
                print("wilcox3")

                auprc_a_median = np.median(auprc_list_a)
                auprc_b_median = np.median(auprc_list_b)

                wilcoxon_auprc[i].append(
                    (auprc_a_median, auprc_b_median, auprc_p_val))


        # Third, write these scores out to disk
        wilcoxon_auprc_file = Path(
            evaluation_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "wilcoxon-auprc-nodes.tsv")

        wilcoxon_auprc_file.parent.mkdir(parents=True, exist_ok=True)

        with wilcoxon_auprc_file.open('w') as f: 
            # Write header
            f.write("Alg1\tAlg2\tMedian1\tMedian2\tp-val (Wilcoxon)\n")

            # Write out matrix2 to the tsv, row by row
            for i, algorithm in enumerate(self.algorithms):

                for j, algorithm2 in enumerate(self.algorithms):
                    name1 = algorithm.get_descriptive_name()
                    name2 = algorithm2.get_descriptive_name()
                
                    median1 = wilcoxon_auprc[i][j][0]
                    median2 = wilcoxon_auprc[i][j][1]
                    p_val = wilcoxon_auprc[i][j][2]

                    f.write("%s\t%s\t%.10f\t%.10f\t%.10f\n" % (
                        name1,
                        name2,
                        median1,
                        median2,
                        p_val))



    def calculate_wilcoxon_scores(self, evaluation_dir=Path()):
        '''
        Write out a TSV relating Wilcoxon scores between algorithms
        '''
        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        avg_prec_algorithm_map = {}
        f1_max_algorithm_map = {}
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            avg_prec_algorithm_map[name] = []
            f1_max_algorithm_map[name] = []
            auprc_algorithm_map[name] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []
                for j, fold in enumerate(test_fold):
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    avg_prec_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "average_precision.txt") 

                    f1_max_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "f1_max.txt") 

                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc.txt") 

                    avg_prec = None
                    f1_max = None
                    auprc = None

                    with avg_prec_file.open('r') as f:
                        line = next(f)
                        avg_prec = float(line.strip())
                        
                    with f1_max_file.open('r') as f:
                        line = next(f)
                        f1_max = float(line.strip())

                    with auprc_file.open('r') as f:
                        line = next(f)
                        auprc = float(line.strip())

                    name = algorithm.get_descriptive_name()

                    avg_prec_algorithm_map[name].append(avg_prec) 
                    f1_max_algorithm_map[name].append(f1_max)
                    auprc_algorithm_map[name].append(auprc)

        # Second, calculate the Wilcoxon stat for each measure
        alpha = .05

        correction = sp.special.comb(len(self.algorithms), 2)

        # This is the Bonferroni-corrected alpha value
        corrected_alpha = alpha / correction

        wilcoxon_avg_prec = []
        wilcoxon_f1_max = []
        wilcoxon_auprc = []

        for i, alg_a in enumerate(self.algorithms):
            wilcoxon_avg_prec.append([])
            wilcoxon_f1_max.append([])
            wilcoxon_auprc.append([])

            for j, alg_b in enumerate(self.algorithms):
                name_a = alg_a.get_descriptive_name()
                name_b = alg_b.get_descriptive_name()

                if name_a == name_b:
                    wilcoxon_avg_prec[i].append((-1, -1, -1))
                    wilcoxon_f1_max[i].append((-1, -1, -1))
                    wilcoxon_auprc[i].append((-1, -1, -1))
                    continue

                avg_prec_list_a = avg_prec_algorithm_map[name_a]
                f1_max_list_a = f1_max_algorithm_map[name_a]
                auprc_list_a = auprc_algorithm_map[name_a]

                avg_prec_list_b = avg_prec_algorithm_map[name_b]
                f1_max_list_b = f1_max_algorithm_map[name_b]
                auprc_list_b = auprc_algorithm_map[name_b]

                print(alg_a.get_descriptive_name())
                print(alg_b.get_descriptive_name())
                print(len(avg_prec_list_a))
                print(len(avg_prec_list_b))

                print(len(f1_max_list_a))
                print(len(f1_max_list_b))

                print(len(auprc_list_a))
                print(len(auprc_list_b))


                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
                avg_prec_stat, avg_prec_p_val = \
                    sp.stats.wilcoxon(avg_prec_list_a, avg_prec_list_b)
                print("wilcox1")

                f1_max_stat, f1_max_p_val = \
                    sp.stats.wilcoxon(f1_max_list_a, f1_max_list_b)
                print("wilcox2")

                auprc_stat, auprc_p_val = \
                    sp.stats.wilcoxon(auprc_list_a, auprc_list_b)
                print("wilcox3")

                avg_prec_a_median = np.median(avg_prec_list_a)
                avg_prec_b_median = np.median(avg_prec_list_b)

                f1_max_a_median = np.median(f1_max_list_a)
                f1_max_b_median = np.median(f1_max_list_b)

                auprc_a_median = np.median(auprc_list_a)
                auprc_b_median = np.median(auprc_list_b)

                wilcoxon_avg_prec[i].append(
                    (avg_prec_a_median, avg_prec_b_median, avg_prec_p_val))

                wilcoxon_f1_max[i].append(
                    (f1_max_a_median, f1_max_b_median, f1_max_p_val))

                wilcoxon_auprc[i].append(
                    (auprc_a_median, auprc_b_median, auprc_p_val))


        # Third, write these scores out to disk
        # TODO: Yes, I know that all of this repetition is screaming for
        # a function. Don't let the perfect be the enemy of the good.

        wilcoxon_avg_prec_file = Path(
            evaluation_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "wilcoxon-avg-prec.tsv")

        wilcoxon_f1_max_file = Path(
            evaluation_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "wilcoxon-f1-max.tsv")

        wilcoxon_auprc_file = Path(
            evaluation_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "wilcoxon-auprc.tsv")

        wilcoxon_avg_prec_file.parent.mkdir(parents=True, exist_ok=True)

        with wilcoxon_avg_prec_file.open('w') as f: 
            # Write header
            f.write("Alg1\tAlg2\tMedian1\tMedian2\tp-val (Wilcoxon)\n")

            # Write out matrix2 to the tsv, row by row
            for i, algorithm in enumerate(self.algorithms):

                for j, algorithm2 in enumerate(self.algorithms):
                    name1 = algorithm.get_descriptive_name()
                    name2 = algorithm2.get_descriptive_name()
                
                    median1 = wilcoxon_avg_prec[i][j][0]
                    median2 = wilcoxon_avg_prec[i][j][1]
                    p_val = wilcoxon_avg_prec[i][j][2]

                    f.write("%s\t%s\t%.10f\t%.10f\t%.10f\n" % (
                        name1,
                        name2,
                        median1,
                        median2,
                        p_val))

        with wilcoxon_f1_max_file.open('w') as f: 
            # Write header
            f.write("Alg1\tAlg2\tMedian1\tMedian2\tp-val (Wilcoxon)\n")

            # Write out matrix2 to the tsv, row by row
            for i, algorithm in enumerate(self.algorithms):

                for j, algorithm2 in enumerate(self.algorithms):
                    name1 = algorithm.get_descriptive_name()
                    name2 = algorithm2.get_descriptive_name()
                
                    median1 = wilcoxon_f1_max[i][j][0]
                    median2 = wilcoxon_f1_max[i][j][1]
                    p_val = wilcoxon_f1_max[i][j][2]

                    f.write("%s\t%s\t%.10f\t%.10f\t%.10f\n" % (
                        name1,
                        name2,
                        median1,
                        median2,
                        p_val))

        with wilcoxon_auprc_file.open('w') as f: 
            # Write header
            f.write("Alg1\tAlg2\tMedian1\tMedian2\tp-val (Wilcoxon)\n")

            # Write out matrix2 to the tsv, row by row
            for i, algorithm in enumerate(self.algorithms):

                for j, algorithm2 in enumerate(self.algorithms):
                    name1 = algorithm.get_descriptive_name()
                    name2 = algorithm2.get_descriptive_name()
                
                    median1 = wilcoxon_auprc[i][j][0]
                    median2 = wilcoxon_auprc[i][j][1]
                    p_val = wilcoxon_auprc[i][j][2]

                    f.write("%s\t%s\t%.10f\t%.10f\t%.10f\n" % (
                        name1,
                        name2,
                        median1,
                        median2,
                        p_val))



    def calculate_wilcoxon_scores_nodes(self, evaluation_dir=Path()):
        '''
        Write out a TSV relating Wilcoxon scores between algorithms
        '''
        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        avg_prec_algorithm_map = {}
        f1_max_algorithm_map = {}
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            avg_prec_algorithm_map[name] = []
            f1_max_algorithm_map[name] = []
            auprc_algorithm_map[name] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []
                for j, fold in enumerate(test_fold):
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    avg_prec_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "average_precision-nodes.txt") 

                    f1_max_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "f1_max-nodes.txt") 

                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc-nodes.txt") 

                    avg_prec = None
                    f1_max = None
                    auprc = None

                    with avg_prec_file.open('r') as f:
                        line = next(f)
                        avg_prec = float(line.strip())
                        
                    with f1_max_file.open('r') as f:
                        line = next(f)
                        f1_max = float(line.strip())

                    with auprc_file.open('r') as f:
                        line = next(f)
                        auprc = float(line.strip())

                    name = algorithm.get_descriptive_name()

                    avg_prec_algorithm_map[name].append(avg_prec) 
                    f1_max_algorithm_map[name].append(f1_max)
                    auprc_algorithm_map[name].append(auprc)

        # Second, calculate the Wilcoxon stat for each measure
        alpha = .05

        correction = sp.special.comb(len(self.algorithms), 2)

        # This is the Bonferroni-corrected alpha value
        corrected_alpha = alpha / correction

        wilcoxon_avg_prec = []
        wilcoxon_f1_max = []
        wilcoxon_auprc = []

        for i, alg_a in enumerate(self.algorithms):
            wilcoxon_avg_prec.append([])
            wilcoxon_f1_max.append([])
            wilcoxon_auprc.append([])

            for j, alg_b in enumerate(self.algorithms):
                name_a = alg_a.get_descriptive_name()
                name_b = alg_b.get_descriptive_name()

                if name_a == name_b:
                    wilcoxon_avg_prec[i].append((-1, -1, -1))
                    wilcoxon_f1_max[i].append((-1, -1, -1))
                    wilcoxon_auprc[i].append((-1, -1, -1))
                    continue

                avg_prec_list_a = avg_prec_algorithm_map[name_a]
                f1_max_list_a = f1_max_algorithm_map[name_a]
                auprc_list_a = auprc_algorithm_map[name_a]

                avg_prec_list_b = avg_prec_algorithm_map[name_b]
                f1_max_list_b = f1_max_algorithm_map[name_b]
                auprc_list_b = auprc_algorithm_map[name_b]

                print(alg_a.get_descriptive_name())
                print(alg_b.get_descriptive_name())
                print(len(avg_prec_list_a))
                print(len(avg_prec_list_b))

                print(len(f1_max_list_a))
                print(len(f1_max_list_b))

                print(len(auprc_list_a))
                print(len(auprc_list_b))


                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
                avg_prec_stat, avg_prec_p_val = \
                    sp.stats.wilcoxon(avg_prec_list_a, avg_prec_list_b)
                print("wilcox1")

                f1_max_stat, f1_max_p_val = \
                    sp.stats.wilcoxon(f1_max_list_a, f1_max_list_b)
                print("wilcox2")

                auprc_stat, auprc_p_val = \
                    sp.stats.wilcoxon(auprc_list_a, auprc_list_b)
                print("wilcox3")

                avg_prec_a_median = np.median(avg_prec_list_a)
                avg_prec_b_median = np.median(avg_prec_list_b)

                f1_max_a_median = np.median(f1_max_list_a)
                f1_max_b_median = np.median(f1_max_list_b)

                auprc_a_median = np.median(auprc_list_a)
                auprc_b_median = np.median(auprc_list_b)

                wilcoxon_avg_prec[i].append(
                    (avg_prec_a_median, avg_prec_b_median, avg_prec_p_val))

                wilcoxon_f1_max[i].append(
                    (f1_max_a_median, f1_max_b_median, f1_max_p_val))

                wilcoxon_auprc[i].append(
                    (auprc_a_median, auprc_b_median, auprc_p_val))


        # Third, write these scores out to disk
        # TODO: Yes, I know that all of this repetition is screaming for
        # a function. Don't let the perfect be the enemy of the good.

        wilcoxon_avg_prec_file = Path(
            evaluation_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "wilcoxon-avg-prec-nodes.tsv")

        wilcoxon_f1_max_file = Path(
            evaluation_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "wilcoxon-f1-max-nodes.tsv")

        wilcoxon_auprc_file = Path(
            evaluation_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "wilcoxon-auprc-nodes.tsv")

        wilcoxon_avg_prec_file.parent.mkdir(parents=True, exist_ok=True)

        with wilcoxon_avg_prec_file.open('w') as f: 
            # Write header
            f.write("Alg1\tAlg2\tMedian1\tMedian2\tp-val (Wilcoxon)\n")

            # Write out matrix2 to the tsv, row by row
            for i, algorithm in enumerate(self.algorithms):

                for j, algorithm2 in enumerate(self.algorithms):
                    name1 = algorithm.get_descriptive_name()
                    name2 = algorithm2.get_descriptive_name()
                
                    median1 = wilcoxon_avg_prec[i][j][0]
                    median2 = wilcoxon_avg_prec[i][j][1]
                    p_val = wilcoxon_avg_prec[i][j][2]

                    f.write("%s\t%s\t%.10f\t%.10f\t%.10f\n" % (
                        name1,
                        name2,
                        median1,
                        median2,
                        p_val))

        with wilcoxon_f1_max_file.open('w') as f: 
            # Write header
            f.write("Alg1\tAlg2\tMedian1\tMedian2\tp-val (Wilcoxon)\n")

            # Write out matrix2 to the tsv, row by row
            for i, algorithm in enumerate(self.algorithms):

                for j, algorithm2 in enumerate(self.algorithms):
                    name1 = algorithm.get_descriptive_name()
                    name2 = algorithm2.get_descriptive_name()
                
                    median1 = wilcoxon_f1_max[i][j][0]
                    median2 = wilcoxon_f1_max[i][j][1]
                    p_val = wilcoxon_f1_max[i][j][2]

                    f.write("%s\t%s\t%.10f\t%.10f\t%.10f\n" % (
                        name1,
                        name2,
                        median1,
                        median2,
                        p_val))

        with wilcoxon_auprc_file.open('w') as f: 
            # Write header
            f.write("Alg1\tAlg2\tMedian1\tMedian2\tp-val (Wilcoxon)\n")

            # Write out matrix2 to the tsv, row by row
            for i, algorithm in enumerate(self.algorithms):

                for j, algorithm2 in enumerate(self.algorithms):
                    name1 = algorithm.get_descriptive_name()
                    name2 = algorithm2.get_descriptive_name()
                
                    median1 = wilcoxon_auprc[i][j][0]
                    median2 = wilcoxon_auprc[i][j][1]
                    p_val = wilcoxon_auprc[i][j][2]

                    f.write("%s\t%s\t%.10f\t%.10f\t%.10f\n" % (
                        name1,
                        name2,
                        median1,
                        median2,
                        p_val))









    def aggregate_pr_over_folds(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Merge the precision/recall results per pathway using folds
        '''

        print("----------------------------------------------------")
        print("Computing Precision/Recall per Fold and Aggregating")
        print("----------------------------------------------------")

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.collection.pathways],
            fold_creators)

        for i, (pathway, fc) in enumerate(creator_pathway_pairs):
            print(pathway.name)
            
            test_folds = self.test_folds[i]

            print(len(test_folds))

            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "aggregate")

            for algorithm in self.algorithms:
                print("    " + algorithm.get_descriptive_name())
                predictions = []
                test_positives = []
                test_negatives = []

                avg_prec = []

                for j, fold in enumerate(test_folds):
                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))
                    
                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        print("ALERT: RECONSTRUCTION FILE NOT FOUND")
                        reconstruction_file.touch()

                    positives = fold[0]
                    negatives = fold[1]

                    print("\nlength of positives", len(positives))
                    print("length of negtives", len(negatives))

                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
                        predictions.append(fold_predictions)

                    test_positives.append(positives)
                    test_negatives.append(negatives)
                
                print("flattening")
                flat_test_pos = set(flatten_fold_aggregate(test_positives))
                flat_test_neg = set(flatten_fold_aggregate(test_negatives))
                flat_pred = flatten_fold_predictions(predictions)
                print("done")


                '''
                BEGIN NEW STUFF
                '''

                # All I need to do is break the edge predictions up into node
                # predictions

                # And break the positives/negatives up in the same way
                 
                # Then send that stuff to prec/rec calculations  

                # Well... no. It depends on what counts as a negative
                # and a positive for nodes.

                # I would figure that the only nodes that count are
                # those that are completely disconnected from the pathway
                # or the interactome

                # That means I'd have to do it per-fold, because each
                # fold will be missing different things.

                def expand_predictions(preds):
                    expansion = []
                    for rank in preds:
                        items = set()
                        for item in rank:
                            items.add((item[0][0], item[1]))
                            items.add((item[0][1], item[1]))

                        expansion.append(items)

                    return expansion
                    
                node_pred = expand_predictions(flat_pred)

                # Read in the pathway and create 10 copies

                pathway_obj = pathway.get_pathway_obj().get_net_from_pathway()

                pathway_objs = []
                for i in range(len(test_folds)):
                    pathway_objs.append(pathway_obj.copy())

                # Read in the interactome and create 10 copies

                # TODO These varible names are AWFUL.
                # dir_neg and undirected return ALL edges in the interactome,
                # including pathway edges.

                # Get directed negatives
                dir_neg = iu.get_directed_interactome_edges(
                    self.interactome.path,
                    self.interactome.direction_file)

                # Get undirected negatives
                undir_neg = iu.get_undirected_interactome_edges(
                    self.interactome.path,
                    self.interactome.direction_file)

                interactome = nx.DiGraph()
                for edge in dir_neg:
                    interactome.add_edge(edge[0], edge[1])

                for edge in undir_neg:
                    interactome.add_edge(edge[0], edge[1])
                    interactome.add_edge(edge[1], edge[0])

                # Remove pathway edges from the interactome here
                for edge in pathway_obj.edges():
                    interactome.remove_edge(edge[0], edge[1])

                # We can't delete pathway nodes because it removes edges from
                # the interactome, and we won't be able to determine which
                # nodes are isolated because of pathway node deletion and which
                # are isolated because we chose that node for deletion.

                #for node in pathway_obj.nodes():
                #    interactome.remove_node(node)

                nodes = interactome.nodes()
                nodes_with_no_neighbors = []

                # I'm not sure this can happen, but it never hurts
                # to be safe. It cannot happen in an S-T pruned pathway.
                for node in nodes:
                    if len(interactome.neighbors(node)) == 0 and len(interactome.predecessors(node)) == 0:
                        print("Hey, I just met you, and this is craaaazy")
                        print("but I have no friends; so delete me baby!")
                        print(node)
                        nodes_with_no_neighbors.append(node)
                        interactome.remove_node(node)

                interactomes = []

                for i in range(len(test_folds)):
                    interactomes.append(interactome.copy())

                # Delete the edges (as obtained from the fold)

                node_test_pos = set()
                node_test_neg = set()

                for j, fold in enumerate(test_folds):
                    positives = fold[0]
                    negatives = fold[1]
                    
                    for p in positives:
                        #print("Deleting edge", p, "in fold", j)
                        pathway_objs[j].remove_edge(p[0], p[1])

                    nodes = set(interactome.nodes())
                    for n in negatives:
                        if n[0] in nodes and n[1] in nodes:
                            interactomes[j].remove_edge(n[0], n[1])

                    lonely_pathway_nodes = []
                    lonely_interactome_nodes = []

                    pathway_nodes = pathway_objs[j].nodes()
                    interactome_nodes = interactomes[j].nodes()

                    for node in pathway_nodes:
                        if len(pathway_objs[j].neighbors(node)) == 0 and len(pathway_objs[j].predecessors(node)) == 0:
                            lonely_pathway_nodes.append(node)

                    for node in interactome_nodes:
                        if len(interactomes[j].neighbors(node)) == 0 and len(interactomes[j].predecessors(node)) == 0:
                            lonely_interactome_nodes.append(node) 

                    lonely_interactome_nodes = set(lonely_interactome_nodes)

                    # Nodes without any edges that are not pathway edges are
                    # the nodes that we marked for deletion, as no node in
                    # the interactome is "lonely" to start
                    lonely_interactome_nodes = \
                        lonely_interactome_nodes - \
                            set(pathway_obj.nodes(data=False)) 

                    pathway_temp = pathway.get_pathway_obj()
                    sources = set(pathway_temp.get_receptors(data=False))
                    targets = set(pathway_temp.get_tfs(data=False))

                    lonely_pathway_nodes = set(lonely_pathway_nodes)
                    lonely_pathway_nodes = \
                        lonely_pathway_nodes - sources.union(targets)

                    for node in lonely_pathway_nodes:
                        node_test_pos.add((node, j))

                    for node in lonely_interactome_nodes:
                        node_test_neg.add((node, j))

                print("computing p/r for nodes")
                # Call existing precrec functions passing these things above
                points = \
                    precrec.compute_precision_recall_curve_negatives_fractions(
                        node_pred, node_test_pos, node_test_neg)
                print("done")

                new_outfile = Path(
                    pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall-nodes.txt") 

                new_outfile.parent.mkdir(parents=True, exist_ok=True)

                with new_outfile.open("w") as f: 
                    precrec.write_precision_recall_fractions(f, points)
                
                '''
                END NEW STUFF
                '''


                print("computing p/r for edges")
                points = \
                    precrec.compute_precision_recall_curve_negatives_fractions(
                        flat_pred, flat_test_pos, flat_test_neg)
                print("done")

                new_outfile = Path(
                    pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 


                with new_outfile.open("w") as f: 
                    precrec.write_precision_recall_fractions(f, points)


    def aggregate_pr_over_pathways(self, evaluation_dir=Path()):
        '''
        Per algorithm, aggregate the precision/recall scores across
        pathways.
        '''

        print("----------------------------------------------------")
        print("Aggregating Precision/Recall Over Pathways")
        print("----------------------------------------------------")

        # Where we will write precision/recall, aggregated over
        # all pathways
        pathway_collection_pr_output_dir = Path(
            evaluation_dir,
            self.interactome.name,
            self.collection.name,
            "aggregate",
            self.get_output_prefix())

        for algorithm in self.algorithms:    
            edge_curves = []
            node_curves = []


            # Where we wrote precision/recall, aggregated over
            # all folds per pathway
            for pathway in self.collection.pathways:
                pathway_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pathway_pr_outfile = Path(
                    pathway_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                pathway_pr_outfile_nodes = Path(
                    pathway_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall-nodes.txt") 

                with pathway_pr_outfile.open('r') as f:
                    curve = precrec.read_precision_recall_fractions(f)
                    edge_curves.append(curve)

                with pathway_pr_outfile_nodes.open('r') as f:
                    curve = precrec.read_precision_recall_fractions(f)
                    node_curves.append(curve)

            print("Aggregating precision/recall curves")

            aggregated_edges = precrec.aggregate_precision_recall_curve_fractions_stretch(
                edge_curves)

            aggregated_nodes = precrec.aggregate_precision_recall_curve_fractions_stretch(
                node_curves)

            # Write aggregated curve back out
            pathway_collection_pr_outfile = Path(
                pathway_collection_pr_output_dir, 
                algorithm.get_output_directory(),
                "precision-recall.txt") 

            # Write aggregated curve back out
            pathway_collection_pr_outfile_nodes = Path(
                pathway_collection_pr_output_dir, 
                algorithm.get_output_directory(),
                "precision-recall-nodes.txt") 

            pathway_collection_pr_outfile.parent.mkdir(
                parents=True, exist_ok=True)

            with pathway_collection_pr_outfile.open("w") as f: 
                precrec.write_precision_recall_fractions(f, aggregated_edges)

            with pathway_collection_pr_outfile_nodes.open("w") as f: 
                precrec.write_precision_recall_fractions(f, aggregated_nodes)

        

    def plot_results(
            self, evaluation_dir=Path(), visualization_dir=Path()):
        '''
        Run all plotting algorithms
        '''

        # Edge Boxplots
        #self.all_algorithm_scores_combined_pathways_boxplot(
        #    evaluation_dir, visualization_dir)

        #self.all_algorithm_scores_individual_pathways_boxplots(
        #    evaluation_dir, visualization_dir)

        #self.individual_algorithm_scores_all_pathways_boxplots(
        #    evaluation_dir, visualization_dir)

        self.all_algorithm_scores_combined_pathways_boxplot_aggregated(
            evaluation_dir, visualization_dir)

        ## Node Boxplots
        #self.all_algorithm_scores_combined_pathways_boxplot_nodes(
        #    evaluation_dir, visualization_dir)

        #self.all_algorithm_scores_individual_pathways_boxplots_nodes(
        #    evaluation_dir, visualization_dir)

        #self.individual_algorithm_scores_all_pathways_boxplots_nodes(
        #    evaluation_dir, visualization_dir)

        self.all_algorithm_scores_combined_pathways_boxplot_aggregated_nodes(
            evaluation_dir, visualization_dir)

        # Precision/Recall
        self.plot_pr_individual_pathways(evaluation_dir, visualization_dir)
        self.plot_pr_all_pathways(evaluation_dir, visualization_dir)

        self.plot_pr_individual_pathways_nodes(
            evaluation_dir, visualization_dir)

        self.plot_pr_all_pathways_nodes(
            evaluation_dir, visualization_dir)


    def all_algorithm_scores_combined_pathways_boxplot_aggregated(
            self, evaluation_dir, visualization_dir): 
        '''
        Make a single boxplot figure, with one boxplot per algorithm, where
        a constituent point is an algorithm's score on a particular fold.
        '''

        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            auprc_algorithm_map[name] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                
                eval_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pr_file = Path(
                    eval_dir,
                    algorithm.get_output_directory(),
                    "precision-recall.txt")

                with pr_file.open('r') as f:
                    pr_points = precrec.read_precision_recall_fractions(f)

                
                # AGHGHHGH I have to convert these to DECIMALS TO USE MY
                # AUPRC FUNCTION

                recall = [float(p[0]) / float(p[2]) for p in pr_points]
                precision = [float(p[0]) / float(p[1]) if p[1] != 0 else 1 for p in pr_points]

                decimals = list(zip(precision, recall))

                # Then get the AUPRC
                auprc = precrec.compute_AUPRC(decimals)

                # Then append it as appropriate
                
                alg_name = algorithm.get_descriptive_name()
                auprc_algorithm_map[alg_name].append(auprc)


        labels = []
        results_auprc = []

        for alg in self.algorithms:
            name = alg.get_descriptive_name()
            labels.append(name)
            results_auprc.append(auprc_algorithm_map[name])

        fig_auprc, ax_auprc = precrec.init_precision_recall_figure()

        ax_auprc.set_title(
            "AUPRC by Algorithm "
            + self.interactome.name + " "
            + self.collection.name + "\n"
            + self.get_details())

        ax_auprc.set_xlabel("AUPRC")
        ax_auprc.set_ylabel("Algorithm")

        auprc_png = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "auprc-aggregated-take-two.png")

        auprc_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "auprc-aggregated-take-two.pdf")

        auprc_png.parent.mkdir(parents=True, exist_ok=True)

        ax_auprc.boxplot(results_auprc, labels=labels, vert=False)

        fig_auprc.savefig(str(auprc_png), bbox_inches='tight')
        fig_auprc.savefig(str(auprc_pdf), bbox_inches='tight')







    def all_algorithm_scores_combined_pathways_boxplot_aggregated_nodes(
            self, evaluation_dir, visualization_dir): 
        '''
        Make a single boxplot figure, with one boxplot per algorithm, where
        a constituent point is an algorithm's score on a particular fold.
        '''

        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            auprc_algorithm_map[name] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                
                eval_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pr_file = Path(
                    eval_dir,
                    algorithm.get_output_directory(),
                    "precision-recall-nodes.txt")

                with pr_file.open('r') as f:
                    pr_points = precrec.read_precision_recall_fractions(f)

                recall = [float(p[0]) / float(p[2]) for p in pr_points]
                precision = [float(p[0]) / float(p[1]) if p[1] != 0 else 1 for p in pr_points]

                decimals = list(zip(precision, recall))


                # Then get the AUPRC
                auprc = precrec.compute_AUPRC(decimals)

                # Then append it as appropriate
                
                alg_name = algorithm.get_descriptive_name()
                auprc_algorithm_map[alg_name].append(auprc)


        labels = []
        results_auprc = []

        for alg in self.algorithms:
            name = alg.get_descriptive_name()
            labels.append(name)
            results_auprc.append(auprc_algorithm_map[name])

        fig_auprc, ax_auprc = precrec.init_precision_recall_figure()

        ax_auprc.set_title(
            "AUPRC by Algorithm "
            + self.interactome.name + " "
            + self.collection.name + "\n"
            + self.get_details())

        ax_auprc.set_xlabel("AUPRC")
        ax_auprc.set_ylabel("Algorithm")

        auprc_png = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "auprc-aggregated-take-two-nodes.png")

        auprc_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "auprc-aggregated-take-two-nodes.pdf")

        auprc_png.parent.mkdir(parents=True, exist_ok=True)

        ax_auprc.boxplot(results_auprc, labels=labels, vert=False)

        fig_auprc.savefig(str(auprc_png), bbox_inches='tight')
        fig_auprc.savefig(str(auprc_pdf), bbox_inches='tight')



    def all_algorithm_scores_combined_pathways_boxplot(
            self, evaluation_dir, visualization_dir): 
        '''
        Make a single boxplot figure, with one boxplot per algorithm, where
        a constituent point is an algorithm's score on a particular fold.
        '''

        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        avg_prec_algorithm_map = {}
        f1_max_algorithm_map = {}
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            avg_prec_algorithm_map[name] = []
            f1_max_algorithm_map[name] = []
            auprc_algorithm_map[name] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []

                for j, fold in enumerate(test_fold):
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    avg_prec_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "average_precision.txt") 

                    f1_max_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "f1_max.txt") 

                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc.txt") 

                    avg_prec = None
                    f1_max = None
                    auprc = None

                    with avg_prec_file.open('r') as f:
                        line = next(f)
                        avg_prec = float(line.strip())
                        
                    with f1_max_file.open('r') as f:
                        line = next(f)
                        f1_max = float(line.strip())

                    with auprc_file.open('r') as f:
                        line = next(f)
                        auprc = float(line.strip())

                    name = algorithm.get_descriptive_name()

                    avg_prec_algorithm_map[name].append(avg_prec) 
                    f1_max_algorithm_map[name].append(f1_max)
                    auprc_algorithm_map[name].append(auprc)

        labels = []
        results_avg_prec = []
        results_f1_max = []
        results_auprc = []

        for alg in self.algorithms:
            name = alg.get_descriptive_name()
            labels.append(name)
            results_avg_prec.append(avg_prec_algorithm_map[name])
            results_f1_max.append(f1_max_algorithm_map[name])
            results_auprc.append(auprc_algorithm_map[name])

        fig_avg_prec, ax_avg_prec = precrec.init_precision_recall_figure()
        fig_f1_max, ax_f1_max = precrec.init_precision_recall_figure()
        fig_auprc, ax_auprc = precrec.init_precision_recall_figure()

        ax_avg_prec.set_title(
            "Average Precision by Algorithm "
            + self.interactome.name + " "
            + self.collection.name + "\n"
            + self.get_details())

        ax_avg_prec.set_xlabel("Average Precision")
        ax_avg_prec.set_ylabel("Algorithm")

        ax_f1_max.set_title(
            "Max F-Measure by Algorithm "
            + self.interactome.name + " "
            + self.collection.name + "\n"
            + self.get_details())

        ax_f1_max.set_xlabel("Max F-Measure")
        ax_f1_max.set_ylabel("Algorithm")

        ax_auprc.set_title(
            "AUPRC by Algorithm "
            + self.interactome.name + " "
            + self.collection.name + "\n"
            + self.get_details())

        ax_auprc.set_xlabel("AUPRC")
        ax_auprc.set_ylabel("Algorithm")

        avg_prec_png = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "avg_prec.png")

        avg_prec_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "avg_prec.pdf")

        f1_max_png = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "f1_max.png")

        f1_max_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "f1_max.pdf")

        auprc_png = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "auprc.png")

        auprc_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "auprc.pdf")

        avg_prec_png.parent.mkdir(parents=True, exist_ok=True)

        ax_avg_prec.boxplot(results_avg_prec, labels=labels, vert=False)
        ax_f1_max.boxplot(results_f1_max, labels=labels, vert=False)
        ax_auprc.boxplot(results_auprc, labels=labels, vert=False)

        fig_avg_prec.savefig(str(avg_prec_png), bbox_inches='tight')
        fig_avg_prec.savefig(str(avg_prec_pdf), bbox_inches='tight')

        fig_f1_max.savefig(str(f1_max_png), bbox_inches='tight')
        fig_f1_max.savefig(str(f1_max_pdf), bbox_inches='tight')

        fig_auprc.savefig(str(auprc_png), bbox_inches='tight')
        fig_auprc.savefig(str(auprc_pdf), bbox_inches='tight')


    def all_algorithm_scores_individual_pathways_boxplots( 
            self, evaluation_dir, visualization_dir): 

        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        avg_prec_algorithm_map = {}
        f1_max_algorithm_map = {}
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)
                avg_prec_algorithm_map[tup] = []
                f1_max_algorithm_map[tup] = []
                auprc_algorithm_map[tup] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []

                for j, fold in enumerate(test_fold):
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    avg_prec_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "average_precision.txt") 

                    f1_max_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "f1_max.txt") 

                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc.txt") 

                    avg_prec = None
                    f1_max = None
                    auprc = None

                    with avg_prec_file.open('r') as f:
                        line = next(f)
                        avg_prec = float(line.strip())
                        
                    with f1_max_file.open('r') as f:
                        line = next(f)
                        f1_max = float(line.strip())

                    with auprc_file.open('r') as f:
                        line = next(f)
                        auprc = float(line.strip())

                    alg_name = algorithm.get_descriptive_name()

                    tup = (alg_name, pathway.name)

                    avg_prec_algorithm_map[tup].append(avg_prec) 
                    f1_max_algorithm_map[tup].append(f1_max)
                    auprc_algorithm_map[tup].append(auprc)

        for pathway in self.collection.pathways:
            labels = []
            results_avg_prec = []
            results_f1_max = []
            results_auprc = []

            for alg in self.algorithms:
                name = alg.get_descriptive_name()
                labels.append(name)
                tup = (name, pathway.name)
                results_avg_prec.append(avg_prec_algorithm_map[tup])
                results_f1_max.append(f1_max_algorithm_map[tup])
                results_auprc.append(auprc_algorithm_map[tup])

            fig_avg_prec, ax_avg_prec = precrec.init_precision_recall_figure()
            fig_f1_max, ax_f1_max = precrec.init_precision_recall_figure()
            fig_auprc, ax_auprc = precrec.init_precision_recall_figure()

            ax_avg_prec.set_title(
                "Average Precision by Algorithm (%s)" % pathway.name
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_avg_prec.set_xlabel("Average Precision")
            ax_avg_prec.set_ylabel("Algorithm")

            ax_f1_max.set_title(
                "Max F-Measure by Algorithm (%s)" % pathway.name
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_f1_max.set_xlabel("Max F-Measure")
            ax_f1_max.set_ylabel("Algorithm")

            ax_auprc.set_title(
                "AUPRC by Algorithm (%s)" % pathway.name
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_auprc.set_xlabel("AUPRC")
            ax_auprc.set_ylabel("Algorithm")

            avg_prec_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "avg_prec.png")

            avg_prec_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "avg_prec.pdf")

            f1_max_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "f1_max.png")

            f1_max_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "f1_max.pdf")

            auprc_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "auprc.png")

            auprc_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "auprc.pdf")

            avg_prec_png.parent.mkdir(parents=True, exist_ok=True)

            ax_avg_prec.boxplot(results_avg_prec, labels=labels, vert=False)
            ax_f1_max.boxplot(results_f1_max, labels=labels, vert=False)
            ax_auprc.boxplot(results_auprc, labels=labels, vert=False)

            fig_avg_prec.savefig(str(avg_prec_png), bbox_inches='tight')
            fig_avg_prec.savefig(str(avg_prec_pdf), bbox_inches='tight')

            fig_f1_max.savefig(str(f1_max_png), bbox_inches='tight')
            fig_f1_max.savefig(str(f1_max_pdf), bbox_inches='tight')

            fig_auprc.savefig(str(auprc_png), bbox_inches='tight')
            fig_auprc.savefig(str(auprc_pdf), bbox_inches='tight')


    def individual_algorithm_scores_all_pathways_boxplots(
            self, evaluation_dir, visualization_dir): 

        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        avg_prec_algorithm_map = {}
        f1_max_algorithm_map = {}
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)
                avg_prec_algorithm_map[tup] = []
                f1_max_algorithm_map[tup] = []
                auprc_algorithm_map[tup] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []

                for j, fold in enumerate(test_fold):
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    avg_prec_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "average_precision.txt") 

                    f1_max_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "f1_max.txt") 

                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc.txt") 

                    avg_prec = None
                    f1_max = None
                    auprc = None

                    with avg_prec_file.open('r') as f:
                        line = next(f)
                        avg_prec = float(line.strip())
                        
                    with f1_max_file.open('r') as f:
                        line = next(f)
                        f1_max = float(line.strip())

                    with auprc_file.open('r') as f:
                        line = next(f)
                        auprc = float(line.strip())

                    alg_name = algorithm.get_descriptive_name()

                    tup = (alg_name, pathway.name)

                    avg_prec_algorithm_map[tup].append(avg_prec) 
                    f1_max_algorithm_map[tup].append(f1_max)
                    auprc_algorithm_map[tup].append(auprc)

        for alg in self.algorithms:
            name = alg.get_descriptive_name()
            labels = []
            results_avg_prec = []
            results_f1_max = []
            results_auprc = []

            for pathway in self.collection.pathways:
                labels.append(pathway.name)
                tup = (name, pathway.name)
                results_avg_prec.append(avg_prec_algorithm_map[tup])
                results_f1_max.append(f1_max_algorithm_map[tup])
                results_auprc.append(auprc_algorithm_map[tup])

            fig_avg_prec, ax_avg_prec = precrec.init_precision_recall_figure()
            fig_f1_max, ax_f1_max = precrec.init_precision_recall_figure()
            fig_auprc, ax_auprc = precrec.init_precision_recall_figure()

            ax_avg_prec.set_title(
                "Average Precision by Pathway"
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_avg_prec.set_xlabel("Average Precision")
            ax_avg_prec.set_ylabel("Pathway")

            ax_f1_max.set_title(
                "Max F-Measure by Pathway"
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_f1_max.set_xlabel("Max F-Measure")
            ax_f1_max.set_ylabel("Pathway")

            ax_auprc.set_title(
                "AUPRC by Pathway"
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_auprc.set_xlabel("AUPRC")
            ax_auprc.set_ylabel("Pathway")

            avg_prec_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-avg_prec.png")

            avg_prec_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-avg_prec.pdf")

            f1_max_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-f1_max.png")

            f1_max_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-f1_max.pdf")

            auprc_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-auprc.png")

            auprc_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-auprc.pdf")

            avg_prec_png.parent.mkdir(parents=True, exist_ok=True)

            ax_avg_prec.boxplot(results_avg_prec, labels=labels, vert=False)
            ax_f1_max.boxplot(results_f1_max, labels=labels, vert=False)
            ax_auprc.boxplot(results_auprc, labels=labels, vert=False)

            fig_avg_prec.savefig(str(avg_prec_png), bbox_inches='tight')
            fig_avg_prec.savefig(str(avg_prec_pdf), bbox_inches='tight')

            fig_f1_max.savefig(str(f1_max_png), bbox_inches='tight')
            fig_f1_max.savefig(str(f1_max_pdf), bbox_inches='tight')

            fig_auprc.savefig(str(auprc_png), bbox_inches='tight')
            fig_auprc.savefig(str(auprc_pdf), bbox_inches='tight')



    def all_algorithm_scores_combined_pathways_boxplot_nodes(
            self, evaluation_dir, visualization_dir): 
        '''
        Make a single boxplot figure, with one boxplot per algorithm, where
        a constituent point is an algorithm's score on a particular fold.
        '''

        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        avg_prec_algorithm_map = {}
        f1_max_algorithm_map = {}
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            avg_prec_algorithm_map[name] = []
            f1_max_algorithm_map[name] = []
            auprc_algorithm_map[name] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []

                for j, fold in enumerate(test_fold):
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    avg_prec_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "average_precision-nodes.txt") 

                    f1_max_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "f1_max-nodes.txt") 

                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc-nodes.txt") 

                    avg_prec = None
                    f1_max = None
                    auprc = None

                    with avg_prec_file.open('r') as f:
                        line = next(f)
                        avg_prec = float(line.strip())
                        
                    with f1_max_file.open('r') as f:
                        line = next(f)
                        f1_max = float(line.strip())

                    with auprc_file.open('r') as f:
                        line = next(f)
                        auprc = float(line.strip())

                    name = algorithm.get_descriptive_name()

                    avg_prec_algorithm_map[name].append(avg_prec) 
                    f1_max_algorithm_map[name].append(f1_max)
                    auprc_algorithm_map[name].append(auprc)

        labels = []
        results_avg_prec = []
        results_f1_max = []
        results_auprc = []

        for alg in self.algorithms:
            name = alg.get_descriptive_name()
            labels.append(name)
            results_avg_prec.append(avg_prec_algorithm_map[name])
            results_f1_max.append(f1_max_algorithm_map[name])
            results_auprc.append(auprc_algorithm_map[name])

        fig_avg_prec, ax_avg_prec = precrec.init_precision_recall_figure()
        fig_f1_max, ax_f1_max = precrec.init_precision_recall_figure()
        fig_auprc, ax_auprc = precrec.init_precision_recall_figure()

        ax_avg_prec.set_title(
            "Average Precision by Algorithm "
            + self.interactome.name + " "
            + self.collection.name + "\n"
            + self.get_details())

        ax_avg_prec.set_xlabel("Average Precision")
        ax_avg_prec.set_ylabel("Algorithm")

        ax_f1_max.set_title(
            "Max F-Measure by Algorithm "
            + self.interactome.name + " "
            + self.collection.name + "\n"
            + self.get_details())

        ax_f1_max.set_xlabel("Max F-Measure")
        ax_f1_max.set_ylabel("Algorithm")

        ax_auprc.set_title(
            "AUPRC by Algorithm "
            + self.interactome.name + " "
            + self.collection.name + "\n"
            + self.get_details())

        ax_auprc.set_xlabel("AUPRC")
        ax_auprc.set_ylabel("Algorithm")

        avg_prec_png = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "avg_prec-nodes.png")

        avg_prec_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "avg_prec-nodes.pdf")

        f1_max_png = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "f1_max-nodes.png")

        f1_max_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "f1_max-nodes.pdf")

        auprc_png = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "auprc-nodes.png")

        auprc_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.collection.name,
            self.get_output_prefix(),
            "auprc-nodes.pdf")

        avg_prec_png.parent.mkdir(parents=True, exist_ok=True)

        ax_avg_prec.boxplot(results_avg_prec, labels=labels, vert=False)
        ax_f1_max.boxplot(results_f1_max, labels=labels, vert=False)
        ax_auprc.boxplot(results_auprc, labels=labels, vert=False)

        fig_avg_prec.savefig(str(avg_prec_png), bbox_inches='tight')
        fig_avg_prec.savefig(str(avg_prec_pdf), bbox_inches='tight')

        fig_f1_max.savefig(str(f1_max_png), bbox_inches='tight')
        fig_f1_max.savefig(str(f1_max_pdf), bbox_inches='tight')

        fig_auprc.savefig(str(auprc_png), bbox_inches='tight')
        fig_auprc.savefig(str(auprc_pdf), bbox_inches='tight')


    def all_algorithm_scores_individual_pathways_boxplots_nodes(
            self, evaluation_dir, visualization_dir): 

        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        avg_prec_algorithm_map = {}
        f1_max_algorithm_map = {}
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)
                avg_prec_algorithm_map[tup] = []
                f1_max_algorithm_map[tup] = []
                auprc_algorithm_map[tup] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []

                for j, fold in enumerate(test_fold):
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    avg_prec_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "average_precision-nodes.txt") 

                    f1_max_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "f1_max-nodes.txt") 

                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc-nodes.txt") 

                    avg_prec = None
                    f1_max = None
                    auprc = None

                    with avg_prec_file.open('r') as f:
                        line = next(f)
                        avg_prec = float(line.strip())
                        
                    with f1_max_file.open('r') as f:
                        line = next(f)
                        f1_max = float(line.strip())

                    with auprc_file.open('r') as f:
                        line = next(f)
                        auprc = float(line.strip())

                    alg_name = algorithm.get_descriptive_name()

                    tup = (alg_name, pathway.name)

                    avg_prec_algorithm_map[tup].append(avg_prec) 
                    f1_max_algorithm_map[tup].append(f1_max)
                    auprc_algorithm_map[tup].append(auprc)

        for pathway in self.collection.pathways:
            labels = []
            results_avg_prec = []
            results_f1_max = []
            results_auprc = []

            for alg in self.algorithms:
                name = alg.get_descriptive_name()
                labels.append(name)
                tup = (name, pathway.name)
                results_avg_prec.append(avg_prec_algorithm_map[tup])
                results_f1_max.append(f1_max_algorithm_map[tup])
                results_auprc.append(auprc_algorithm_map[tup])

            fig_avg_prec, ax_avg_prec = precrec.init_precision_recall_figure()
            fig_f1_max, ax_f1_max = precrec.init_precision_recall_figure()
            fig_auprc, ax_auprc = precrec.init_precision_recall_figure()

            ax_avg_prec.set_title(
                "Average Precision by Algorithm (%s)" % pathway.name
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_avg_prec.set_xlabel("Average Precision")
            ax_avg_prec.set_ylabel("Algorithm")

            ax_f1_max.set_title(
                "Max F-Measure by Algorithm (%s)" % pathway.name
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_f1_max.set_xlabel("Max F-Measure")
            ax_f1_max.set_ylabel("Algorithm")

            ax_auprc.set_title(
                "AUPRC by Algorithm (%s)" % pathway.name
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_auprc.set_xlabel("AUPRC")
            ax_auprc.set_ylabel("Algorithm")

            avg_prec_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "avg_prec-nodes.png")

            avg_prec_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "avg_prec-nodes.pdf")

            f1_max_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "f1_max-nodes.png")

            f1_max_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "f1_max-nodes.pdf")

            auprc_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "auprc-nodes.png")

            auprc_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "auprc-nodes.pdf")

            avg_prec_png.parent.mkdir(parents=True, exist_ok=True)

            ax_avg_prec.boxplot(results_avg_prec, labels=labels, vert=False)
            ax_f1_max.boxplot(results_f1_max, labels=labels, vert=False)
            ax_auprc.boxplot(results_auprc, labels=labels, vert=False)

            fig_avg_prec.savefig(str(avg_prec_png), bbox_inches='tight')
            fig_avg_prec.savefig(str(avg_prec_pdf), bbox_inches='tight')

            fig_f1_max.savefig(str(f1_max_png), bbox_inches='tight')
            fig_f1_max.savefig(str(f1_max_pdf), bbox_inches='tight')

            fig_auprc.savefig(str(auprc_png), bbox_inches='tight')
            fig_auprc.savefig(str(auprc_pdf), bbox_inches='tight')


    def individual_algorithm_scores_all_pathways_boxplots_nodes(
            self, evaluation_dir, visualization_dir): 

        pairs = zip(self.collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        avg_prec_algorithm_map = {}
        f1_max_algorithm_map = {}
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)
                avg_prec_algorithm_map[tup] = []
                f1_max_algorithm_map[tup] = []
                auprc_algorithm_map[tup] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []

                for j, fold in enumerate(test_fold):
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))

                    avg_prec_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "average_precision-nodes.txt") 

                    f1_max_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "f1_max-nodes.txt") 

                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc-nodes.txt") 

                    avg_prec = None
                    f1_max = None
                    auprc = None

                    with avg_prec_file.open('r') as f:
                        line = next(f)
                        avg_prec = float(line.strip())
                        
                    with f1_max_file.open('r') as f:
                        line = next(f)
                        f1_max = float(line.strip())

                    with auprc_file.open('r') as f:
                        line = next(f)
                        auprc = float(line.strip())

                    alg_name = algorithm.get_descriptive_name()

                    tup = (alg_name, pathway.name)

                    avg_prec_algorithm_map[tup].append(avg_prec) 
                    f1_max_algorithm_map[tup].append(f1_max)
                    auprc_algorithm_map[tup].append(auprc)

        for alg in self.algorithms:
            name = alg.get_descriptive_name()
            labels = []
            results_avg_prec = []
            results_f1_max = []
            results_auprc = []

            for pathway in self.collection.pathways:
                labels.append(pathway.name)
                tup = (name, pathway.name)
                results_avg_prec.append(avg_prec_algorithm_map[tup])
                results_f1_max.append(f1_max_algorithm_map[tup])
                results_auprc.append(auprc_algorithm_map[tup])

            fig_avg_prec, ax_avg_prec = precrec.init_precision_recall_figure()
            fig_f1_max, ax_f1_max = precrec.init_precision_recall_figure()
            fig_auprc, ax_auprc = precrec.init_precision_recall_figure()

            ax_avg_prec.set_title(
                "Average Precision by Pathway"
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_avg_prec.set_xlabel("Average Precision")
            ax_avg_prec.set_ylabel("Pathway")

            ax_f1_max.set_title(
                "Max F-Measure by Pathway"
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_f1_max.set_xlabel("Max F-Measure")
            ax_f1_max.set_ylabel("Pathway")

            ax_auprc.set_title(
                "AUPRC by Pathway"
                + self.interactome.name + " "
                + self.collection.name + "\n"
                + self.get_details())

            ax_auprc.set_xlabel("AUPRC")
            ax_auprc.set_ylabel("Pathway")

            avg_prec_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-avg_prec-nodes.png")

            avg_prec_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-avg_prec-nodes.pdf")

            f1_max_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-f1_max-nodes.png")

            f1_max_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-f1_max-nodes.pdf")

            auprc_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-auprc-nodes.png")

            auprc_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                name + "-auprc-nodes.pdf")

            avg_prec_png.parent.mkdir(parents=True, exist_ok=True)

            ax_avg_prec.boxplot(results_avg_prec, labels=labels, vert=False)
            ax_f1_max.boxplot(results_f1_max, labels=labels, vert=False)
            ax_auprc.boxplot(results_auprc, labels=labels, vert=False)

            fig_avg_prec.savefig(str(avg_prec_png), bbox_inches='tight')
            fig_avg_prec.savefig(str(avg_prec_pdf), bbox_inches='tight')

            fig_f1_max.savefig(str(f1_max_png), bbox_inches='tight')
            fig_f1_max.savefig(str(f1_max_pdf), bbox_inches='tight')

            fig_auprc.savefig(str(auprc_png), bbox_inches='tight')
            fig_auprc.savefig(str(auprc_pdf), bbox_inches='tight')



    def plot_pr_individual_pathways(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        for pathway in self.collection.pathways:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " +
                self.collection.name + " " +
                pathway.name)

            # Where we wrote precision/recall results
            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "aggregate")

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "precision-recall.pdf")
    
            # PNG file we will write 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "precision-recall.png")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                
                pr_file = Path(
                    pr_output_dir,
                    algorithm.get_output_directory(),
                    "precision-recall.txt")

                points = []

                with pr_file.open('r') as f:
                    points = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    points, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

            fig.savefig(str(vis_file_png), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')


    def plot_pr_all_pathways(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        for pathway in self.collection.pathways:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " +
                self.collection.name + " " +
                self.get_details())

            # Where we wrote precision/recall results
            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.collection.name,
                "aggregate",
                self.get_output_prefix())

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                "precision-recall.pdf")
    
            # PNG file we will write 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                "precision-recall.png")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                
                pr_file = Path(
                    pr_output_dir,
                    algorithm.get_output_directory(),
                    "precision-recall.txt")

                points = []

                with pr_file.open('r') as f:
                    points = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    points, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

            fig.savefig(str(vis_file_png), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')


    def plot_pr_individual_pathways_nodes(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        for pathway in self.collection.pathways:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " +
                self.collection.name + " " +
                pathway.name)

            # Where we wrote precision/recall results
            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "aggregate")

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "precision-recall-nodes.pdf")
    
            # PNG file we will write 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                pathway.name,
                self.get_output_prefix(),
                "precision-recall-nodes.png")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                
                pr_file = Path(
                    pr_output_dir,
                    algorithm.get_output_directory(),
                    "precision-recall-nodes.txt")

                points = []

                with pr_file.open('r') as f:
                    points = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    points, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

            fig.savefig(str(vis_file_png), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')


    def plot_pr_all_pathways_nodes(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        for pathway in self.collection.pathways:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " +
                self.collection.name + " " +
                self.get_details())

            # Where we wrote precision/recall results
            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.collection.name,
                "aggregate",
                self.get_output_prefix())

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                "precision-recall-nodes.pdf")
    
            # PNG file we will write 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                "precision-recall-nodes.png")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                
                pr_file = Path(
                    pr_output_dir,
                    algorithm.get_output_directory(),
                    "precision-recall-nodes.txt")

                points = []

                with pr_file.open('r') as f:
                    points = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    points, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

            fig.savefig(str(vis_file_png), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')


    def successively_aggregate_pr_by_pathway_auprc(self, evaluation_dir=Path(),
            visualization_dir="outputs/pathway-reconstruction/visualization"):

        '''
        For each algorithm:
            - Order pathways by median AUPRC
            - Plot aggregated precision recall for 1st, 2nd pathway
            - Then plot for 1st + 2nd + 3rd
            - Etc.

        Note: this method will have to be called subsequent to the calculation
        of AUPRC per pathway.
        '''

        print("Aggregating successively by AUPRC per pathway")

        # First, read in the AUPRC scores per pathway

        # Initialize map of algorithms to lists of points
        auprc_algorithm_map = {}

        median_map = {}
        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            median_map[name] = []

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)
                auprc_algorithm_map[tup] = []

        # Now read in the fold AUPRC per pathway

        pairs = zip(self.collection.pathways, self.test_folds)

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                
                eval_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pr_file = Path(
                    eval_dir,
                    algorithm.get_output_directory(),
                    "precision-recall.txt")

                with pr_file.open('r') as f:
                    pr_points = precrec.read_precision_recall_fractions(f)

                
                # AGHGHHGH I have to convert these to DECIMALS TO USE MY
                # AUPRC FUNCTION

                recall = [float(p[0]) / float(p[2]) for p in pr_points]
                precision = [float(p[0]) / float(p[1]) if p[1] != 0 else 1 for p in pr_points]

                decimals = list(zip(precision, recall))

                # Then get the AUPRC
                auprc = precrec.compute_AUPRC(decimals)

                # Then append it as appropriate
                
                alg_name = algorithm.get_descriptive_name()
                tup = (alg_name, pathway.name)
                auprc_algorithm_map[tup].append(auprc)


        # Now make a sorted list of pathways by median AUPRC per algorithm

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)

                #median = np.median(auprc_algorithm_map[tup])
                val = auprc_algorithm_map[tup]
                median_map[name].append((pathway.name, val))

        
        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            median_map[name] = sorted(
                median_map[name], 
                key=lambda x: x[1],
                reverse=True)

            print(median_map[name])


        print("-----------------------------------------------------------")
        print("Writing algorithm plots to show performance on each pathway")
        print("-----------------------------------------------------------")

        for algorithm in self.algorithms:

            alg_name = algorithm.get_descriptive_name()

            # Create a plot here
            pathways_fig, pathways_ax = precrec.init_precision_recall_figure()

            pathways_ax.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" + alg_name
                + self.get_details())

            # PDF file we will write
            overall_vis = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name,
                "per-pathway.pdf")

            for pathway in self.collection.pathways:
                
                curve = []

                # Where we wrote precision/recall, aggregated over
                # all folds per pathway
                pathway_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pathway_pr_outfile = Path(
                    pathway_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                with pathway_pr_outfile.open('r') as f:
                    curve = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    curve, label=pathway.name, ax=pathways_ax)

            pathways_ax.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            overall_vis.parent.mkdir(parents=True, exist_ok=True)

            handles, labels = pathways_ax.get_legend_handles_labels()

            lgd = pathways_ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            pathways_fig.savefig(str(overall_vis), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')
                




        # Now I have exactly what I need to determine which pathways to do
        # things with

        alg_to_groupings = {}


        for algorithm in self.algorithms:
            alg_name = algorithm.get_descriptive_name()
            
            ordered_pathways = median_map[alg_name]
            
            # Now I can get a list of lists of (name, median) tuples...

            # Given: [1 2 3 4]
            # I want: [[1], [1, 2], [1, 2, 3], [1, 2, 3, 4]]
            
            # Uncomment this line to order by best pathway per algorithm
            successive_groupings = [ordered_pathways[:i] for i in range(1, len(ordered_pathways) + 1)]

            # Uncomment these lines to order by pathway edge count
            #pathways_by_edge_count = [
            #    ("EGFR1", "placeholder"),
            #    ("TNFalpha", "placeholder"),
            #    ("TGF_beta_Receptor", "placeholder"),
            #    ("TCR", "placeholder"),
            #    ("Wnt", "placeholder"),
            #    ("IL2", "placeholder"),
            #    ("Prolactin", "placeholder"),
            #    ("KitReceptor", "placeholder"),
            #    ("IL1", "placeholder"),
            #    ("IL6", "placeholder"),
            #    ("IL3", "placeholder"),
            #    ("Leptin", "placeholder"),
            #    ("RANKL", "placeholder"),
            #    ("BDNF", "placeholder"),
            #    ("IL-7", "placeholder"),
            #    ]

            #successive_groupings = [pathways_by_edge_count[:i] for i in range(1, len(pathways_by_edge_count) + 1)]

            # Uncomment these lines to order by pathway node count
            #

            alg_to_groupings[alg_name] = successive_groupings

            name_file = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name,
                "pathway-order.txt")

            name_file.parent.mkdir(exist_ok=True, parents=True)

            with name_file.open('w') as f:
                last_grouping = successive_groupings[-1]
                names = [tup[0] for tup in last_grouping]
                for name in names:
                    f.write(name)
                    f.write("\n")

            pairs = list(zip(
                self.collection.pathways, self.test_folds))

            fig2, ax2 = precrec.init_precision_recall_figure()

            ax2.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            # PDF file we will write
            overall_vis = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name,
                "overall.pdf")

            for i, grouping in enumerate(successive_groupings):
                names = [tup[0] for tup in grouping]
                print(names)
    
                pair_subset = []

                for name in names:
                    # This list comprehension should only return a list with
                    # a single item. Poor man's filter to preserve order of
                    # items in the names list
                    pair_subset = pair_subset + \
                        [tup for tup in pairs if tup[0].name == name]


                # Where we will write precision/recall, aggregated over
                # pathways
                pathway_collection_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    "aggregate",
                    self.get_output_prefix(),
                    alg_name)

                edge_curves = []

                # I probably don't actually need the test_fold part of this
                # tuple
                for pathway, test_fold in pair_subset:

                    # Where we wrote precision/recall, aggregated over
                    # all folds per pathway
                    pathway_pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "aggregate")

                    pathway_pr_outfile = Path(
                        pathway_pr_output_dir, 
                        algorithm.get_output_directory(),
                        "precision-recall.txt") 

                    with pathway_pr_outfile.open('r') as f:
                        curve = precrec.read_precision_recall_fractions(f)
                        edge_curves.append(curve)

                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions(
                    edge_curves)
                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions_stretch(
                    edge_curves)

                # Write aggregated curve back out
                pathway_collection_pr_outfile = Path(
                    pathway_collection_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                pathway_collection_pr_outfile.parent.mkdir(
                    parents=True, exist_ok=True)

                with pathway_collection_pr_outfile.open("w") as f: 
                    precrec.write_precision_recall_fractions(f, aggregated_edges)

                ################### FIRST plot the edges
                fig, ax = precrec.init_precision_recall_figure()

                ax.set_title(
                    self.interactome.name + " " + self.collection.name 
                    + " "
                    + "\n" 
                    + alg_name
                    + self.get_details())

                # PDF file we will write
                vis_file_pdf = Path(
                    visualization_dir,
                    self.interactome.name,
                    self.collection.name,
                    self.get_output_prefix(),
                    alg_name,
                    "precision-recall-top-" + str(i) + ".pdf")
                
                vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

                precrec.plot_precision_recall_curve_fractions(
                    aggregated_edges, label=algorithm.get_descriptive_name(), ax=ax)

                handles, labels = ax.get_legend_handles_labels()

                lgd = ax.legend(handles, labels, loc='upper center', 
                    bbox_to_anchor=(0.5,-0.1))

                fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')
                
                precrec.plot_precision_recall_curve_fractions(
                    aggregated_edges, label="+" + pathway.name, ax=ax2)
            
            overall_vis.parent.mkdir(parents=True, exist_ok=True)

            handles, labels = ax2.get_legend_handles_labels()

            lgd = ax2.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig2.savefig(str(overall_vis), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

        

        # THIS PLOTS ALL ALGORITHMS ON THEIR BEST, THEN BEST + SECOND BEST,
        # ETC.
            
        # There are 15 pathways, so 15 groupings
        for grouping in range(len(self.collection.pathways)):

            # First, initialize a figure for the grouping. We want 15 figures
            # when we are done
        
            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                "aggregated",
                "precision-recall-top-" + str(grouping) + ".pdf")
            
            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                alg_name = algorithm.get_descriptive_name()

                pairs = list(zip(
                    self.collection.pathways, self.test_folds))

                group = alg_to_groupings[alg_name][grouping] 
                names = [tup[0] for tup in group]
                pair_subset = []

                for name in names:
                    # This list comprehension should only return a list with
                    # a single item. Poor man's filter to preserve order of
                    # items in the names list
                    pair_subset = pair_subset + \
                        [tup for tup in pairs if tup[0].name == name]

                edge_curves = []

                for pathway, _ in pair_subset:

                    # Where we wrote precision/recall, aggregated over
                    # all folds per pathway
                    pathway_pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "aggregate")

                    pathway_pr_outfile = Path(
                        pathway_pr_output_dir, 
                        algorithm.get_output_directory(),
                        "precision-recall.txt") 

                    with pathway_pr_outfile.open('r') as f:
                        curve = precrec.read_precision_recall_fractions(f)
                        edge_curves.append(curve)

                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions(
                    edge_curves)
                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions_stretch(
                    edge_curves)

                precrec.plot_precision_recall_curve_fractions(
                    aggregated_edges, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')








    def successively_aggregate_pr_by_pathway_auprc_nodes(self, evaluation_dir=Path(),
            visualization_dir="outputs/pathway-reconstruction/visualization"):

        '''
        For each algorithm:
            - Order pathways by median AUPRC
            - Plot aggregated precision recall for 1st, 2nd pathway
            - Then plot for 1st + 2nd + 3rd
            - Etc.

        Note: this method will have to be called subsequent to the calculation
        of AUPRC per pathway.
        '''

        print("Aggregating successively by AUPRC per pathway")

        # First, read in the AUPRC scores per pathway

        # Initialize map of algorithms to lists of points
        auprc_algorithm_map = {}

        median_map = {}
        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            median_map[name] = []

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)
                auprc_algorithm_map[tup] = []

        # Now read in the fold AUPRC per pathway

        pairs = zip(self.collection.pathways, self.test_folds)

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                
                eval_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pr_file = Path(
                    eval_dir,
                    algorithm.get_output_directory(),
                    "precision-recall-nodes.txt")

                with pr_file.open('r') as f:
                    pr_points = precrec.read_precision_recall_fractions(f)

                
                # AGHGHHGH I have to convert these to DECIMALS TO USE MY
                # AUPRC FUNCTION

                recall = [float(p[0]) / float(p[2]) for p in pr_points]
                precision = [float(p[0]) / float(p[1]) if p[1] != 0 else 1 for p in pr_points]

                decimals = list(zip(precision, recall))

                # Then get the AUPRC
                auprc = precrec.compute_AUPRC(decimals)

                # Then append it as appropriate
                
                alg_name = algorithm.get_descriptive_name()
                tup = (alg_name, pathway.name)
                auprc_algorithm_map[tup].append(auprc)



        # Now make a sorted list of pathways by median AUPRC per algorithm

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)

                #median = np.median(auprc_algorithm_map[tup])
                val = auprc_algorithm_map[tup]
                median_map[name].append((pathway.name, val))

        
        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            median_map[name] = sorted(
                median_map[name], 
                key=lambda x: x[1],
                reverse=True)

            print(median_map[name])


        print("-----------------------------------------------------------")
        print("Writing algorithm plots to show performance on each pathway")
        print("-----------------------------------------------------------")

        for algorithm in self.algorithms:

            alg_name = algorithm.get_descriptive_name()

            # Create a plot here
            pathways_fig, pathways_ax = precrec.init_precision_recall_figure()

            pathways_ax.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" + alg_name
                + self.get_details())

            # PDF file we will write
            overall_vis = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name + "-nodes",
                "per-pathway.pdf")

            for pathway in self.collection.pathways:
                
                curve = []

                # Where we wrote precision/recall, aggregated over
                # all folds per pathway
                pathway_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pathway_pr_outfile = Path(
                    pathway_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall-nodes.txt") 

                with pathway_pr_outfile.open('r') as f:
                    curve = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    curve, label=pathway.name, ax=pathways_ax)

            pathways_ax.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            overall_vis.parent.mkdir(parents=True, exist_ok=True)

            handles, labels = pathways_ax.get_legend_handles_labels()

            lgd = pathways_ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            pathways_fig.savefig(str(overall_vis), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')
                




        # Now I have exactly what I need to determine which pathways to do
        # things with

        alg_to_groupings = {}


        for algorithm in self.algorithms:
            alg_name = algorithm.get_descriptive_name()
            
            ordered_pathways = median_map[alg_name]
            
            # Now I can get a list of lists of (name, median) tuples...

            # Given: [1 2 3 4]
            # I want: [[1], [1, 2], [1, 2, 3], [1, 2, 3, 4]]
            
            # Uncomment this line to order by best pathway per algorithm
            successive_groupings = [ordered_pathways[:i] for i in range(1, len(ordered_pathways) + 1)]

            # Uncomment these lines to order by pathway edge count
            #pathways_by_edge_count = [
            #    ("EGFR1", "placeholder"),
            #    ("TNFalpha", "placeholder"),
            #    ("TGF_beta_Receptor", "placeholder"),
            #    ("TCR", "placeholder"),
            #    ("Wnt", "placeholder"),
            #    ("IL2", "placeholder"),
            #    ("Prolactin", "placeholder"),
            #    ("KitReceptor", "placeholder"),
            #    ("IL1", "placeholder"),
            #    ("IL6", "placeholder"),
            #    ("IL3", "placeholder"),
            #    ("Leptin", "placeholder"),
            #    ("RANKL", "placeholder"),
            #    ("BDNF", "placeholder"),
            #    ("IL-7", "placeholder"),
            #    ]

            #successive_groupings = [pathways_by_edge_count[:i] for i in range(1, len(pathways_by_edge_count) + 1)]

            # Uncomment these lines to order by pathway node count
            #

            alg_to_groupings[alg_name] = successive_groupings

            name_file = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name + "-nodes",
                "pathway-order.txt")

            name_file.parent.mkdir(exist_ok=True, parents=True)

            with name_file.open('w') as f:
                last_grouping = successive_groupings[-1]
                names = [tup[0] for tup in last_grouping]
                for name in names:
                    f.write(name)
                    f.write("\n")

            pairs = list(zip(
                self.collection.pathways, self.test_folds))

            fig2, ax2 = precrec.init_precision_recall_figure()

            ax2.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            # PDF file we will write
            overall_vis = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name + "-nodes",
                "overall.pdf")

            for i, grouping in enumerate(successive_groupings):
                names = [tup[0] for tup in grouping]
                print(names)
    

                # pair_subset = [tup for tup in pairs if tup[0].name in names]

                # No... this will use the order of pathways in pairs

                # Instead, I need to somehow order the things by their position
                # in the names array. There are lots of easy ways to do this!

                # Heck, let's do it the dumb way

                pair_subset = []

                for name in names:
                    # This list comprehension should only return a list with
                    # a single item. Poor man's filter to preserve order of
                    # items in the names list
                    pair_subset = pair_subset + \
                        [tup for tup in pairs if tup[0].name == name]


                # Where we will write precision/recall, aggregated over
                # pathways
                pathway_collection_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    "aggregate",
                    self.get_output_prefix(),
                    alg_name)

                edge_curves = []

                # I probably don't actually need the test_fold part of this
                # tuple
                for pathway, test_fold in pair_subset:

                    # Where we wrote precision/recall, aggregated over
                    # all folds per pathway
                    pathway_pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "aggregate")

                    pathway_pr_outfile = Path(
                        pathway_pr_output_dir, 
                        algorithm.get_output_directory(),
                        "precision-recall-nodes.txt") 

                    with pathway_pr_outfile.open('r') as f:
                        curve = precrec.read_precision_recall_fractions(f)
                        edge_curves.append(curve)

                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions(
                    edge_curves)
                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions_stretch(
                    edge_curves)

                # Write aggregated curve back out
                pathway_collection_pr_outfile = Path(
                    pathway_collection_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall-nodes.txt") 

                pathway_collection_pr_outfile.parent.mkdir(
                    parents=True, exist_ok=True)

                with pathway_collection_pr_outfile.open("w") as f: 
                    precrec.write_precision_recall_fractions(f, aggregated_edges)

                ################### FIRST plot the edges
                fig, ax = precrec.init_precision_recall_figure()

                ax.set_title(
                    self.interactome.name + " " + self.collection.name 
                    + " "
                    + "\n" 
                    + alg_name
                    + self.get_details())

                # PDF file we will write
                vis_file_pdf = Path(
                    visualization_dir,
                    self.interactome.name,
                    self.collection.name,
                    self.get_output_prefix(),
                    alg_name + "-nodes",
                    "precision-recall-top-" + str(i) + ".pdf")
                
                vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

                precrec.plot_precision_recall_curve_fractions(
                    aggregated_edges, label=algorithm.get_descriptive_name(), ax=ax)

                handles, labels = ax.get_legend_handles_labels()

                lgd = ax.legend(handles, labels, loc='upper center', 
                    bbox_to_anchor=(0.5,-0.1))

                fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')
                
                precrec.plot_precision_recall_curve_fractions(
                    aggregated_edges, label="+" + pathway.name, ax=ax2)
            
            overall_vis.parent.mkdir(parents=True, exist_ok=True)

            handles, labels = ax2.get_legend_handles_labels()

            lgd = ax2.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig2.savefig(str(overall_vis), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

        

        # THIS PLOTS ALL ALGORITHMS ON THEIR BEST, THEN BEST + SECOND BEST,
        # ETC.
            
        # There are 15 pathways, so 15 groupings
        for grouping in range(len(self.collection.pathways)):

            # First, initialize a figure for the grouping. We want 15 figures
            # when we are done
        
            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                "aggregated-nodes",
                "precision-recall-top-" + str(grouping) + ".pdf")
            
            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                alg_name = algorithm.get_descriptive_name()

                pairs = list(zip(
                    self.collection.pathways, self.test_folds))

                group = alg_to_groupings[alg_name][grouping] 
                names = [tup[0] for tup in group]
                pair_subset = []

                for name in names:
                    # This list comprehension should only return a list with
                    # a single item. Poor man's filter to preserve order of
                    # items in the names list
                    pair_subset = pair_subset + \
                        [tup for tup in pairs if tup[0].name == name]

                edge_curves = []

                for pathway, _ in pair_subset:

                    # Where we wrote precision/recall, aggregated over
                    # all folds per pathway
                    pathway_pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "aggregate")

                    pathway_pr_outfile = Path(
                        pathway_pr_output_dir, 
                        algorithm.get_output_directory(),
                        "precision-recall-nodes.txt") 

                    with pathway_pr_outfile.open('r') as f:
                        curve = precrec.read_precision_recall_fractions(f)
                        edge_curves.append(curve)

                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions(
                    edge_curves)
                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions_stretch(
                    edge_curves)

                precrec.plot_precision_recall_curve_fractions(
                    aggregated_edges, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')



    def successively_aggregate_pr_by_auprc(self, evaluation_dir=Path(),
            visualization_dir="outputs/pathway-reconstruction/visualization"):

        '''
        For each algorithm:
            - Order pathways by median AUPRC
            - Plot aggregated precision recall for 1st, 2nd pathway
            - Then plot for 1st + 2nd + 3rd
            - Etc.

        Note: this method will have to be called subsequent to the calculation
        of AUPRC per pathway.
        '''

        print("Aggregating successively by AUPRC per pathway")

        # First, read in the AUPRC scores per pathway

        # Initialize map of algorithms to lists of points
        auprc_algorithm_map = {}

        median_map = {}
        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            median_map[name] = []

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)
                auprc_algorithm_map[tup] = []

        # Now read in the fold AUPRC per pathway

        pairs = zip(self.collection.pathways, self.test_folds)

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:

                for j, fold in enumerate(test_fold):
                    
                    # Where we wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))


                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc.txt") 

                    auprc = None

                    with auprc_file.open('r') as f:
                        line = next(f)
                        auprc = float(line.strip())

                    alg_name = algorithm.get_descriptive_name()

                    tup = (alg_name, pathway.name)

                    auprc_algorithm_map[tup].append(auprc)

        # Now make a sorted list of pathways by median AUPRC per algorithm

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)

                median = np.median(auprc_algorithm_map[tup])
                median_map[name].append((pathway.name, median))

        
        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            median_map[name] = sorted(
                median_map[name], 
                key=lambda x: x[1],
                reverse=True)

            print(median_map[name])


        print("-----------------------------------------------------------")
        print("Writing algorithm plots to show performance on each pathway")
        print("-----------------------------------------------------------")

        for algorithm in self.algorithms:

            alg_name = algorithm.get_descriptive_name()

            # Create a plot here
            pathways_fig, pathways_ax = precrec.init_precision_recall_figure()

            pathways_ax.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" + alg_name
                + self.get_details())

            # PDF file we will write
            overall_vis = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name,
                "per-pathway.pdf")

            for pathway in self.collection.pathways:
                
                curve = []

                # Where we wrote precision/recall, aggregated over
                # all folds per pathway
                pathway_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pathway_pr_outfile = Path(
                    pathway_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                with pathway_pr_outfile.open('r') as f:
                    curve = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    curve, label=pathway.name, ax=pathways_ax)

            pathways_ax.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            overall_vis.parent.mkdir(parents=True, exist_ok=True)

            handles, labels = pathways_ax.get_legend_handles_labels()

            lgd = pathways_ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            pathways_fig.savefig(str(overall_vis), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')
                




        # Now I have exactly what I need to determine which pathways to do
        # things with

        alg_to_groupings = {}


        for algorithm in self.algorithms:
            alg_name = algorithm.get_descriptive_name()
            
            ordered_pathways = median_map[alg_name]
            
            # Now I can get a list of lists of (name, median) tuples...

            # Given: [1 2 3 4]
            # I want: [[1], [1, 2], [1, 2, 3], [1, 2, 3, 4]]
            
            # Uncomment this line to order by best pathway per algorithm
            successive_groupings = [ordered_pathways[:i] for i in range(1, len(ordered_pathways) + 1)]

            # Uncomment these lines to order by pathway edge count
            #pathways_by_edge_count = [
            #    ("EGFR1", "placeholder"),
            #    ("TNFalpha", "placeholder"),
            #    ("TGF_beta_Receptor", "placeholder"),
            #    ("TCR", "placeholder"),
            #    ("Wnt", "placeholder"),
            #    ("IL2", "placeholder"),
            #    ("Prolactin", "placeholder"),
            #    ("KitReceptor", "placeholder"),
            #    ("IL1", "placeholder"),
            #    ("IL6", "placeholder"),
            #    ("IL3", "placeholder"),
            #    ("Leptin", "placeholder"),
            #    ("RANKL", "placeholder"),
            #    ("BDNF", "placeholder"),
            #    ("IL-7", "placeholder"),
            #    ]

            #successive_groupings = [pathways_by_edge_count[:i] for i in range(1, len(pathways_by_edge_count) + 1)]

            # Uncomment these lines to order by pathway node count
            #

            alg_to_groupings[alg_name] = successive_groupings

            name_file = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name,
                "pathway-order.txt")

            name_file.parent.mkdir(exist_ok=True, parents=True)

            with name_file.open('w') as f:
                last_grouping = successive_groupings[-1]
                names = [tup[0] for tup in last_grouping]
                for name in names:
                    f.write(name)
                    f.write("\n")

            pairs = list(zip(
                self.collection.pathways, self.test_folds))

            fig2, ax2 = precrec.init_precision_recall_figure()

            ax2.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            # PDF file we will write
            overall_vis = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name,
                "overall.pdf")

            for i, grouping in enumerate(successive_groupings):
                names = [tup[0] for tup in grouping]
                print(names)
    

                # pair_subset = [tup for tup in pairs if tup[0].name in names]

                # No... this will use the order of pathways in pairs

                # Instead, I need to somehow order the things by their position
                # in the names array. There are lots of easy ways to do this!

                # Heck, let's do it the dumb way

                pair_subset = []

                for name in names:
                    # This list comprehension should only return a list with
                    # a single item. Poor man's filter to preserve order of
                    # items in the names list
                    pair_subset = pair_subset + \
                        [tup for tup in pairs if tup[0].name == name]


                # Where we will write precision/recall, aggregated over
                # pathways
                pathway_collection_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    "aggregate",
                    self.get_output_prefix(),
                    alg_name)

                edge_curves = []

                # I probably don't actually need the test_fold part of this
                # tuple
                for pathway, test_fold in pair_subset:

                    # Where we wrote precision/recall, aggregated over
                    # all folds per pathway
                    pathway_pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "aggregate")

                    pathway_pr_outfile = Path(
                        pathway_pr_output_dir, 
                        algorithm.get_output_directory(),
                        "precision-recall.txt") 

                    with pathway_pr_outfile.open('r') as f:
                        curve = precrec.read_precision_recall_fractions(f)
                        edge_curves.append(curve)

                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions(
                    edge_curves)
                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions_stretch(
                    edge_curves)

                # Write aggregated curve back out
                pathway_collection_pr_outfile = Path(
                    pathway_collection_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                pathway_collection_pr_outfile.parent.mkdir(
                    parents=True, exist_ok=True)

                with pathway_collection_pr_outfile.open("w") as f: 
                    precrec.write_precision_recall_fractions(f, aggregated_edges)

                ################### FIRST plot the edges
                fig, ax = precrec.init_precision_recall_figure()

                ax.set_title(
                    self.interactome.name + " " + self.collection.name 
                    + " "
                    + "\n" 
                    + alg_name
                    + self.get_details())

                # PDF file we will write
                vis_file_pdf = Path(
                    visualization_dir,
                    self.interactome.name,
                    self.collection.name,
                    self.get_output_prefix(),
                    alg_name,
                    "precision-recall-top-" + str(i) + ".pdf")
                
                vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

                precrec.plot_precision_recall_curve_fractions(
                    aggregated_edges, label=algorithm.get_descriptive_name(), ax=ax)

                handles, labels = ax.get_legend_handles_labels()

                lgd = ax.legend(handles, labels, loc='upper center', 
                    bbox_to_anchor=(0.5,-0.1))

                fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')
                
                precrec.plot_precision_recall_curve_fractions(
                    aggregated_edges, label="+" + pathway.name, ax=ax2)
            
            overall_vis.parent.mkdir(parents=True, exist_ok=True)

            handles, labels = ax2.get_legend_handles_labels()

            lgd = ax2.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig2.savefig(str(overall_vis), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

        

        # THIS PLOTS ALL ALGORITHMS ON THEIR BEST, THEN BEST + SECOND BEST,
        # ETC.
            
        # There are 15 pathways, so 15 groupings
        for grouping in range(len(self.collection.pathways)):

            # First, initialize a figure for the grouping. We want 15 figures
            # when we are done
        
            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                "aggregated",
                "precision-recall-top-" + str(grouping) + ".pdf")
            
            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                alg_name = algorithm.get_descriptive_name()

                pairs = list(zip(
                    self.collection.pathways, self.test_folds))

                group = alg_to_groupings[alg_name][grouping] 
                names = [tup[0] for tup in group]
                pair_subset = []

                for name in names:
                    # This list comprehension should only return a list with
                    # a single item. Poor man's filter to preserve order of
                    # items in the names list
                    pair_subset = pair_subset + \
                        [tup for tup in pairs if tup[0].name == name]

                edge_curves = []

                for pathway, _ in pair_subset:

                    # Where we wrote precision/recall, aggregated over
                    # all folds per pathway
                    pathway_pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "aggregate")

                    pathway_pr_outfile = Path(
                        pathway_pr_output_dir, 
                        algorithm.get_output_directory(),
                        "precision-recall.txt") 

                    with pathway_pr_outfile.open('r') as f:
                        curve = precrec.read_precision_recall_fractions(f)
                        edge_curves.append(curve)

                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions(
                    edge_curves)
                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions_stretch(
                    edge_curves)

                precrec.plot_precision_recall_curve_fractions(
                    aggregated_edges, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')








    def successively_aggregate_pr_by_auprc_nodes(self, evaluation_dir=Path(),
            visualization_dir="outputs/pathway-reconstruction/visualization"):

        '''
        For each algorithm:
            - Order pathways by median AUPRC
            - Plot aggregated precision recall for 1st, 2nd pathway
            - Then plot for 1st + 2nd + 3rd
            - Etc.

        Note: this method will have to be called subsequent to the calculation
        of AUPRC per pathway.
        '''

        print("Aggregating successively by AUPRC per pathway")

        # First, read in the AUPRC scores per pathway

        # Initialize map of algorithms to lists of points
        auprc_algorithm_map = {}

        median_map = {}
        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            median_map[name] = []

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)
                auprc_algorithm_map[tup] = []

        # Now read in the fold AUPRC per pathway

        pairs = zip(self.collection.pathways, self.test_folds)

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:

                for j, fold in enumerate(test_fold):
                    
                    # Where we wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))


                    auprc_file = Path(
                        eval_dir, 
                        algorithm.get_output_directory(),
                        "auprc-nodes.txt") 

                    auprc = None

                    with auprc_file.open('r') as f:
                        line = next(f)
                        auprc = float(line.strip())

                    alg_name = algorithm.get_descriptive_name()

                    tup = (alg_name, pathway.name)

                    auprc_algorithm_map[tup].append(auprc)

        # Now make a sorted list of pathways by median AUPRC per algorithm

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            for pathway in self.collection.pathways:
                tup = (name, pathway.name)

                median = np.median(auprc_algorithm_map[tup])
                median_map[name].append((pathway.name, median))

        
        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            median_map[name] = sorted(
                median_map[name], 
                key=lambda x: x[1],
                reverse=True)

            print(median_map[name])

        # Now I have exactly what I need to determine which pathways to do
        # things with

        alg_to_groupings = {}

        for algorithm in self.algorithms:
            alg_name = algorithm.get_descriptive_name()
            
            ordered_pathways = median_map[alg_name]
            
            # Now I can get a list of lists of (name, median) tuples...

            # Given: [1 2 3 4]
            # I want: [[1], [1, 2], [1, 2, 3], [1, 2, 3, 4]]
            successive_groupings = [ordered_pathways[:i] for i in range(1, len(ordered_pathways) + 1)]

            alg_to_groupings[alg_name] = successive_groupings

            name_file = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name + "-nodes",
                "pathway-order-node-auc.txt")

            name_file.parent.mkdir(exist_ok=True, parents=True)

            with name_file.open('w') as f:
                last_grouping = successive_groupings[-1]
                names = [tup[0] for tup in last_grouping]
                for name in names:
                    f.write(name)
                    f.write("\n")

            fig2, ax2 = precrec.init_precision_recall_figure()

            ax2.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            # PDF file we will write
            overall_vis = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                alg_name + "-nodes",
                "overall.pdf")

            pairs = list(zip(
                self.collection.pathways, self.test_folds))

            for i, grouping in enumerate(successive_groupings):
                names = [tup[0] for tup in grouping]
                print(names)
    

                # pair_subset = [tup for tup in pairs if tup[0].name in names]

                # No... this will use the order of pathways in pairs

                # Instead, I need to somehow order the things by their position
                # in the names array. There are lots of easy ways to do this!

                # Heck, let's do it the dumb way

                pair_subset = []

                for name in names:
                    # This list comprehension should only return a list with
                    # a single item. Poor man's filter to preserve order of
                    # items in the names list
                    pair_subset = pair_subset + \
                        [tup for tup in pairs if tup[0].name == name]


                # Where we will write precision/recall, aggregated over
                # pathways
                pathway_collection_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.collection.name,
                    "aggregate-nodes",
                    self.get_output_prefix(),
                    alg_name)

                node_curves = []

                # I probably don't actually need the test_fold part of this
                # tuple
                for pathway, test_fold in pair_subset:

                    # Where we wrote precision/recall, aggregated over
                    # all folds per pathway
                    pathway_pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "aggregate")

                    pathway_pr_outfile_nodes = Path(
                        pathway_pr_output_dir, 
                        algorithm.get_output_directory(),
                        "precision-recall-nodes.txt") 

                    with pathway_pr_outfile_nodes.open('r') as f:
                        curve = precrec.read_precision_recall_fractions(f)
                        node_curves.append(curve)

                '''
                aggregated_nodes = precrec.aggregate_precision_recall_curve_fractions(
                    node_curves)
                '''
                aggregated_nodes = precrec.aggregate_precision_recall_curve_fractions_stretch(
                    node_curves)

                # Write aggregated curve back out
                pathway_collection_pr_outfile = Path(
                    pathway_collection_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall-nodes.txt") 

                pathway_collection_pr_outfile.parent.mkdir(
                    parents=True, exist_ok=True)

                with pathway_collection_pr_outfile.open("w") as f: 
                    precrec.write_precision_recall_fractions(f, aggregated_nodes)

                ######################## Plot node pr/
                fig, ax = precrec.init_precision_recall_figure()

                ax.set_title(
                    self.interactome.name + " " + self.collection.name 
                    + " "
                    + "\n" 
                    + alg_name
                    + self.get_details())

                # PDF file we will write
                vis_file_pdf = Path(
                    visualization_dir,
                    self.interactome.name,
                    self.collection.name,
                    self.get_output_prefix(),
                    alg_name + "-nodes",
                    "precision-recall-top-" + str(i) + "-nodes.pdf")
                
                vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

                precrec.plot_precision_recall_curve_fractions(
                    aggregated_nodes, label=algorithm.get_descriptive_name(), ax=ax)

                handles, labels = ax.get_legend_handles_labels()

                lgd = ax.legend(handles, labels, loc='upper center', 
                    bbox_to_anchor=(0.5,-0.1))

                fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')

                precrec.plot_precision_recall_curve_fractions(
                    aggregated_nodes, label=str(i), ax=ax2)
            
            overall_vis.parent.mkdir(parents=True, exist_ok=True)

            handles, labels = ax2.get_legend_handles_labels()

            lgd = ax2.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig2.savefig(str(overall_vis), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')
       

        # THIS PLOTS ALL ALGORITHMS ON THEIR BEST, THEN BEST + SECOND BEST,
        # ETC.

        # There are 15 pathways, so 15 groupings
        for grouping in range(len(self.collection.pathways)):

            # First, initialize a figure for the grouping. We want 15 figures
            # when we are done
        
            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " + self.collection.name 
                + " "
                + "\n" 
                + alg_name
                + self.get_details())

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.collection.name,
                self.get_output_prefix(),
                "aggregated-nodes",
                "precision-recall-top-" + str(grouping) + ".pdf")
            
            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                alg_name = algorithm.get_descriptive_name()

                pairs = list(zip(
                    self.collection.pathways, self.test_folds))

                group = alg_to_groupings[alg_name][grouping] 
                names = [tup[0] for tup in group]
                pair_subset = []

                for name in names:
                    # This list comprehension should only return a list with
                    # a single item. Poor man's filter to preserve order of
                    # items in the names list
                    pair_subset = pair_subset + \
                        [tup for tup in pairs if tup[0].name == name]

                edge_curves = []

                for pathway, _ in pair_subset:

                    # Where we wrote precision/recall, aggregated over
                    # all folds per pathway
                    pathway_pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "aggregate")

                    pathway_pr_outfile = Path(
                        pathway_pr_output_dir, 
                        algorithm.get_output_directory(),
                        "precision-recall-nodes.txt") 

                    with pathway_pr_outfile.open('r') as f:
                        curve = precrec.read_precision_recall_fractions(f)
                        edge_curves.append(curve)

                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions(
                    edge_curves)
                '''
                aggregated_edges = precrec.aggregate_precision_recall_curve_fractions_stretch(
                    edge_curves)

                precrec.plot_precision_recall_curve_fractions(
                    aggregated_edges, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')


    # TODO: I don't think it's worth it to split this function into
    # "evaluation" and "plotting" components right this very moment.
    def calculate_and_plot_precision_per_rank(
            self, reconstruction_dir=Path(), evaluation_dir=Path(),
            visualization_dir="outputs/pathway-reconstruction/visualization"):
        '''
        Look at each rank and plot the precision of that rank
        '''

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.collection.pathways],
            fold_creators)

        for i, (pathway, fc) in enumerate(creator_pathway_pairs):
            print(pathway.name)

            test_folds = self.test_folds[i] #fc.get_test_folds()

            for algorithm in self.algorithms:
                # Initialize figure
                fig, ax = precrec.init_precision_recall_figure()

                fig_just_one, ax_just_one = \
                    precrec.init_precision_recall_figure()

                fig_greater_than_half, ax_greater_than_half = \
                    precrec.init_precision_recall_figure()

                ax.set_xlim(auto=True)
                ax_just_one.set_xlim(auto=True)
                ax_greater_than_half.set_xlim(auto=True)

                # PDF file we will write
                vis_file_pdf = Path(
                    visualization_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    algorithm.get_descriptive_name() + 
                        "-precision-per-rank.pdf")

                # PDF file we will write
                vis_file_pdf_just_one = Path(
                    visualization_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    algorithm.get_descriptive_name() + 
                        "-precision-per-rank-just-one.pdf")

                # PDF file we will write
                vis_file_pdf_greater_than_half = Path(
                    visualization_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    algorithm.get_descriptive_name() + 
                        "-precision-per-rank-greater-than-half.pdf")
        
                # PNG file we will write 
                vis_file_png = Path(
                    visualization_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    algorithm.get_descriptive_name() + 
                        "-precision-per-rank.png")

                vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

                # Calculate points
                print("    " + algorithm.get_descriptive_name())
                predictions = []
                test_positives = []
                test_negatives = []

                avg_prec = []

                for j, fold in enumerate(test_folds):
                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "fold-" + str(j))
                    
                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        print("ALERT: RECONSTRUCTION FILE NOT FOUND")
                        reconstruction_file.touch()

                    # Where we will write precision/recall results
                    pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "aggregate")

                    positives = fold[0]
                    negatives = fold[1]

                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
                        predictions.append(fold_predictions)

                    test_positives.append(positives)
                    test_negatives.append(negatives)
                
                print("flattening")
                flat_test_pos = set(flatten_fold_aggregate(test_positives))
                flat_test_neg = set(flatten_fold_aggregate(test_negatives))
                flat_pred = flatten_fold_predictions(predictions)
                print("done")

                print("computing precision per rank")
                #pr_points = \
                #    precrec.compute_precision_recall_curve_negatives_decimals(
                #        flate_pred, flat_test_pos, flat_test_neg)

                points = \
                    precrec.compute_precision_per_rank_negatives(
                        flat_pred, flat_test_pos, flat_test_neg)

                points_just_one = \
                    precrec.compute_precision_per_rank_negatives_just_one(
                        flat_pred, flat_test_pos, flat_test_neg)

                points_greater_than_half = \
                    precrec.compute_precision_per_rank_negatives_greater_than_half(
                        flat_pred, flat_test_pos, flat_test_neg)
                print("done")

                def grouper(n, iterable, fillvalue=None):
                    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
                    args = [iter(iterable)] * n
                    return itertools.zip_longest(*args, fillvalue=fillvalue)

                # We'll do 20 groups to start
                #group_size = int(len(points) / 20)
                group_size = 10
                points = grouper(group_size, points)
                points = [[y for y in x if y != None] for x in points]
                points = [np.average(ls) for ls in points]

                points_just_one = grouper(group_size, points_just_one)
                points_just_one = [[y for y in x if y != None] for x in points_just_one]
                points_just_one = [np.average(ls) for ls in points_just_one]

                points_greater_than_half = grouper(group_size, points_greater_than_half)
                points_greater_than_half = [[y for y in x if y != None] for x in points_greater_than_half]
                points_greater_than_half = [np.average(ls) for ls in points_greater_than_half]

                points = points[:1000]
                points_just_one = points_just_one[:1000]
                points_greater_than_half = points_greater_than_half[:1000]

                ###############################################################
                ax.set_title(
                    self.interactome.name + " " +
                    self.collection.name + " " +
                    pathway.name + " " + 
                    self.get_details() + " " + 
                    "Group size: " +  str(group_size))

                ax.plot(points, label=algorithm.get_descriptive_name())

                handles, labels = ax.get_legend_handles_labels()

                lgd = ax.legend(handles, labels, loc='upper center', 
                    bbox_to_anchor=(0.5, -0.1))

                print("saving", (str(vis_file_pdf)))

                fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')

                fig.savefig(str(vis_file_png), bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')
            
                ###############################################################
                ax_just_one.set_title(
                    self.interactome.name + " " +
                    self.collection.name + " " +
                    pathway.name + " " + 
                    self.get_details() + " " + 
                    "Group size: " +  str(group_size))

                ax_just_one.plot(
                    points_just_one, label=algorithm.get_descriptive_name())

                handles, labels = ax_just_one.get_legend_handles_labels()

                lgd = ax_just_one.legend(handles, labels, loc='upper center', 
                    bbox_to_anchor=(0.5, -0.1))

                print("saving", (str(vis_file_pdf_just_one)))

                fig_just_one.savefig(str(vis_file_pdf_just_one), 
                    bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')

                ###############################################################
                ax_greater_than_half.set_title(
                    self.interactome.name + " " +
                    self.collection.name + " " +
                    pathway.name + " " + 
                    self.get_details() + " " + 
                    "Group size: " +  str(group_size))

                ax_greater_than_half.plot(points_greater_than_half, 
                    label=algorithm.get_descriptive_name())

                handles, labels = \
                    ax_greater_than_half.get_legend_handles_labels()

                lgd = ax_greater_than_half.legend(
                    handles, 
                    labels, 
                    loc='upper center', 
                    bbox_to_anchor=(0.5, -0.1))

                print("saving", (str(vis_file_pdf_greater_than_half)))

                fig_greater_than_half.savefig(
                    str(vis_file_pdf_greater_than_half), 
                    bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')






    def purge_results(self, reconstruction_dir=Path()):
        '''
        Delete previously-computed pathway reconstructions 
        for the algorithms specified in the config file.
        '''

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            training_folds = fc.get_training_folds()

            for j, fold in enumerate(training_folds):
                output_dir = Path(
                    reconstruction_dir,
                    self.interactome.name,
                    self.collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "fold-" + str(j))

                for algorithm in self.algorithms:
                    alg_dir = algorithm.get_full_output_directory(output_dir)
                    print(str(alg_dir))
                    if os.path.exists(str(alg_dir)):
                        shutil.rmtree(str(alg_dir))
    

    def run(self, *args, **kwargs): 
        '''
        0) Remove reconstructions created during previous runs of algorithms
           (this does not remove evaluations or plots at this point)

        1) Run each algorithm over each pathway in the pathway collection for
           each fold defined by the algorithm evaluator's get_fold_creator
           method

        2) Run evaluations (like precision/recall, AUPRC) for each run

        3) Plot the results of the above evaluations
        '''

        output_dir = kwargs.pop('output_dir')
        output_dir = Path(output_dir, "pathway-reconstruction")

        # TODO: Add as paramaters, and override with config-file specified 
        # directories in the pipeline itself
        reconstruction_dir = Path(output_dir, "reconstruction")
        evaluation_dir = Path(output_dir, "evaluation")
        visualization_dir = Path(output_dir, "visualization")

        print("Beginning evaluation of:\n"
            + "    interactome: %s\n" % self.interactome.name
            + "    pathway collection: %s\n" % self.collection.name
            + "    procedure: %s\n" % self.get_name())

        print("Running reconstructions...")
        self.run_reconstructions(reconstruction_dir)
        print("Finished running reconstructions!")
        
        print("Evaluating reconstructions...")
        self.evaluate_reconstructions(reconstruction_dir, evaluation_dir)
        print("Finished evaluating")

        print("Plotting results...")
        self.plot_results(evaluation_dir, visualization_dir)
        print("Finished plotting")
