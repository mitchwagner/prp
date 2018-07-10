import time
from pathlib import Path

import scipy as sp
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors

import src.external.pathlinker.parse as pl_parse
import src.external.utils.precision_recall.precision_recall as precrec

from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph 

import src.algorithms.RankingAlgorithm as RankingAlgorithm

from src.evaluators.Evaluator import Evaluator

# TODO: In the past, I had this "write tp_fp with folds" function to 
# write out whether or not an edge was a true positive or a false positive

# I also plotted out the weights of these edges and divided them into two
# bins: true positives or false positives: to see what their overall weights
# are

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
        #group = groups.get(x[0][1], set())
        #group.add(x)

    regrouped = [] 
    for k, v in sorted(groups.items(), reverse=True):
        regrouped.append(v)

    final = [{(x[0][0], x[1]) for x in xs} for xs in regrouped]

    return final


class AlgorithmEvaluator(Evaluator):
    '''
    Base class for an object that has a "fold" creation procedure and
    runs a set of supplied algorithms over each fold.
    '''

    def __init__(
            self, interactome, pathway_collection, algorithms, options={},
                graphspace_settings=None):
        '''
        :param interactome: on-disk interactome object
        :param pathway_collection: PathwayCollection object
        :param algorithms: list of RankingAlgorithms
        :param options: map of options for the evaluator
        '''
        self.interactome = interactome
        self.pathway_collection = pathway_collection
        self.algorithms = algorithms
        self.options = options
        
        print("initializing folds")
        
        self.training_folds = self.get_training_folds()
        self.test_folds = self.get_test_folds()

        #if graphspace_settings != None:
        #    self.post_folds_to_graphspace(graphspace_settings)

        print("done initializing folds")


    def get_graphspace_name(self, pathway, fold):
        return str(self.get_output_prefix()) + "-" +  self.interactome.name \
            + "-" + pathway.name + "-" + str(fold)
        

    def post_folds_to_graphspace(self, graphspace_settings):
        '''
        Iterate through all pathways/folds and post the splits to GraphSpace
        '''

        print("Posting folds to GraphSpace")

        graphspace_instance = GraphSpace(
            graphspace_settings.email,
            graphspace_settings.password)


        triples = zip(self.pathway_collection.pathways, 
            self.training_folds, self.test_folds)

        for pathway, training_folds, test_folds in triples:

            pathway_obj = pathway.get_pathway_obj()

            sources = set(pathway_obj.get_receptors(data=False))
            targets = set(pathway_obj.get_tfs(data=False))

            for i, pair in enumerate(zip(training_folds, test_folds)):

                G = GSGraph()
                name = self.get_graphspace_name(pathway, i)
                G.set_name(name)

                pos_training = set(pair[0][0])
                pos_test = set(pair[1][0])
                neg_training = set(pair[0][1])
                neg_test = set(pair[1][1])

                all_edges = pos_training.union(
                    pos_test).union(neg_training).union(neg_test)
                
                nodes = set((node for edge in all_edges for node in edge))

                for node in nodes:

                    G.add_node(node, label=node)

                    if node in sources:
                        G.add_node_style(node, color='red', shape="triangle")

                    elif node in targets:
                        G.add_node_style(node, color='blue', shape="rectangle")

                    else:
                        G.add_node_style(node, color='blue', shape="ellipse")

                for edge in pos_training:
                    G.add_edge(edge[0], edge[1])
                    G.add_edge_style(edge[0], edge[1], directed=True, color='#FFDAB9')

                for edge in pos_test:
                    G.add_edge(edge[0], edge[1])
                    G.add_edge_style(edge[0], edge[1], directed=True, color='red')

                for edge in neg_training:
                    G.add_edge(edge[0], edge[1])
                    G.add_edge_style(edge[0], edge[1], directed=True, color='black')

                for edge in neg_test:
                    G.add_edge(edge[0], edge[1])
                    G.add_edge_style(edge[0], edge[1], directed=True, color='blue')
                    
                print("POSTING")
                graph_id = graphspace_instance.get_graph(name)

                if graph_id == None:    
                    graphspace_instance.post_graph(G)
                else:
                    graphspace_instance.update_graph(G)


    def get_name(self):
        raise NotImplementedError()


    def get_details(self):
        raise NotImplementedError()


    def get_output_prefix(self):
        raise NotImplementedError()


    def get_fold_creator(self, pathway):
        raise NotImplementedError()


    def get_fold_creators(self):
        fcs = []

        for pathway in self.pathway_collection.pathways:
            fcs.append(self.get_fold_creator(pathway))

        return fcs


    def get_training_folds(self):
        '''
        Returns a list of tuples (positives, negatives) detailing the 
        set of training positives/negatives for each fold
        '''
        fold_creators = self.get_fold_creators()
        folds = []

        for fc in fold_creators:
            training_folds = fc.get_training_folds()
            folds.append(training_folds)

        return folds


    def get_test_folds(self):
        '''
        Returns a list of tuples (positives, negatives) detailing the
        set of test positives/negatives for each fold
        '''
        fold_creators = self.get_fold_creators()
        folds = []

        for fc in fold_creators:
            test_folds = fc.get_test_folds()
            folds.append(test_folds)

        return folds


    def run_reconstructions(self, output_dir=Path()):
        '''
        Run each algorithm over each fold of each pathway in the 
        pathway collection.
        '''

        pairs = zip(self.pathway_collection.pathways, self.training_folds)

        for pathway, training_fold in pairs: 
            for fold in training_fold:
                print(fold[2])
                for algorithm in self.algorithms:
                    # First, write output directory
                    full_output_dir = Path(
                        output_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    alg_dir = algorithm.get_full_output_directory(
                        full_output_dir)
                    
                    # TODO: Is this step even necessary?
                    alg_dir.mkdir(parents=True, exist_ok=True)
                    
                    # TODO: Leaving this here in case we want to revisit
                    # it later. This would limit the sources and targets
                    # provided to only those incident on edges in the set
                    # of training positives. Needs to be finished.
                    '''                    
                    ###########################################################
                    # Create a new pathway_nodes_file 

                    # 1) Get a list of training nodes
                    training_nodes = set()
                    
                    positive_edges = fold[0]

                    for edge in positive_edgeS: 
                        training_nodes.add(edge[0])
                        training_nodes.add(edge[1])

                    # 2) Define a new nodes file
                    new_nodes_file = Path(alg_dir, "new_nodes.txt")    

                    with pathway.get_nodes_file.open('r') as f1,
                            new_nodes_file.open('w') as f2:
                        for line in f1:
                            if line.startswith("#"):
                                f2.write(line)
                            else:
                                toks = line.split("\t")

                    ###########################################################
                    '''

                    # Second, run the algorithms         
                    alg_input = RankingAlgorithm.PathwayReconstructionInput(
                        self.interactome.path,
                        fold[0], # training positive edges
                        pathway.get_nodes_file(),
                        full_output_dir,
                        pathway.get_edges_file(), # for sanity checks
                        fold[1]) # training negative edges

                    print("Running algorithm: ")
                    print("    " + self.interactome.name)
                    print("    " + self.pathway_collection.name)
                    print("    " + pathway.name)
                    print("    " + self.get_name())
                    print("    " + str(fold[2]))
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
        
        # AUPRC, Average Precision, and Max F1-measure plots 
        self.calculate_metrics(reconstruction_dir, evaluation_dir)
        self.calculate_wilcoxon_scores(evaluation_dir)

        # Precision/Recall plots
        self.aggregate_pr_over_folds(reconstruction_dir, evaluation_dir)
        self.aggregate_pr_over_pathways(evaluation_dir)

        # Precision per rank 
        self.calculate_and_plot_precision_per_rank(
            reconstruction_dir, evaluation_dir)

        #self.s_t_paths_analysis(reconstruction_dir, evaluation_dir,
        #    Path(evaluation_dir.parent, "visualization"))


    def calculate_metrics(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Calculate the precision recall curve and average precison
        for each fold independently, writing it to disk.
        '''

        pairs = zip(self.pathway_collection.pathways, self.test_folds)

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                for fold in test_fold:
                    
                    # TODO: get a directory from a function instead

                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

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
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    positives = set(fold[0])
                    negatives = set(fold[1])

                    fold_predictions = None
                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
               
                    fold_predictions = [
                        set([tup[0] for tup in s])
                        for s in fold_predictions]

                    print("Calculating precision/recall points")
                    print("Pathway:", pathway.name)
                    print("Algorithm:", algorithm.get_descriptive_name())

                    print(positives)
                    print("--------")
                    print(negatives)

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
   

    def calculate_wilcoxon_scores(self, evaluation_dir=Path()):
        '''
        Write out a TSV relating Wilcoxon scores between algorithms
        '''
        pairs = zip(self.pathway_collection.pathways, self.test_folds)

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
                for fold in test_fold:
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

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
            self.pathway_collection.name,
            self.get_output_prefix(),
            "wilcoxon-avg-prec.tsv")

        wilcoxon_f1_max_file = Path(
            evaluation_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "wilcoxon-f1-max.tsv")

        wilcoxon_auprc_file = Path(
            evaluation_dir,
            self.interactome.name,
            self.pathway_collection.name,
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
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            print(pathway.name)
            test_folds = fc.get_test_folds()
            for algorithm in self.algorithms:
                print("    " + algorithm.get_descriptive_name())
                predictions = []
                test_positives = []
                test_negatives = []

                avg_prec = []

                # Where we will write precision/recall results
                pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.pathway_collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                for fold in test_folds:
                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])
                    
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

                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
                        predictions.append(fold_predictions)
                        
                    test_positives.append(positives)
                    test_negatives.append(negatives)



                    
                    # Write individual fold prec/rec curves
                    fold_output_file = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2],
                        algorithm.get_output_directory(),
                        "precision-recall.txt") 

                    just_edges = [set([x[0] for x in y]) 
                        for y in fold_predictions]
                    
                    fold_curve = \
                        precrec.compute_precision_recall_curve_negatives_fractions(
                            just_edges, set(positives), set(negatives))

                    fold_output_file.parent.mkdir(parents=True, exist_ok=True) 

                    with fold_output_file.open("w") as f: 
                        precrec.write_precision_recall_fractions(f, fold_curve)
                
                print("flattening")
                flat_test_pos = set(flatten_fold_aggregate(test_positives))
                flat_test_neg = set(flatten_fold_aggregate(test_negatives))
                flat_pred = flatten_fold_predictions(predictions)
                print("done")

                print("computing p/r")
                # Call existing precrec functions passing these things above
                points = \
                    precrec.compute_precision_recall_curve_negatives_fractions(
                        flat_pred, flat_test_pos, flat_test_neg)
                print("done")

                new_outfile = Path(
                    pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                new_outfile.parent.mkdir(parents=True, exist_ok=True)

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
            self.pathway_collection.name,
            "aggregate",
            self.get_output_prefix())

        for algorithm in self.algorithms:    
            curves = []
            
            # Where we wrote precision/recall, aggregated over
            # all folds per pathway
            for pathway in self.pathway_collection.pathways:
                pathway_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.pathway_collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "aggregate")

                pathway_pr_outfile = Path(
                    pathway_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                with pathway_pr_outfile.open('r') as f:
                    curve = precrec.read_precision_recall_fractions(f)
                    curves.append(curve)

            print("Aggregating precision/recall curves")
            aggregated = precrec.aggregate_precision_recall_curve_fractions(
                curves)

            # Write aggregated curve back out
            pathway_collection_pr_outfile = Path(
                pathway_collection_pr_output_dir, 
                algorithm.get_output_directory(),
                "precision-recall.txt") 

            pathway_collection_pr_outfile.parent.mkdir(
                parents=True, exist_ok=True)

            with pathway_collection_pr_outfile.open("w") as f: 
                precrec.write_precision_recall_fractions(f, aggregated)

        
    """

    # TODO: This needs to be updated to reflect Aditya's most recent changes!
    def s_t_paths_analysis(
            self, reconstruction_dir=Path(), evaluation_dir=Path(),
            visualization_dir=Path()):
        '''
        See how quickly the algorithms connect all receptors to all
        transcription factors.
        '''
        
        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:

            # Start with the p-graph
            train_folds = fc.get_training_folds()
            test_folds = fc.get_test_folds()

            fig, ax = plt.subplots()

            ax.set_title(
                "Fraction s-t pairs Connected vs. Rank"
                + self.interactome.name + " "
                + self.pathway_collection.name + " "
                + pathway.name + "\n"
                + self.get_details())

            ax.set_xlabel("Rank")
            ax.set_ylabel("Fraction s-t Pairs Connected")

            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "st-pairs.png")

            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "st-pairs.pdf")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:

                rank_st_path_points = [] 

                for (train_fold, test_fold) in zip(train_folds, test_folds):
                    points = []

                    # 1) Get the list of edges for the algorithm               

                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        train_fold[2])

                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        reconstruction_file.touch()
        
                    # 2) Read the list of edges in by rank
                    fold_predictions = None
                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)

                    # 3) Read the lists of sources and targets

                    pathway_obj = pathway.get_pathway_obj()

                    sources = pathway_obj.get_receptors(data=False)
                    targets = pathway_obj.get_tfs(data=False)
                    
                    # 4) Construct the base graph

                    base_graph = nx.DiGraph()

                    for source in sources:
                        base_graph.add_node(source)

                    for target in targets:
                        base_graph.add_node(target)

                    for edge in set(train_fold[0]):
                        base_graph.add_edge(edge[0], edge[1])

                    # 5) For reach rank...
                    for i, rank in enumerate(fold_predictions):

                        # 5) Add the edges with the rank to the graph
                        for edge in rank:
                            base_graph.add_edge(edge[0][0], edge[0][1])
                        if i%20 == 0:
                            # 6) For each source/target pair, try to find a
                            # path. Keep track of the number of successes

                            path_sum = 0
                            
                            for source in sources:
                                for target in targets:
                                    if nx.has_path(
                                        base_graph, source,target):
                                        path_sum += 1
                            
                            total_possible = len(sources) * len(targets)
                            frac = float(path_sum) / float(total_possible)
 
                            # 7) Create a tuple (rank, sum)
                            tup = (i, frac)

                        # 8) Append tuple to a list
                            points.append(tup)
                            #if frac == 1:
                    # 9) Append list to overall list
                    rank_st_path_points.append(points)

                # 10) Avg. across folds
                maxLen = max([len(x) for x in rank_st_path_points])
                AvgYs = [0] * maxLen
                Counts = [0] * maxLen
                StddevList = [[] for i in range(maxLen)]

                # Make all lists the same length by padding the end of the
                # shorter ones with their last element
                for x in rank_st_path_points:
                    xLen =  len(x)
                    lastElement = x[-1]
                    
                    for i in range(maxLen - xLen):
                        x.append(lastElement)
                
                for i, ls in enumerate(rank_st_path_points):

                    for j in range(len(ls)):
                        AvgYs[j] += ls[j][1]
                        StddevList[j].append(ls[j][1])
                        Counts[j] += 1
                    
                AvgYs = [x/Counts[i] for i,x in enumerate(AvgYs)]
                Ebars = [np.std(x) for i,x in enumerate(StddevList)]

                xs = [i*20 for i in range(len(AvgYs))]
                    
                label = algorithm.get_descriptive_name() 
                ax.plot(xs, AvgYs, label=label)
                Ebars_minus = [max(0,x-Ebars[i]) for i,x in enumerate(AvgYs)]
                Ebars_plus = [min(1,x+Ebars[i]) for i,x in enumerate(AvgYs)]
                ax.fill_between(xs, Ebars_minus, Ebars_plus,alpha=0.5,)

                ax.set_ybound(0,1)
                ax.set_xbound(0,1000)
            ax.legend(loc='best')

            fig.savefig(str(vis_file_pdf), bbox_inches='tight')
            fig.savefig(str(vis_file_png), bbox_inches='tight')

    """


    def plot_results(
            self, evaluation_dir=Path(), visualization_dir=Path()):
        '''
        Run all plotting algorithms
        '''

        # Boxplots
        self.all_algorithm_scores_combined_pathways_boxplot(
            evaluation_dir, visualization_dir)

        self.all_algorithm_scores_individual_pathways_boxplots( 
            evaluation_dir, visualization_dir)

        self.individual_algorithm_scores_all_pathways_boxplots(
            evaluation_dir, visualization_dir)

        # Precision/Recall
        self.plot_pr_individual_pathways(evaluation_dir, visualization_dir)
        self.plot_pr_all_pathways(evaluation_dir, visualization_dir)


    def all_algorithm_scores_combined_pathways_boxplot(
            self, evaluation_dir, visualization_dir): 
        '''
        Make a single boxplot figure, with one boxplot per algorithm, where
        a constituent point is an algorithm's score on a particular fold.
        '''

        pairs = zip(self.pathway_collection.pathways, self.test_folds)

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
                for fold in test_fold:
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

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
            + self.pathway_collection.name + "\n"
            + self.get_details())

        ax_avg_prec.set_xlabel("Average Precision")
        ax_avg_prec.set_ylabel("Algorithm")

        ax_f1_max.set_title(
            "Max F-Measure by Algorithm "
            + self.interactome.name + " "
            + self.pathway_collection.name + "\n"
            + self.get_details())

        ax_f1_max.set_xlabel("Max F-Measure")
        ax_f1_max.set_ylabel("Algorithm")

        ax_auprc.set_title(
            "AUPRC by Algorithm "
            + self.interactome.name + " "
            + self.pathway_collection.name + "\n"
            + self.get_details())

        ax_auprc.set_xlabel("AUPRC")
        ax_auprc.set_ylabel("Algorithm")

        avg_prec_png = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "avg_prec.png")

        avg_prec_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "avg_prec.pdf")

        f1_max_png = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "f1_max.png")

        f1_max_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "f1_max.pdf")

        auprc_png = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "auprc.png")

        auprc_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
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

        pairs = zip(self.pathway_collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        avg_prec_algorithm_map = {}
        f1_max_algorithm_map = {}
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            for pathway in self.pathway_collection.pathways:
                tup = (name, pathway.name)
                avg_prec_algorithm_map[tup] = []
                f1_max_algorithm_map[tup] = []
                auprc_algorithm_map[tup] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []
                for fold in test_fold:
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

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

        for pathway in self.pathway_collection.pathways:
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
                + self.pathway_collection.name + "\n"
                + self.get_details())

            ax_avg_prec.set_xlabel("Average Precision")
            ax_avg_prec.set_ylabel("Algorithm")

            ax_f1_max.set_title(
                "Max F-Measure by Algorithm (%s)" % pathway.name
                + self.interactome.name + " "
                + self.pathway_collection.name + "\n"
                + self.get_details())

            ax_f1_max.set_xlabel("Max F-Measure")
            ax_f1_max.set_ylabel("Algorithm")

            ax_auprc.set_title(
                "AUPRC by Algorithm (%s)" % pathway.name
                + self.interactome.name + " "
                + self.pathway_collection.name + "\n"
                + self.get_details())

            ax_auprc.set_xlabel("AUPRC")
            ax_auprc.set_ylabel("Algorithm")

            avg_prec_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "avg_prec.png")

            avg_prec_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "avg_prec.pdf")

            f1_max_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "f1_max.png")

            f1_max_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "f1_max.pdf")

            auprc_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "auprc.png")

            auprc_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
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

        pairs = zip(self.pathway_collection.pathways, self.test_folds)

        # First, read in the scores 

        # Map of algorithms to lists of points
        avg_prec_algorithm_map = {}
        f1_max_algorithm_map = {}
        auprc_algorithm_map = {}

        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()

            for pathway in self.pathway_collection.pathways:
                tup = (name, pathway.name)
                avg_prec_algorithm_map[tup] = []
                f1_max_algorithm_map[tup] = []
                auprc_algorithm_map[tup] = []

        for pathway, test_fold in pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []
                for fold in test_fold:
                    
                    # Where we will wrote scores to
                    eval_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

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

            for pathway in self.pathway_collection.pathways:
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
                + self.pathway_collection.name + "\n"
                + self.get_details())

            ax_avg_prec.set_xlabel("Average Precision")
            ax_avg_prec.set_ylabel("Pathway")

            ax_f1_max.set_title(
                "Max F-Measure by Pathway"
                + self.interactome.name + " "
                + self.pathway_collection.name + "\n"
                + self.get_details())

            ax_f1_max.set_xlabel("Max F-Measure")
            ax_f1_max.set_ylabel("Pathway")

            ax_auprc.set_title(
                "AUPRC by Pathway"
                + self.interactome.name + " "
                + self.pathway_collection.name + "\n"
                + self.get_details())

            ax_auprc.set_xlabel("AUPRC")
            ax_auprc.set_ylabel("Pathway")

            avg_prec_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                name + "-avg_prec.png")

            avg_prec_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                name + "-avg_prec.pdf")

            f1_max_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                name + "-f1_max.png")

            f1_max_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                name + "-f1_max.pdf")

            auprc_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                name + "-auprc.png")

            auprc_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
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

    """
        # Snippets of old Wilcoxon plotting code
        for i, algorithm1 in enumerate(self.algorithms):
            for algorithm2 in self.algorithms:

                # Greater and significant: green
                # Lesser and significant: red
                # Not significant: black
                # Colormap defined below
                if (p_val < corrected_alpha):
                    if median1 > median2: 
                        matrix[i].append(2) 
                    else:
                        matrix[i].append(1)
                else:
                    matrix[i].append(0)

        fig, ax = plt.subplots()

        ax.set_title(
            "Wilcoxon Rank Sum Test"
            + self.interactome.name + " "
            + self.pathway_collection.name + "\n"
            + "Node Percent Kept: " + str(
                self.options["percent_nodes_to_keep"]) + " "
            + "Edge Percent Kept: " + str(
                self.options["percent_edges_to_keep"]) + " " 
            + "Iterations: " + str(
                self.options["iterations"]))

        ax.set_xlabel("Algorithm")
        ax.set_ylabel("Algorithm")

        vis_file_png = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "wilcoxon.png")

        vis_file_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "wilcoxon.pdf")

        array = np.array(matrix)
        print(array)

        cmap = colors.ListedColormap([[0, 0, 0], [1, 0 ,0], [0, 1, 0]])
        ax.matshow(array, cmap=cmap)
        
        plt.xticks(range(0, len(array)), labels, rotation="vertical")
        ax.xaxis.tick_bottom()

        plt.yticks(range(0, len(array)), labels)

        #ax.set_xticklabels(['']+labels)
        #ax.set_yticklabels(['']+labels)
        
        vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)


        fig.savefig(str(vis_file_pdf), bbox_inches='tight')
        fig.savefig(str(vis_file_png), bbox_inches='tight')
    """


    def plot_pr_individual_pathways(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        for pathway in self.pathway_collection.pathways:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " +
                self.pathway_collection.name + " " +
                pathway.name)

            # Where we wrote precision/recall results
            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "aggregate")

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "precision-recall.pdf")
    
            # PNG file we will write 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
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

        for pathway in self.pathway_collection.pathways:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " +
                self.pathway_collection.name + " " +
                self.get_details())

            # Where we wrote precision/recall results
            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.pathway_collection.name,
                "aggregate",
                self.get_output_prefix())

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                "precision-recall.pdf")
    
            # PNG file we will write 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
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
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            print(pathway.name)
            test_folds = fc.get_test_folds()


            for algorithm in self.algorithms:
                # Initialize figure
                fig, ax = precrec.init_precision_recall_figure()

                ax.set_title(
                    self.interactome.name + " " +
                    self.pathway_collection.name + " " +
                    pathway.name + " " + 
                    self.get_details())

                ax.set_xlim(auto=True)

                # PDF file we will write
                vis_file_pdf = Path(
                    visualization_dir,
                    self.interactome.name,
                    self.pathway_collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    algorithm.get_descriptive_name() + 
                        "-precision-per-rank.pdf")
        
                # PNG file we will write 
                vis_file_png = Path(
                    visualization_dir,
                    self.interactome.name,
                    self.pathway_collection.name,
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

                for fold in test_folds:
                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])
                    
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
                        self.pathway_collection.name,
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
                points = \
                    precrec.compute_precision_per_rank_negatives(
                        flat_pred, flat_test_pos, flat_test_neg)
                print("done")

                ax.plot(points, label=algorithm.get_descriptive_name())

                handles, labels = ax.get_legend_handles_labels()

                lgd = ax.legend(handles, labels, loc='upper center', 
                    bbox_to_anchor=(0.5,-0.1))

                print("saving", (str(vis_file_pdf)))

                fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')

                fig.savefig(str(vis_file_png), bbox_extra_artists=(lgd,), 
                    bbox_inches='tight')


    def purge_results(self, reconstruction_dir=Path()):
        '''
        Delete previously-computed pathway reconstructions 
        for the algorithms specified in the config file.
        '''

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            training_folds = fc.get_training_folds()
            for fold in training_folds:
                output_dir = Path(
                    reconstruction_dir,
                    self.interactome.name,
                    self.pathway_collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    fold[2])

                for algorithm in self.algorithms:
                    alg_dir = algorithm.get_full_output_directory(output_dir)
                    print(str(alg_dir))
                    if os.path.exists(str(alg_dir)):
                        shutil.rmtree(str(alg_dir))
    

    def run(self, output_dir=Path(), purge_results=False):
        '''
        0) Remove reconstructions created during previous runs of algorithms
           (this does not remove evaluations or plots at this point)

        1) Run each algorithm over each pathway in the pathway collection for
           each fold defined by the algorithm evaluator's get_fold_creator
           method

        2) Run evaluations (like precision/recall, AUPRC) for each run

        3) Plot the results of the above evaluations
        '''

        output_dir = Path(output_dir, "pathway-reconstruction")

        # TODO: Add as paramaters, and override with config-file specified 
        # directories in the pipeline itself
        reconstruction_dir = Path(output_dir, "reconstruction")
        evaluation_dir = Path(output_dir, "evaluation")
        visualization_dir = Path(output_dir, "visualization")

        if purge_results:
            self.purge_results(reconstruction_dir)

        print("Beginning evaluation of:\n"
            + "    interactome: %s\n" % self.interactome.name
            + "    pathway collection: %s\n" % self.pathway_collection.name
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
