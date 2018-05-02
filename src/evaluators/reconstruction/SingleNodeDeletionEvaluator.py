class SingleNodeDeletionEvaluator(AlgorithmEvaluator):
    '''
    1) Rank pathway nodes by centrality

    2) Delete the one with highest centrality, then replace it
       and delete the one with the second-highest, etc.
    
        3) After each deletion, run QuickLinker 
    
    4) After everything is done running, evaluate AUPRC for each "fold"
    
    5) Plot fold # vs. AUC. If hypothesis is correct, line should go
       up and to the right
    '''

    # To evaluate the results we build a chart, per pathway
    # 1) x-axis: degree, or order of node deleted?
    #    Okay, I can have nodes with the same degree have error bars
    # 2) y-axis: AUPRC for that point

    def get_fold_creator(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = SingleNodeDeletionFoldCreator(
            self.interactome, pathway, self.options)

        return fc

    # TODO: Apparently I'm not using these functions elsewhere? I can
    # probably ignore it here then, too...
    def get_name(self):
        return "single-node-deletion-eval"


    def get_output_prefix(self):
        return Path("single-node-deletion")


    def evaluate_reconstructions(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):

        self.calculate_auc_per_fold(reconstruction_dir, evaluation_dir)

        self.calculate_and_plot_lineplot(reconstruction_dir, evaluation_dir,
            Path(evaluation_dir.parent, "visualization"))


    # TODO: It looks like I've written this method generally enough to abstract
    # it out of each evaluator that uses this
    def calculate_auc_per_fold(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Calculate the precision recall curve and average precison
        for each fold independently, writing it to disk.
        '''

        print("----------------------------------------------------")
        print("Calculating area under the PR curve for each fold")
        print("----------------------------------------------------")

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            test_folds = fc.get_test_folds()
            for algorithm in self.algorithms:
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
                        reconstruction_file.touch()
        
                    # Where we will write precision/recall results
                    pr_output_dir = Path(
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

                    points = \
                        precrec.compute_precision_recall_curve_negatives_decimals(
                            fold_predictions, positives, negatives)
                    
                    auc = precrec.compute_area_under_precision_recall_curve(
                        points)

                    new_outfile = Path(
                        pr_output_dir, 
                        algorithm.get_output_directory(),
                        "auprc.txt") 

                    new_outfile.parent.mkdir(parents=True, exist_ok=True)

                    with new_outfile.open("w") as f: 
                        f.write(str(auc))
   

    def calculate_and_plot_lineplot(
            self, reconstruction_dir=Path(), evaluation_dir=Path(),
            visualization_dir=Path()):
        
        # TODO:
        # 1) Per pathway..
        #   2) Per algorithm...
        #      3) Per fold...
        #      - Read in the AUPRC for each fold
        #      - Place tue AUPRCs in bins according to the degree of the node 
        #        (degree in the pathway, not the interactome)
        #
        #            
        print("----------------------------------------------------")
        print("Plotting Node Degree vs. AUPRC")
        print("----------------------------------------------------")

        # {(pathway,algorithm) : {degree:[AUPRC scores]}}
        pathway_algorithm_map = {}

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = list(zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators))
        
        # Initialize the map with empty list
        for pathway, _ in creator_pathway_pairs:
            for algorithm in self.algorithms:
                name = algorithm.get_descriptive_name()
                pathway_algorithm_map[(pathway.name, name)] = {}

        print("----------------------------------------------------")
        print("First, read in AUPRCs")
        print("----------------------------------------------------")

        for pathway, fc in creator_pathway_pairs:
            test_folds = fc.get_test_folds()
            for algorithm in self.algorithms:
                for i, fold in enumerate(test_folds):

                    # fold[0] is set of test positives 
                    degree = len(set(fold[0]))

                    # Already-written AUPRC 
                    auprc_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    auprc_file= Path(
                        auprc_dir, 
                        algorithm.get_output_directory(),
                        "auprc.txt") 

                    point = None

                    with auprc_file.open('r') as f:
                        line = next(f)
                        point = float(line.strip())

                    tup = (pathway.name, algorithm.get_descriptive_name())

                    # {(pathway,algorithm) : {degree:[AUPRC scores]}}
                    degree_list = pathway_algorithm_map[tup].get(degree, [])

                    degree_list.append(point)

                    pathway_algorithm_map[tup][degree] = degree_list

        print("----------------------------------------------------")
        print("Second, plot AUPRCs")
        print("----------------------------------------------------")

        for pathway, _ in creator_pathway_pairs:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                "AUPRC vs. Degree of Node Deleted (%s)" % pathway.name
                + self.interactome.name + " "
                + self.pathway_collection.name)

            ax.set_xlabel("Degree of Node Deleted")
            ax.set_xlim(auto=True)

            ax.set_ylabel("AUPRC")

            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "auc.png")

            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "auc.pdf")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                alg_name = algorithm.get_descriptive_name()

                tup = (pathway.name, alg_name)

                x = sorted(list(pathway_algorithm_map[tup].keys()))

                y = []  

                stddev = []  

                for key in x:
                    # List of values associated with the degree (key)
                    vals = pathway_algorithm_map[tup].get(key,[])
                    y.append(sum(vals) / len(vals))
                    stddev.append(np.std(vals))

                print(x)
                print(y)
                print(stddev)
                ax.errorbar(x, y, yerr=stddev,label=alg_name)

            fig.savefig(str(vis_file_pdf), bbox_inches='tight')
            fig.savefig(str(vis_file_png), bbox_inches='tight')


    def plot_results(
            self, evaluation_dir=Path(), visualization_dir=Path()):
        # If I plot in the above function (hackier but whatever) this
        # function can do nothing
        None
