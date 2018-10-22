import csv

from pathlib import Path

# TODO: This code is obsolete, as the import below has been removed.
# Leaving this here for now if anybody would like to adapt it in the 
# future.
from src.evaluators.Evaluator import Evaluator

class RemovalEvaluator(Evaluator):
    '''
    Analyze the impact of creating folds in a given way
    '''

    def __init__(
            self, interactome, pathway_collection, options={}):
        '''
        :param interactome: on-disk interactome object
        :param pathway_collection: PathwayCollection object
        :param options: map of options for the evaluator
        '''
        self.interactome = interactome
        self.pathway_collection = pathway_collection
        self.options = options


    def get_fold_creator(self, pathway):
        raise NotImplementedError()


    def analyze_fold_creation(self, pathway):
        ''' 
        Return a list of lists, where each interior list can be viewed
        as a column corresponding to statistics on a particular row
        '''
        raise NotImplementedError()


    def get_output_prefix(self):
        '''
        Return a directory name indicating the options passed in to
        parameterize the fold creation process.
        '''
        raise NotImplementedError()


    def get_outfile(self, output_dir, pathway):
        return Path(
            output_dir, self.get_output_prefix(), pathway.name + "-stats.txt")
   

    def average_stats(self, stats):
        '''
        This function averages statistics gathered for a particular pathway.
        By putting this function here, I'm putting particular constraints on
        all of the rest of the removal evaluators.

        In particular, I'm assuming that they all return stats with a 
        string in the first column, and a header row.
        '''
        avg_row = []

        for col in stats[0]:
            avg_row.append(0)
        
        # Preserve the pathway name, which can't be averaged
        avg_row[0] = stats[1][0] + "-average"

        for row in stats[1:]:
            # Skip the first column, as it is alphanumeric
            for i, col in enumerate(row[1:]):
                avg_row[i + 1] += col
        
        for i, col in enumerate(avg_row[1:]):
            avg_row[i+1] = col / (len(stats) - 1)

        stats.append(avg_row)


    def write_stats(self, stats, output_dir, pathway):
        outfile = self.get_outfile(output_dir, pathway)

        outfile.parent.mkdir(parents=True, exist_ok=True)

        with outfile.open('w') as f:
            writer = csv.writer(f, delimiter='\t') 

            for row in stats:
                writer.writerow(row)


    # TODO: the purge_results option is not implemented
    def run(self, output_dir=Path(), purge_results=False):
        output_dir = Path(output_dir, "fold_stats")

        for pathway in self.pathway_collection.pathways: 
            stats = self.analyze_fold_creation(pathway)
            self.average_stats(stats)
            self.write_stats(stats, output_dir, pathway)
