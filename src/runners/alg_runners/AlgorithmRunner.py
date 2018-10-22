import time
import itertools
from pathlib import Path

from src.runners.Runner import Runner

class AlgorithmRunner(Runner):
    '''
    Runs algorithms over... what? What should be the base level of
    running, here? Before, we ran it over each fold, for each 
    pathway in the pathway collection

    Really, I could just paramaterize this with a fold directory and
    a name for that fold.

    This means I'll need a way to read the fold in, but that doesn't
    have to be defined here.
    '''


    def __init__(self, *args, **kwargs):
        '''
        Initialization sets instance variables
        '''
        self.interactome = kwargs.pop('interactome')
        self.pathway = kwargs.pop('pathway_collection')
        self.algorithms = kwargs.pop('algorithms')
        
        self.output_path = kwargs.pop('output_path')
        self.fold_path = kwargs.pop('fold_path')
   

    def run_internal(self, *args, **kwargs):
        '''
        oooooh I can write names out to folds. yay!
        '''

        # read in the fold

        for algorithm in self.algorithms:
            # figure out where to write output and mkdir it
            # create alg_input
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
