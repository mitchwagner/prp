from pathlib import Path
from typing import List

# Local imports
from src.algorithms.RankingAlgorithm import RankingAlgorithm


class Runner(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''

    def __init__(self, *args, **kwargs):
        self.interactome = kwargs.pop('interactome')
        self.collection = kwargs.pop('collection')
        self.algorithms = kwargs.pop('algorithms')
        self.graphspace = kwargs.pop('graphspace')
        
    
    def get_name(self) -> str:
        '''
        Return a name for the evaluation performed by the Evaluator
        '''
        raise NotImplementedError()


    def run(self, *args, **kwargs) -> None:
        '''
        Task to be performed.
        '''
        raise NotImplementedError()
