from pathlib import Path

class Evaluator(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''
    
    def get_name(self):
        '''
        Return a name for the evaluation performed by the Evaluator
        '''
        raise NotImplementedError()


    def run(self, output_dir=Path(), purge_results=False):
        '''
        Glues the entire Evaluator's suite of functionality into a
        single runnable function
        '''
        raise NotImplementedError()
