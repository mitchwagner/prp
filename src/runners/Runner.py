from pathlib import Path

class Runner(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''
    
    def get_name(self):
        '''
        Return a name for the evaluation performed by the Evaluator
        '''
        raise NotImplementedError()


    def run(self, *args, **kwargs):
        '''
        Glues the entire Evaluator's suite of functionality into a
        single runnable function
        '''
        should_run = kwargs.pop('should_run')

        if should_run:
            self.run_internal(*args, **kwargs)
                

    def run_internal(self, *args, **kwargs):
        raise NotImplementedError()
