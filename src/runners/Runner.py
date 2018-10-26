from pathlib import Path

class Runner(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''
    
    def get_name(self) -> str:
        '''
        Return a name for the evaluation performed by the Evaluator
        '''
        raise NotImplementedError()


    def run(self, *args, **kwargs) -> None:
        '''
        Wrapper that first checks if the Runner's task should in fact
        be performed.
        '''
        should_run = kwargs.pop('should_run')

        if should_run:
            self.run_internal(*args, **kwargs)
                

    def run_internal(self, *args, **kwargs) -> None:
        '''
        The task to be performed by the runner.
        '''
        raise NotImplementedError()
