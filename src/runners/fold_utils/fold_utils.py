'''
A module for reading/writing folds to disk and representing them in
memory.
'''

class EdgeTypeSeparatedFold:
    '''
    This class abstracts the notion of a machine learning "fold",
    keeping track of training and test positive negatives. In particular,
    this class is meant for folds of edges in a graph, where some edges
    are directed and others are undirected. 
    '''

    def  __init__(self, 
            train_dir_pos, train_undir_pos, train_dir_neg, train_undir_neg):
            test_dir_pos, test_undir_pos, test_dir_neg, test_undir_neg):
        '''
        Set instance variables.

        Edges are either in the train/test sets, are either
        directed/undirected, and are either positives or negatives. The
        algorithm constructor requires all 8 sets of edges to be specified.
        '''

        self.train_dir_pos = list(train_dir_pos)
        self.train_undir_pos = list(train_undir_pos)
        self.train_dir_neg = list(train_dir_neg)
        self.train_undir_neg = list(train_undir_neg)

        self.test_dir_pos = list(test_dir_pos)
        self.test_undir_pos = list(test_undir_pos)
        self.test_dir_neg = list(test_dir_neg)
        self.test_undir_neg = list(test_undir_neg)


    def write_to_disk(self, path):
        None


    def read_from_disk(path):
        None


    @property
    def train_dir_pos(self):
        return self.train_dir_pos


    @property
    def train_undir_pos(self):
        return self.train_dir_pos

    
    @property
    def train_dir_neg(self):
        return self.train_dir_pos

   
   @property
    def test_undir_pos(self):
        return self.test_dir_pos


    @property
    def test_dir_pos(self):
        return self.test_dir_pos


    @property
    def test_undir_pos(self):
        return self.test_dir_pos

    
    @property
    def test_dir_neg(self):
        return self.test_dir_pos

   
   @property
    def test_undir_pos(self):
        return self.test_dir_pos
