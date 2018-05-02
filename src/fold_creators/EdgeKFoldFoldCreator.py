from pathlib import Path

from sklearn.model_selection import KFold

import src.fold_creators.FoldCreator as fc

class EdgeKFoldFoldCreator(fc.FoldCreator):
    '''
    Create folds by dividing edges into folds.
    '''
    
    def __init__(self, interactome, pathway, options):
        self.interactome = interactome
        self.pathway = pathway
        self.num_folds = options["num_folds"]


    def create_positive_folds(self):
        '''
        Divide the positives into folds of train and test sets
        '''
        pathway_obj = self.pathway.get_pathway_obj()

        edges = fc.get_filtered_pathway_edges(pathway_obj, self.interactome)
        edges.sort(key=lambda edge:(edge[0], edge[1]))

        return split_items_into_folds(edges, self.num_folds)


    def create_negative_folds(self): 
        '''
        Divide the negatives into folds of train and test sets
        '''
        interactome_edges = set((x, y) 
            for x, y, line in self.interactome.get_interactome_edges())

        pathway_edges = self.pathway.get_pathway_obj().get_edges(data=False)
        pathway_edges = set(pathway_edges)

        negatives = list(interactome_edges.difference(pathway_edges)) 
        negatives.sort(key = lambda edge:(edge[0], edge[1]))

        return split_items_into_folds(negatives, self.num_folds)
        
        
    def get_fold_prefix(self, fold):
        return Path("fold-%d" % fold)
    

    def get_training_folds(self):
        '''
        Returns a list of tuples:
            (traing_positives, train_negatives, fold_name)
        '''
        positive_folds = self.create_positive_folds()
        negative_folds = self.create_negative_folds()

        folds = []

        for i, pair in enumerate(zip(positive_folds, negative_folds)):
            fold_name = self.get_fold_prefix(i)
            folds.append((pair[0][0], pair[1][0], fold_name))

        return folds


    def get_test_folds(self):
        '''
        Returns a list of tuples:
            (test_positives, test_negatives, fold_name)
        '''
        positive_folds = self.create_positive_folds()
        negative_folds = self.create_negative_folds()

        folds = []

        for i, pair in enumerate(zip(positive_folds, negative_folds)):
            fold_name = self.get_fold_prefix(i)
            folds.append((pair[0][1], pair[1][1], fold_name))

        return folds


def get_folds_from_split(items, split):
    """
    Scikit-learn returns a "split" structure that stores the indices
    of items in the train and test set of particular fold. This 
    takes that structure and the items that were divided up, to 
    return a structure of items instead of indices.
    """
    folds = []

    for i, (train, test) in enumerate(split):
        train_items = [items[x] for x in train]
        test_items = [items[y] for y in test]

        folds.append((train_items, test_items))
    return folds


def split_items_into_folds(items, num_folds):
    """
    Use Scikit-learn k-fold cross validation functions to divide
    the items supplied to the function into num_folds folds of 
    train and test sets.
    """
    kf = KFold(n_splits=num_folds, shuffle=True, random_state=1800)

    split = kf.split(items)
    
    folds = []

    return get_folds_from_split(items, split) 
