from pathlib import Path
import os

import src.external.pathlinker.parse as pl_parse

class PathwayReconstructionInput(object):
    """
    Input necessary for a RankingAlgorithm to create a reconstruction.
    """

    def __init__(self, interactome, training_edges, pathway_nodes_file, 
            output_dir, all_edges_file, training_negatives):
        self.interactome = interactome
        
        # The ENTIRE set of edges in a pathway (not just positives for a given
        # fold
        self.all_edges_file = all_edges_file

        # A pathway edge with JUST positives for a given fold
        # TODO: should probably rename, but most if not all algorithms have
        # hardcoded dependence on this parameter, so it would take a few 
        # minutes
        self.training_edges = training_edges 
        self.training_negatives = training_negatives

        self.pathway_nodes_file = pathway_nodes_file
        self.output_dir = output_dir


    def label_interactome_file(
            self, in_handle, out_handle, sets, default="none"):
        '''
        Sets is a list of tuples, where the first element in each
        tuple is a label, and the second is the set of elements tto 
        give that label.

        The sets of elements are intended to be disjoint. 

        :param default: the default label to give to anything that does
            not match anything in any of the provided sets
        '''

        edge_label_map = {}

        for s in sets:
            for edge in s[1]:
                # Don't overwrite earlier label
                if edge not in edge_label_map:
                    edge_label_map[edge] = s[0]

        for line in in_handle:
            if pl_parse.is_comment_line(line):
                out_handle.write(line)
            else:
                tokens = pl_parse.tokenize(line)
                edge = (tokens[0], tokens[1])

                label = edge_label_map.get(edge, default)

                '''
                flag = False
                label = ""
                for s in sets:
                    # If the edge is in this set, label it with the set's 
                    # label and stop checking sets
                    if edge in s[1]:
                        label = s[0]
                        flag = True
                        break;

                # If the edge was not in any set, give the label "none"
                if flag == False:
                    label = default
                '''

                out_handle.write(line.rstrip() + "\t" + label + "\n")


class RankingAlgorithm(object):
    """
    Abstract representation of a RankingAlgorithm in the pipeline.
    """

    def run_wrapper(self, reconstruction_input, should_force=False):
        """
        Call the RankingAlgorithm's run function, first performing 
        some bookkeeping, including making sure the output directory 
        exists, and checking to make sure if it should overwrite
        previously-written results 

        :param reconstruction_input: an object containing the input 
            parameters required by the reconstruction algorithm

        :param should_force: whether or not to force overwriting of 
            previous results
        """
        self.ensure_output_directory_exists(reconstruction_input)

        if (should_force):
            self.run(reconstruction_input)
        else:
            if not self.output_previously_written(reconstruction_input): 
                self.run(reconstruction_input)
            else:
                print("Skipping (already run):", self.get_descriptive_name())

        self.conform_output(reconstruction_input.output_dir)


    def run(self, reconstruction_input):
        raise NotImplementedError() 


    def output_previously_written(self, file_location_context):
        """
        Return a boolean inidicating if this algorithm has previously
        written its output to the location specified by the 
        file location object passed in.

        Can be overwritten to check more files, if necessary.
        """
        return os.path.exists(
            str(self.get_full_output_file(
                file_location_context.output_dir)))


    def conform_output(self, output_dir):
        raise NotImplementedError()


    def get_name(self):
        """
        Return the full name of the algorithm
        """
        raise NotImplementedError()


    def get_descriptive_name(self):
        return self.get_name()


    def get_output_file(self):
        return "ranked-edges.txt"


    def get_output_directory(self):
        """
        Return the name of the folder to store results in. By default, the
        name of the folder is just the name of the algorithm.
        """
        return self.get_name()


    def get_full_output_file(self, output_dir):
        return Path(
            self.get_full_output_directory(output_dir),
            self.get_output_file())


    def get_full_output_directory(self, output_dir):
        """
        Return the full output directory the RankingAlgorithm 
        should write the output of its run function to.
        """
        return Path(output_dir, self.get_output_directory())


    def ensure_output_directory_exists(self, reconstruction_input):
        outdir = self.get_full_output_directory(
            reconstruction_input.output_dir)

        outdir.mkdir(parents=True, exist_ok=True)
