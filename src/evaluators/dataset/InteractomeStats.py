from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np


from src.evaluators.Evaluator import Evaluator
import src.external.pathlinker.PathLinker as pl

class InteractomeStats(Evaluator):
    '''
    Runnable analysis that outputs several statistics about the 
    interactomes.
    '''


    def __init__(self, interactome):
        self.interactome = interactome 

    
    def get_name(self):
        '''
        Return a name for the evaluation performed by this Evaluator 
        '''
        return "Interactome Statistics Evaluation"


    def get_interactome_net(self):
        with self.interactome.path.open('r') as f:
            return pl.readNetworkFile(f)

    
    def get_output_dir(self):
        return Path("interactome-stats")


    def write_summary_statistics(self, output_dir):
        name = self.interactome.name
        net = self.get_interactome_net()

        outfile = Path(output_dir, self.get_output_dir(), name, "summary.txt")

        outfile.parent.mkdir(parents=True, exist_ok=True)

        with outfile.open('w') as f:
            f.write("Number of nodes: %d\n" % len(net.nodes())) 
            f.write("Number of edges: %d\n" % len(net.edges()))


    def init_edge_weight_figure(self):
        fig, ax = plt.subplots()
        ax.set_title("Interactome Edge Weight Distribution")
        ax.set_xlabel("Edge weight (bin size = .05)")
        ax.set_ylabel("# of edges in bin")

        return fig, ax


    def save_edge_weight_distribution_fig(self, fig, output_dir):
        name = self.interactome.name

        outfile_pdf = \
            Path(output_dir, self.get_output_dir(), name, "edge-weight-distribution.pdf")

        outfile_png = \
            Path(output_dir, self.get_output_dir(), name, "edge-weight-distribution.png")

        outfile_pdf.parent.mkdir(parents=True, exist_ok=True)

        fig.tight_layout()
        fig.savefig(str(outfile_png))
        fig.savefig(str(outfile_pdf))
        

    def plot_edge_weight_distribution(self, output_dir):
        net = self.get_interactome_net()

        weights = [edge[2]["weight"] for edge in net.edges(data=True)]

        fig, ax = self.init_edge_weight_figure()

        ax.hist(weights, bins=np.arange(0, 1, .05), histtype='bar')

        self.save_edge_weight_distribution_fig(fig, output_dir)


    # TODO: purge_results is not implemented
    def run(self, output_dir=Path(), purge_results=False):
        '''
        Gathers statistics on the following:
        1) The number of nodes and edges in the interactome
        2) Interactome edge-weight distributions
        '''
        self.write_summary_statistics(output_dir)
        self.plot_edge_weight_distribution(output_dir)
