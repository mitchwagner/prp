'''
The purpose of this file is to use the pathway edges files to create entries
to add to the interactome file.
'''
import csv
from pathlib import Path


pathway_dir = Path("filtered-pathways-s-t-pruned")

# The list of pathway file names 

pathway_names = \
    ["BDNF"
    ,"EGFR1"
    ,"IL1"
    ,"IL2"
    ,"IL3"
    ,"IL6"
    ,"IL-7"
    ,"KitReceptor"
    ,"Leptin"
    ,"Prolactin"
    ,"RANKL"
    ,"TCR"
    ,"TGF_beta_Receptor"
    ,"TNFalpha"
    ,"Wnt"
    ]

# bi-directional interactions
physicalKeywords = set(['Interaction', 'physical']) 

# directed interactions                                                         
catalysisKeywords = set(['Phosphorylation',
                         'Dephosphorylation',
                         'Proteolytic',
                         'Ubiquitination',
                         'Deubiquitination',
                         'Glycosylation',
                         'Acetylation',
                         'Deacetylation',
                         'Palmitoylation',
                         'Sumoylation',
                         'Desumoylation', 
                         'Methylation',
                         'protein cleavage',
                         'protein_cleavage',
                         'n-palmitoyl-cysteine',
                        ])

NetPathNames = {
    #"NetPath_137" : "Advanced glycation end-products (AGE/RAGE)",
    #"NetPath_1"   : "Alpha6 Beta4 Integrin",
    #"NetPath_2"   : "Androgen receptor (AR)",
    #"NetPath_12"  : "B cell receptor (BCR)",
    "NetPath_76"  : "Brain-derived neurotrophic factor (BDNF)",
    #"NetPath_129" : "Corticotropin-releasing hormone (CRH)",
    "NetPath_4"   : "Epidermal growth factor receptor (EGFR)",
    #"NetPath_25"  : "Follicle-stimulating hormone (FSH)",
    #"NetPath_134" : "Fibroblast growth factor-1 (FGF1)",
    #"NetPath_154" : "Gastrin",
    #"NetPath_118" : "Ghrelin",
    #"NetPath_10"  : "Hedgehog",
    #"NetPath_5"   : "Inhibitor of differentiation (ID)",
    "NetPath_13"  : "Interleukin-1 (IL-1)",
    "NetPath_14"  : "Interleukin-2 (IL-2)",
    "NetPath_15"  : "Interleukin-3 (IL-3)",
    #"NetPath_16"  : "Interleukin-4 (IL-4)",
    #"NetPath_17"  : "Interleukin-5 (IL-5)",
    "NetPath_18"  : "Interleukin-6 (IL-6)",
    "NetPath_19"  : "Interleukin-7 (IL-7)",
    #"NetPath_20"  : "Interleukin-9 (IL-9)",
    #"NetPath_132" : "Interleukin-10 (IL-10)",
    #"NetPath_147" : "Interleukin-11 (IL-11)",
    "NetPath_6"   : "Kit Receptor",
    "NetPath_22"  : "Leptin",
    #"NetPath_3"   : "Notch",
    #"NetPath_114" : "Oncostatin-M (OSM)",
    "NetPath_56"  : "Prolactin",
    "NetPath_21"  : "Receptor activator of nuclear factor kappa-B ligand (RANKL)",
    "NetPath_11"  : "T Cell Receptor (TCR)",
    "NetPath_7"   : "Transforming growth factor beta (TGF-beta) receptor",
    #"NetPath_138" : "Tyrosine kinase with angiopoietins (TIE2/TEK)",
    "NetPath_9"   : "Tumor necrosis factor (TNF) alpha",
    #"NetPath_23"  : "Thyroid-stimulating hormone (TSH)",
    #"NetPath_24"  : "Thymic stromal lymphopoietin (TSLP)",
    #"NetPath_26"  : "TNF-related weak inducer of apoptosis (TWEAK)",
    "NetPath_8"   : "Wnt",
}

def get_is_directed(keyword): 
    if keyword in physicalKeywords:
        return False
    else:
        return True

# First, read the evidence file to get a list of nodes in the interactome

evidence_file = Path("2018_01pathlinker-no-kegg-spike.tsv")

# Read in the evidence file and get a set of nodes


nodes = set()
counter = 0
with evidence_file.open('r') as f:
    for line in f:
        counter += 1
        # Ignore comments
        if line.startswith("#"):
            continue
        else:
            toks = [col.strip() for col in line.split("\t")]
            tail = toks[0]
            head = toks[1]
            nodes.add(tail)
            nodes.add(head)

print("Number of lines:", counter)
print("Number of nodes in the evidence file:", len(nodes))

# Next, add pathway edges to an "extra evidence" file if both nodes are
# part of the interactome proper

# The file we will write extra lines of evidence to
extra_evidence_file = Path("extra-evidence.txt")

bad_matches = 0
with extra_evidence_file.open("w") as f1:
    writer = csv.writer(f1, delimiter="\t")

    for name in pathway_names:
        
        old_edge_file = Path(pathway_dir, name + "-edges.txt")

        with old_edge_file.open('r') as f2:
            for line in f2:
                # Ignore comments
                if line.startswith("#"):
                    continue
                else:
                    toks = line.split("\t")
                    tail = toks[0]
                    head = toks[1]

                    if tail in nodes and head in nodes:

                        interaction_type = toks[5]
                        is_directed = get_is_directed(interaction_type)

                        pathway_id = toks[4]
                        full_pathway_id = "netpath:NetPath_" + pathway_id

                        pathway_name = NetPathNames["NetPath_" + pathway_id]

                        # The evidence file has 7 columns. Write out the values
                        # for these columns, for this line in an edge file,
                        # to the list of extra lines of evidence.

                        col = [tail
                              ,head
                              ,str(is_directed)
                              ,interaction_type
                              ,pathway_name
                              ,full_pathway_id
                              ,"NetPath"
                              ]
                        
                        writer.writerow(col)

                    else:
                        bad_matches += 1

print(bad_matches, "lines where nodes where not in the interactome")

# We can then use the resulting evidence file to get a list of directions
# for every edge in the interactome

# Finally, add the extra evidence lines to the evidence file and rebuild and
# reweight the interactome
