# Algorithm Runners
from src.runners.reconstruction.NodeAndEdgeWithholdingRunner \
    import NodeAndEdgeWithholdingRunner

# Post-hoc Evaluators
from src.runners.post_hoc.FullPathwayEvaluatorV4 \
    import FullPathwayRunner

# RWR "q" Estimators 
from src.runners.stats.qestimator.NodeEdgeQEstimatorV3 import NodeEdgeQEstimatorV3

RUNNERS = {
    "node-and-edge-removal-reconstruction": 
        NodeAndEdgeWithholdingRunner,

    "node-edge-q-estimator-v3": 
        NodeEdgeQEstimatorV3,

    "full-pathway":
        FullPathwayRunner,
}


