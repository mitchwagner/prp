# Algorithm Runners
from src.runners.reconstruction.NodeAndEdgeWithholdingEvaluatorV4 \
    import NodeAndEdgeWithholdingEvaluatorV4

# Post-hoc Evaluators
from src.runners.post_hoc.FullPathwayEvaluatorV4 \
    import FullPathwayEvaluatorV4

# RWR "q" Estimators 
from src.runners.qestimator.NodeEdgeQEstimatorV3 import NodeEdgeQEstimatorV3

RUNNERS = {
    "node-and-edge-removal-reconstruction-v4": 
        NodeAndEdgeWithholdingEvaluatorV4,

    "node-edge-q-estimator-v3": 
        NodeEdgeQEstimatorV3,

    "full-pathway-v4":
        FullPathwayEvaluatorV4,
}


