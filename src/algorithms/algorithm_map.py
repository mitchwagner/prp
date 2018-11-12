'''
Exports a dictionary that maps strings to Python classes.
'''

from typing import Dict, Type

# Bookkeeping/sanity checks 
import src.algorithms.RankingAlgorithm as RankingAlgorithm
import src.algorithms.QuickRegLinkerSanityCheck as SanityCheck 

# PathLinker
import src.algorithms.PathLinker as PathLinker
import src.algorithms.PathLinkerRWER as PathLinkerRWER

# Induced Subgraph
import src.algorithms.InducedSubgraph as InducedSubgraph
import src.algorithms.GenInduced as GenInduced
import src.algorithms.GenInducedRWR as GenInducedRWR
import src.algorithms.GenInducedRWER as GenInducedRWER
 
import src.algorithms.InducedRWER as InducedRWER
import src.algorithms.InducedRWR as InducedRWR

# Shortcuts and generalized shortcuts
import src.algorithms.Shortcuts as Shortcuts
import src.algorithms.ShortcutsRWER as ShortcutsRWER
import src.algorithms.ShortcutsRWR as ShortcutsRWR

import src.algorithms.GeneralizedShortcuts as GeneralizedShortcuts
import src.algorithms.GeneralizedShortcutsRWER as GeneralizedShortcutsRWER

# ZeroQuickLinker
import src.algorithms.ZeroQuickLinkerLabelNegatives as \
    ZeroQuickLinkerLabelNegatives

# Random Walks
import src.algorithms.RWR as RWR
import src.algorithms.RWER as RWER

# Final version of RegLinker 
import src.algorithms.RegLinker as RegLinker
import src.algorithms.RegLinkerPaths as RegLinkerPaths
import src.algorithms.RegLinkerRWER as RegLinkerRWER
import src.algorithms.RegLinkerRWERPaths as RegLinkerRWERPaths
import src.algorithms.RegLinkerRWERNoLoops as RegLinkerRWERNoLoops
import src.algorithms.RegLinkerRWR as RegLinkerRWR

import src.algorithms.RegLinkerBetter as RegLinkerBetter

RANKING_ALGORITHMS: Dict[str, Type[RankingAlgorithm.RankingAlgorithm]] = {
    "quickreglinker-sanity": SanityCheck.QuickRegLinkerSanityCheck,

    "induced-subgraph": InducedSubgraph.InducedSubgraph,
    "InducedRWR": InducedRWR.InducedRWR,
    "InducedRWER": InducedRWER.InducedRWER,

    "GenInduced": GenInduced.GenInduced,
    "GenInducedRWR": GenInducedRWR.GenInducedRWR,
    "GenInducedRWER": GenInducedRWER.GenInducedRWER,

    "Shortcuts": Shortcuts.Shortcuts,
    "ShortcutsRWR": ShortcutsRWR.ShortcutsRWR,
    "ShortcutsRWER": ShortcutsRWER.ShortcutsRWER,

    "GeneralizedShortcuts": GeneralizedShortcuts.GeneralizedShortcuts,
    "GeneralizedShortcutsRWER": 
        GeneralizedShortcutsRWER.GeneralizedShortcutsRWER,

    "PathLinker": PathLinker.PathLinker,
    "PathLinkerRWER": PathLinkerRWER.PathLinkerRWER,

    "ZeroQuickLinkerLabelNegatives": 
        ZeroQuickLinkerLabelNegatives.ZeroQuickLinkerLabelNegatives,

    "RWR": RWR.RWR,
    "RWER": RWER.RWER,
    
    "RegLinker":  RegLinker.RegLinker,
    "RegLinkerPaths":  RegLinkerPaths.RegLinkerPaths,
    "RegLinkerRWR":  RegLinkerRWR.RegLinkerRWR,
    "RegLinkerRWER":  RegLinkerRWER.RegLinkerRWER,
    "RegLinkerRWERPaths":  RegLinkerRWERPaths.RegLinkerRWERPaths,
    "RegLinkerRWERNoLoops":  RegLinkerRWERNoLoops.RegLinkerRWERNoLoops,

    "RegLinkerBetter":  RegLinkerBetter.RegLinkerBetter,
}
