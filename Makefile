.RECIPEPREFIX = >

netpath-wnt-all-algorithms:
> python master-script.py \
  --ppiversion pathlinker-signaling-children-reg \
  --weightedppi \
  --onlynetpathwnt \
  --pathlinker \
  --shortestpaths \
  --pagerank \
  --responsenet \
  --pcsf \
  --anat \
  --ipa \
  --bowtiebuilder \
  --computeprecrec \
  --precrecviz \
  --forceviz;\
  python master-script.py \
  --ppiversion pathlinker-signaling-children-reg \
  --weightedppi  \
  --onlynetpathwnt \
  --pathlinker \
  --shortestpaths \
  --pagerank \
  --responsenet \
  --pcsf \
  --anat \
  --ipa \
  --bowtiebuilder \
  --computeprecrec \
  --precrecviz \
  --forceviz \
  --ignorekeggpositives

netpath-limited:
> python master-script.py \
  --ppiversion pathlinker-signaling-children-reg \
  --weightedppi \
  --onlynetpathwnt \
  --pathlinker \
  --shortestpaths \
  --computeprecrec \
  --precrecviz\
  --forceviz
