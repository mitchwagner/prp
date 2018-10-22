# Pathway Reconstruction Pipeline 

This repository tracks code used for evaluating algorithms for signaling 
pathway reconstruction.

## Setup

Setup has been mostly automated in the script `setup`. Currently, the pipeline
must be setup somewhere on the Murali group computing cluster, for access to
the SVN repository thereon.

Things that need to be added to the setup script:

    - One must build the quicklinker jar in the subrepo. For this, one needs
      to install Gradle and run "gradle build" in the subrepo TLD. 

    - The subrepositories seem to detach themselves (e.g., are not tracking
      any particular remote branch). It would be good to fix this.

## Running

The pipeline should be executed as follows:

```bash
source venv/bin/activate
python pipeline.py
```

A limited set of command line parameters can be passed as well. Run 
`pipeline.py --help` for more information.


## Configuration 

The pipeline is configured via the config file. By default, the pipeline reads
from config.yaml. Instructions for setting various parameters can be found in
that file.

## Inputs

- *Interactomes*: stored as edge lists, as constructed from CSBDB.

- *Pathways*: stored as node and edge lists:
    - \<pathway-name>-nodes.txt
    - \<pathway-name>-edges.txt

- *Pathway Collections*: pathway collections consist of two things: 
    - A folder, where the pathways in the collection are stored
    - A list of pathway names, whose node and edge lists can be found in the
      pathway collection's folder
