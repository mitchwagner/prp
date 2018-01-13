'''
This script's purpose is multiprocess parallelization of the existing pipeline 

Current, arduous steps for running meta-pipeline:

1) Go through each config file in config and select the correct
   subnetwork creation method

2) Run this script, commenting out the final step

3) Comment out the first two steps, comment in the final step,
   and run this script again.

In order to create folds properly, config.yaml and the config files
in ./config need to use the same subnetwork creation method. Otherwise,
the folds won't be created and the algorithms will fail.

Yes, I know this is bad, but it doesn't seem worth it yet to make it
too much better.

TODO: 
Murali had the idea of passing in pathways/algorithms as command-line
parameters, and running only the pathway passed.  That might actually be much
better here. Getting the algorithms to work the same way would be a bit of
a stretch, but it's preferable to creating a million config files.

'''
import subprocess

# local imports
import src.external.utils.baobab_utils as baobab

'''
# First, run the normal, unparallelized pipeline to set things up.
# Don't perform any of the latter steps.
subprocess.call(
    ["python"
    ,"pipeline.py"
    ,"--config=config.yaml"
    ,"--run-reconstructions-off"
    ,"--compute-precision-recall-off"
    ,"--aggregate-precision-recall-folds-off"
    ,"--plot-aggregate-precision-recall-folds-off"
    ,"--aggregate-precision-recall-pathways-off"
    ,"--plot-aggregate-precision-recall-pathways-off"
    ])

# Second, run all of the algorithms on the pathways in parallel on baobab.
# ONLY run the algorithms.

# 1) Write qsub file

# 2) Submit qsub file

for i in range(15):
    config_file = "config/config" + str(i) + ".yaml"
    qsub_file = "/home/mitchw94/Desktop/pathway-reconstruction-pipeline/qsub/qsub" + str(i)

    baobab.writeQsubFile(
        ["cd ~/Desktop/pathway-reconstruction-pipeline"
        ,"source baobab-venv/bin/activate"
        ,"baobab-venv/bin/python pipeline.py " 
        "--config=" + config_file + " "
        "--pathway-specific-interactomes-off "
        "--create-folds-off "
        "--compute-precision-recall-off "
        "--aggregate-precision-recall-folds-off "
        "--plot-aggregate-precision-recall-folds-off "
        "--aggregate-precision-recall-pathways-off "
        "--plot-aggregate-precision-recall-pathways-off "
        ],
        qsub_file,
        name="Run algs for: %s" % config_file,
        ppn=8
        )

    baobab.submitQsubFile(qsub_file)
'''
# Third, create precision/recall measurements and plots in serial
subprocess.call(
    ["python"
    ,"pipeline.py"
    ,"--config=config.yaml"
    ,"--pathway-specific-interactomes-off"
    ,"--create-folds-off"
    ,"--run-reconstructions-off"
    ])
