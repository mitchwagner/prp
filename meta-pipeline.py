'''
This script's purpose is multiprocess parallelization of the existing pipeline 
'''
import subprocess

# local imports
import src.external.utils.baobab_utils as baobab

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
        ,"source venv/bin/activate"
        ,"python pipeline.py " 
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
'''
