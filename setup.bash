###############################################################################
# Pull submodules

git submodule update --init --recursive

# Build gradle after downloading QuickLinker
gradle build -p src/external/quicklinker/

# TODO: Pypy!

###############################################################################
# Create Python virtual environment for the pipeline

virtualenv venv -p python3

source venv/bin/activate

pip install -r requirements.txt

deactivate

# Create Python virtual environment for algorithms in the pipeline

virtualenv venv-regpathlinker -p python2

source venv-regpathlinker/bin/activate

pip install -r requirements-regpathlinker.txt

deactivate

###############################################################################
# TODO: This likely won't work unless this is actually run on baobab
# Create Python virtual environment for running pipeline on baobab 
# (local compute cluster)

#python3.6 -m venv baobab-venv
#
#source baobab-venv/bin/activate
#
#pip install -r requirements.txt
#
#deactivate
#
## Create Python virtual environment for running algorithms in the pipeline
## on baobab (local compute cluster)
#
## NOTE: Right now, to run on baobab, you will have to manually change the
## venv sourced in the code to run this algorithm...far less than ideal...
#
#python3.6 -m venv baobab-venv-regpathlinker -p python2
#
#source baobab-venv-regpathlinker/bin/activate
#
#pip install -r requirements-regpathlinker.txt
#
#deactivate

###############################################################################
# Fetch the appropriate input from SVN

svn checkout file:///home/murali/svnrepo/data/interactions inputs/interactions

svn checkout file:///home/murali/svnrepo/data/interactomes inputs/interactomes

svn checkout file:///home/murali/svnrepo/data/pathways inputs/pathways
