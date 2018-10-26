###############################################################################
# Pull submodules

git submodule update --init --recursive

# Build gradle after downloading QuickLinker
gradle build -p src/external/quicklinker/

###############################################################################
# Create Python virtual environment for the pipeline

virtualenv venv -p python3

source venv/bin/activate

pip install -r requirements/requirements.txt

deactivate

# Create Python virtual environment for algorithms in the pipeline

virtualenv venv-regpathlinker -p python2

source venv-regpathlinker/bin/activate

pip install -r requirements/requirements-regpathlinker.txt

deactivate

###############################################################################
# Fetch the appropriate input from SVN

svn checkout file:///home/murali/svnrepo/data/interactions inputs/interactions

svn checkout file:///home/murali/svnrepo/data/interactomes inputs/interactomes

svn checkout file:///home/murali/svnrepo/data/pathways inputs/pathways

svn checkout file:///home/murali/svnrepo/data/transcription-factors inputs/transcription-factors

svn checkout file:///home/murali/svnrepo/data/receptors inputs/receptors

