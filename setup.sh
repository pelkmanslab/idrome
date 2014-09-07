#! /bin/bash

if command -v conda >/dev/null 2>&1; then
	echo "Anaconda installed! Let's go"
else
	echo "Anaconda needs to be in your path. I recommend using the Miniconda installer. See http://conda.pydata.org/miniconda.html for details."
	exit 1;
fi

echo "~~~~~~~~~~~~~~~~~~~~~~~"
echo "Creating new virtual environment named 'idrome'"
echo "~~~~~~~~~~~~~~~~~~~~~~~"
conda create -n idrome python=2.7 numpy scipy cython pandas matplotlib ipython-notebook pip --yes
echo "Done!"

echo "Switching to virtual environment."
echo echo "~~~~~~~~~~~~~~~~~~~~~~~"
echo "You'll need to do this each time you want to use the notebooks."
echo "This can be done via the 'source activate idrome' command"
echo "~~~~~~~~~~~~~~~~~~~~~~~"

source activate idrome

echo "~~~~~~~~~~~~~~~~~~~~~~~"
echo "Installing non-conda packages via pip"
echo "~~~~~~~~~~~~~~~~~~~~~~~"
pip install seaborn yahmm matplotlib-venn ete2

echo ""
echo "Done! The virtual environment 'idrome' has been successfully created."
echo "~~~~~~~~~~~~~~~~~~~~~~~"
echo "Running the IPython notebook server. This can be done at any time from within the environment, using the 'ipython notebook' command"

ipython notebook
