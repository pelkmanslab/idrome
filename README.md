# idrome

[![BSD License](http://img.shields.io/badge/license-BSD 3--clause-blue.svg?style=flat-square)](http://opensource.org/licenses/BSD-3-Clause)
[![Github Issues](http://img.shields.io/github/issues/badges/shields.svg?style=flat-square))](https://github.com/pelkmanslab/idrome/issues)

This project aims to systematically analyze the functional attributes of proteins containing intrinsically disordered regions (IDRs) and low complexity regions (LCRs) in the human proteome and the human virus proteomes (potentially others as well).

It makes extensive use of [IPython notebooks](http://www.ipython.org), which are a convenient way of showing both the code and thinking process/formulas/logic that went into creating said code. The notebooks all have `.ipynb` extensions and are simply JSON files. The best way to interact with the notebooks is to setup the computing environment associated with this project and run the IPython notebook server locally. 

See [USAGE.md](USAGE.md) for usage information.

*Note: This code has been tested on Linux and OS X, it has not been tested on Windows. It is very likely that there will be several PATH related issues. Additionally, there is only a OS X compatible binary for CAST at the moment, which is required in the data acquisition step.*

## Project Structure

- [Index.ipynb](Index.ipynb) contains all the information regarding the organization and structure of the code files.
- [RAWDATASETS.md](RAWDATASETS.md) describes the layout of the raw data files
- [GENDATASETS.md](GENDATASETS.md) describes the layout of generated data files
- [config.json](config.json) keeps track of important files and directories used in the project

## Setup

This project is built on the scientific Python environment and makes heavy use of several libraries (i.e. toolboxes in MATLAB) that will all need to be installed. The [Anaconda package manager](http://docs.continuum.io/anaconda/) can greatly ease the setup and compilation process and is the recommended way of reproducing the project's computing environment. First, clone the project by using either using Github's GUI client (for [Mac](http://mac.github.com) or [Windows](http://windows.github.com)) or running the following command (assuming `git` is installed):
```
git clone https://github.com/pelkmanslab/idrome.git
```
Next, change directories into the newly created "idrome" folder and run 
```
git submodule init
```
and then
```
git submodule update
```
to clone the necessary [submodules](http://www.git-scm.com/book/en/Git-Tools-Submodules). 

### Install Anaconda 
The default Anaconda install includes a very large number of scientific Python packages and is very large. Miniconda, a barebones version that includes only the Python interpreter and the "conda" package, is recommended for quicker installation. Miniconda is available at <http://conda.pydata.org/miniconda.html>. The Python 3.x version is recommended, but the Python 2.7 version will also work because a virtual environment with the correct version of interpreter will be created in either case. Once Miniconda is installed, continue to the next section.

### Setup project environment using setup.sh

Change directories into the folder where the project was cloned and run the "setup.sh" script (e.g. `bash setup.sh`). This creates a virtual environment, called "idrome", that is pre-populated with all the necessary packages. The script will then switch to the new environment and start up the notebook server. In future sessions, you will need to switch to this new environment manually using the `source activate idrome` command. Once that is done, the IPython notebook server can be run by executing `ipython notebook` in the project directory.

### Setup project environment manually

In case the previous script did not work, the process can also be carried out manually. First, we will need to create a virtual environment (called "idrome" in this example) and populate it with the most important packages. Run the following command:

```
conda create -n idrome python=2.7 numpy scipy cython pandas matplotlib ipython-notebook pip pyqt
```
Anaconda does not have binaries for all of the required packages. We will need to install the remainder via the `pip` package manager. But before we can do that, we need to actually switch into the new virtual environment. As mentioned before, this needs to be done every time you want to start the ipython notebook environment. Next, run `conda info -e` to show all environments (which are simply folders with packages). "idrome" should be in the list. Next, run:
```
source activate idrome
```
We are now in the virtual environment. Next, run
```
pip install seaborn yahmm matplotlib-venn ete2
```
And finally we can run `ipython notebook` in the project directory to actually view all the notebooks.
