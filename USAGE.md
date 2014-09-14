# Usage Information

Once you have set everything up, open up the terminal and run `source activate idrome` followed by `ipython notebook`. You should see a list of all the IPython notebooks. The [Index.ipynb](Index.ipynb) details which notebook covers which topic. Now, if you want to re-run any of the analyses then you will need to either

- run [Data_Prep.ipynb](Data_Prep.ipynb), which will download the latest Uniprot release, run the disordered analysis, clean up the data, and update any file links. This process can take around 30 minutes to 1 hour.
- if you do not want to re-run the analysis then it's possible to copy the datasets from my folder on share4 and modify the [config.json](config.json) file to point to the location of the datasets. If you choose this approach then I recommend you make a `human/` and `viral/` folders and copy the datasets there. 

Once that is complete, you can re-run any of the analysis by clicking on the notebook and then running it by clicking on the `Run All` option under the `Cell` menu option in the menubar. The notebook will then systematically execute the separate cells and any output can be viewed within the notebook itself. Cells that are currently executing have an * inside the square brackets on the left. All the datasets can be used independently of the code itself and [RAWDATASETS.md](RAWDATASETS.md) and [GENDATASETS.md](GENDATASETS.md) contain the necessary data on that front.

I will try and make sure to keep the [issues page](https://github.com/pelkmanslab/idrome/issues) up-to-date with any current issues. Running `git pull` in the root directory of the idrome project will retrieve any fixes that I have made since you have last updated. If you run into any problems/bugs with the code itself or with the documentation, please create a new issue and I will try to get to it as quickly as possible.
