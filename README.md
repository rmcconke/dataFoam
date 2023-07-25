# `dataFoam`: a set of tools for converting OpenFOAM fields into data for machine learning

The goal of this repository is to:
1. Write certain fields of interest to a set of OpenFOAM cases
2. Read these fields in python, and save them as numpy binaries
3. Interpolate fields from fine meshes (e.g. LES, DNS) to coarse meshes (e.g. RANS)
4. (Optional) Save a csv file containing the final numpy fields as columns.

An example is provided in `example.py`. The example data in this repo comes from the dataset by Xiao et al. https://github.com/xiaoh/para-database-for-PIML, https://doi.org/10.1016/j.compfluid.2020.104431. 

You need OpenFOAM installed to use this repository. There are three OpenFOAM applications which need to be compiled. The source codes are in `foam_applications/`. The script `compile_foam_applications.sh` should be able to compile these applications for you.

A series of unit tests have been written to verify the fields are being read correctly and certain properties of calculated fields are satisfied.

If you are just looking for extracting basic OpenFOAM fields, you can also use this repository. The extra fields written are input features and target quantities for turbulence modelling via machine learning. For more information on these fields, see https://www.nature.com/articles/s41597-021-01034-2, https://www.kaggle.com/datasets/ryleymcconkey/ml-turbulence-dataset. Note: some of the scalings used for input features are different in this repository. To use these fields for machine learning, make sure you understand the fields being written by the corresponding writeFields_* utility. For more information, please contact rmcconke@uwaterloo.ca, or open an issue here.

This is research code, and I have tried where necessary to document and comment it. However, if you are interested in helping me with this repository, I would appreciate the help. 



