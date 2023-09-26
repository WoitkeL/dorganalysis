# dorganalysis

Dorganalyses is a tool to find all distributed organizations of a given reaction networks. A distributed organization is a 
set of molecules (or, equivalently, a subnetwork) that can potentially persist forever [1]. The tool and the algorithms used by the toolare intrduced in Ref. [2].

Deorganalysis is inplemented in python (`code` folder) and uses Gurobi for solving mixed linear integer problems. 

Scripts and results from applying the tool to models from the database BioModels can be found in the directory `àppliance_Database`. 

Additional data contains the exel-file as result of the iteration over the Biomodels database, as well as a ipynb file to create meaningful figures.

The following packages have to be installed: 
-LibSBML
-Graphviz 
-regex
-Gurobi: Note that to handle large models, Gurobi has to be installed via the website with a fitting license instead of just by pip.

A step by step showcase of the most applications can be seen in the example.ipynb.

# References

[1] Peter, S., Ibrahim, B., & Dittrich, P. (2021). Linking network structure and dynamics to describe the set of persistent species in reaction diffusion systems. SIAM Journal on Applied Dynamical Systems, 20(4), 2037-2076.

[2] Computing all persistent subspaces of a reaction-diffusion system, submitted, under review, 2023 

