# Introduction

**Dorganalyses** is a tool to find all distributed organizations of a given reaction networks. A distributed organization is a 
set of molecules (or, equivalently, a subnetwork) that can potentially persist forever [1]. The tool and the algorithms used by the toolare intrduced in Ref. [2].

Deorganalysis is inplemented in python (`code` folder) and uses Gurobi for solving mixed linear integer problems. 

Scripts and results from applying the tool to models from the database BioModels can be found in the directory `àppliance_Database`. 

# Installation

## Requirements
Deorganalysis requires **python** (e.g., by using the annaconda framework) including **Juypter Notebook** and the following packages: 
- LibSBML: `pip install python-libsbml`
- Graphviz: `pip install graphviz`
- regex: `pip install regex`
- Gurobi: `pip install gurobipy`

Note that to handle large models, Gurobi has to be installed via the website with a fitting license instead of just by pip.

## Download Dorganalysis


`git clone https://github.com/WoitkeL/dorganalysis`




# Example

A step by step showcase of the most applications can be seen in the example.ipynb.

# References

[1] Peter, S., Ibrahim, B., & Dittrich, P. (2021). Linking network structure and dynamics to describe the set of persistent species in reaction diffusion systems. SIAM Journal on Applied Dynamical Systems, 20(4), 2037-2076.

[2] Computing all persistent subspaces of a reaction-diffusion system, submitted, under review, 2023 

