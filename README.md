# Introduction

**Dorganalyses** is a tool to find all distributed organizations of a given reaction networks [1]. A distributed organization is a 
set of molecules (or, equivalently, a subnetwork) that can potentially persist forever [2]. 

Deorganalysis is inplemented in python (`code` folder) and uses Gurobi for solving mixed linear integer problems.  

# Installation

## Requirements
Deorganalysis requires **python** (e.g., by using the annaconda framework) including **Jupyter Notebook** and the following packages: 
- LibSBML: `pip install python-libsbml`
- Graphviz: `pip install graphviz`
- regex: `pip install regex`
- Gurobi: `pip install gurobipy`

Note that to handle large models, Gurobi has to be installed via the website with a fitting license instead of just by pip.

## Download Dorganalysis


`git clone https://github.com/WoitkeL/dorganalysis`




# Example

Go to the `code` directory and open ` example.ipynb` using `Jupyter Notebook`.  It contains a step by step walkthrough of a typical analysis.

# Analysis of BioModels

1. Download models from [BioModels](https://www.ebi.ac.uk/biomodels/) and place them in a dedicated directory.
1. Change the path in `make_biomodels_csv.py` to that directory and run the script.

Note that the computation of a particular model can take long (several hours or days). There are several timeouts to limit the compute time. See `model.params.TimeLimit`

# Citation
If you find the tool useful, please cite the companion paper:

```
Peter, S., Woitke, L., Dittrich, P., Ibrahim, B. (2023). Computing all persistent subspaces of a reaction-diffusion system, JOURNAL, submitted, under review
```


# References

[1] Peter, S., Woitke, L., Dittrich, P., Ibrahim, B. (2023). Computing all persistent subspaces of a reaction-diffusion system, JOURNAL, submitted, under review

[2] Peter, S., Ibrahim, B., & Dittrich, P. (2021). Linking network structure and dynamics to describe the set of persistent species in reaction diffusion systems. SIAM Journal on Applied Dynamical Systems, 20(4), 2037-2076.

