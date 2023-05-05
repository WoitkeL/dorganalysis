# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 15:18:10 2023

@author: linus
"""

from dorganalysis import *
from reactionnetwork import *


example=get_example()

analysis=Analysis(example)
print(example)

allDOs=analysis.all_DOs()

print(len(allDOs))
print(allDOs)
allSORs=analysis.all_SORs()
print(allSORs)
#x.draw_hasse("DOs")
largestSOR=allSORs[-1]
DOs_of_largest_SOR=analysis.get_DOs_of_SOR(largestSOR)
print(DOs_of_largest_SOR)
analysis.draw_DOs()
analysis.draw_SORs(show_compartments=True)
#analysis.draw_DO()
#analysis.draw_DO(DOs_of_largest_SOR)
print("omegalul")
ads=analysis.get_compartmentalization_of_SOR_DO_pair(allSORs[5])
print("butwhy")
print(ads)
ads=analysis.get_compartmentalization_of_SOR_DO_pair(largestSOR)
print(ads)
print("omegalul")
print(analysis.print_ERC())
for ele in analysis.ERCs:
    print(analysis.ERCs[ele])
print(analysis.SORs)
print("equal")
print(analysis.DOs)
print(DOs_of_largest_SOR)
SOR=allSORs[-1]
DOs_of_largest_SOR=analysis.get_DOs_of_SOR(SOR)
#analysis.draw_DO_subset(DOs_of_largest_SOR)
analysis.draw_DOs(DOs_of_largest_SOR)
analysis.draw_SORs()
analysis.draw_DOs()