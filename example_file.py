# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:40:35 2023

@author: linus
"""
import os
from analyse_class import Analyse
from reaction_network import Reaction,get_example

"""def remove_duplicates(sorted_list):
    previous_element = None
    i = 0
    while i < len(sorted_list):
        current_element = sorted_list[i]
        if current_element == previous_element:
            del sorted_list[i]
        else:
            previous_element = current_element
            i += 1
    return sorted_list"""
def remove_duplicates(lst):
    unique_lst = []
    for item in lst:
        if item not in unique_lst:
            unique_lst.append(item)
    return unique_lst

bla1=get_example()
print("teetet")
print(bla1)
print(str(bla1))
print("teetet")
x=Analyse(bla1)
x2=Reaction("r14",["b21"],["b1" , "a1"],[1],[1,1])
bla1.add_reaction(x2)
#example1=Analyse()
x.all_SARs()
z=x.SAR_solution
x.all_SARs()
x.all_DOs()
y=x.DOs

#y=computeDO(x)
#z=computeSOR(x)
print(y)
print(z)
set_xbla=set()
print("fasle")
for ele in z:
    print(ele)
    print(get_DOs_of_SOR(x,ele))
    #set_xbla.update(get_DOs_of_SOR(x,ele))
#xbla=get_DOs_of_SOR(x,z[-1])
print(set_xbla)
print("wut?")
print(y)
print("rly?")
print(len(set_xbla))
xa=set()
for ele in y:
    xa.add(frozenset(ele))

print(len(xa))
zaza=[tuple(i) for i in y]
print("mecl")
for i in xa:
    if i not in set_xbla:
        print(i)
print(os.getcwd())
x.draw_hasse()
# or initializing with path to sbml file
#path=.../...xml
#example1=Analyse(path)

#you can change the handling of reverse reactions in the following ways:
#example1=Analyse(path,consider_reverse=False)

#or make it one directional: 
#but note that you need to use naive algorithm to draw the hasse diagram
#example1=Analyse(path,alternative_reverse=True)

#you can also create inflow reactions for constant species or speices with a concentration higher than 0:
#example1=Analyse(path,consider_constant=True,consider_init_ammount=True)


#options to change reaction network
#example1.change_network(remove_species_from_reactions="EmptySet")
#example1.change_network(change_inflow=True)
"""
#print reaction network
example1.print_RN()



#get list of ERCs
example1.print_ERC()

#get list of species
example1.get_species()
print(example1.species)



#calculating DOS and SARs respectively

example1.SARs()
example1.DOs(second_optimize=False,use_gurobi=True)

print()
#solve for SARs with overproduction
#example1.SARs(second_optimize=True,use_gurobi=False)

#example1.DOs()
numb=calculate_os(example1.DO_solution,example1.reaction_network)
print("guszav")
print(numb)
example1.DO_solution=remove_duplicates(example1.DO_solution)
#print(example1.SAR_solution)
print("gnampf")
print(type(example1.DO_solution[0][0]))
print(len(example1.DO_solution))
print(example1.DO_solution)
aba=example1.DO_solution.copy()
example1.DOs(second_optimize=False,use_gurobi=False)
print(len(example1.DO_solution))
print(example1.DO_solution)

print("chick")
for ele in aba:
    if ele not in example1.DO_solution:
        print(ele)
#print(example1.evaluate_solution())
# calculating all minimal compartments
#mrc numbers can give number of MRCs before using get_min_compartments()
# mrc_numbers=[]

# for i in example1.SAR_solution:
#     bol,boos=check_ORG_LOR_expanded(i[0],example1.reaction_network)
#     if not bol:
#         solution_gc=get_compartments(example1.reaction_network,i[0])
#         mrc_numbers.append(len(solution_gc[1]))
#         compartments=get_min_compartments(*solution_gc)
#         print(compartments)
        
#draw functions here you can see all optional parameters:
#example1.draw_hasse(solution="SARs",use_naive=False,second_value=False,show_species=False,show_compartments=False,show_new=False,shortform=False)
#this command line is the same as
#example1.draw_hasse()
#calculate_number_DO(example1.species,example1.SAR_solution,example1.reaction_network)
#use this to draw DOs:
example1.draw_hasse(solution="DOs",use_naive=True,show_compartments=False, second_value=False)

#used to get partial data of exel file
#analyze_dict= analyze_DO(example1.SAR_solution, example1.reaction_network)
#print(analyze_dict)
"""

