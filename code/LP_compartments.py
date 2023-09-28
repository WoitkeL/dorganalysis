# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:08:37 2023

@author: linus
"""

import gurobipy as gp
from gurobipy import GRB
import regex as re
def get_MCs(RN,SOR, setOfSpecies=set()):
    
    #function to delete subsets
    def delete_smaller_sets(candidate_set):
        candidate_set_remover=set()
        for element in candidate_set:
            for element2 in candidate_set:
                if element!=element2:
                    if element.issubset(element2):
                        candidate_set_remover.add(element)
                        break
        candidate_set.difference_update(candidate_set_remover)
        return(candidate_set)
    
    #function to split speciesset for a not proper reaction
    def splitt_speciesset(element, reaction):
        extender=set()
        for reactant in reaction.listOfReactants:
            #for each reaction we remove each reactant once and add the
            #remaining set to the extender, which is then returned
            extender_element=set(element.copy())
            extender_element.remove(reactant)
            extender.add(frozenset(extender_element))
        return(extender)
    
    
    lOR_gc=RN.reactions[:]
    LOR_dictionary={}
    setOfSpecies_gc=set()
    candidate_set=set()
    inactive_reactions=set()
    SOR_r_class=set()

    #create LOR-dictionary to link name to reaction-class-objekt
    for reaction in lOR_gc:
        LOR_dictionary[reaction.defined_name]=reaction
    
    #set species
    for reaction in SOR:
        setOfSpecies_gc.update(LOR_dictionary[reaction].listOfReactants)  
        setOfSpecies_gc.update(LOR_dictionary[reaction].listOfProducts)
    setOfSpecies_gc.update(setOfSpecies)
    
    #remove uninvolved reactions
    for reaction in lOR_gc:
        if not set(reaction.listOfReactants).issubset(setOfSpecies_gc):
            lOR_gc.remove(reaction)
            continue

    #create list of active and inactive reactions
    for reaction in lOR_gc:
        if reaction.defined_name in SOR:
            SOR_r_class.add(reaction)
        else:
            inactive_reactions.add(reaction)
    
    #start MCs with set of all species as candidate set
    candidate_set.add(frozenset(setOfSpecies_gc.copy()))
    counter=0
    #iterate over every inactive reaction to separate the reactants
    for inactive_reaction in inactive_reactions:
        counter+=1      
        candidate_set_copy=candidate_set.copy()
        #iterate over all candidate sets to seperate reactants of inactive reactions
        for element_set in candidate_set_copy:
            if set(inactive_reaction.listOfReactants).issubset(element_set):

                #add all smaller candidates, which obey the absence of the support of the current inactive reaction   
                candidate_set.update(splitt_speciesset(element_set,inactive_reaction))
                
                #delete all old candidates, which contain the support of the current inactive reaction
                candidate_set.remove(element_set)
 

        if counter%5==0:
            candidate_set= delete_smaller_sets(candidate_set)

    
    #delete subsets in set
    candidate_set= delete_smaller_sets(candidate_set)
    #check for closed
    counter=0

   
    #iterate over every active reaction to check for closedness of every compartment
    nochange=False
    while nochange==False:
        nochange_current=True
        for reaction in SOR_r_class:
            
            counter+=1
            candidate_set_copy=candidate_set.copy()
            for element_set in candidate_set_copy:    
                if set(reaction.listOfReactants).issubset(element_set):
    
                    if not set(reaction.listOfProducts).issubset(element_set):
                        #add all smaller candidates, which obey the absence of the support of the current non-closed reaction 
                        candidate_set.update(splitt_speciesset(element_set, reaction))
                        nochange_current=False
                        #delete all old candidates, which contain the support of the current non-closed reaction 
                        candidate_set.remove(element_set)
        if nochange_current==True:
            nochange=True
      #  if counter%50==0:
       #     candidate_set= delete_smaller_sets(candidate_set)

            
    candidate_set= delete_smaller_sets(candidate_set)

    #setup for LP_solver
    candidate_reactions_dictionary={}
    reaction_compartment_dictionary={}
    species_compartment_dictionary={}
    for reaction in SOR:
        reaction_compartment_dictionary[reaction]=[]
    for species in setOfSpecies_gc:
        species_compartment_dictionary[species]=[]
    for element in candidate_set:
        reactionlist=[]
        for reaction in SOR_r_class:
            if set(reaction.listOfReactants).issubset(element):
                reactionlist.extend(reaction.defined_name)
                reaction_compartment_dictionary[reaction.defined_name].append(element)
        for species in element:
            species_compartment_dictionary[species].append(element)
    
    
 
    
    return(candidate_set)

def get_min_compartments(RN, MCs, SOR, setOfSpecies_input=""): 
    SOR_r_class=set()
    
    for reaction in RN.reactions:
        if reaction.defined_name in SOR:
            SOR_r_class.add(reaction)
    
    setOfSpecies=set()    
    for reaction in SOR_r_class:
        setOfSpecies.update(reaction.listOfReactants)  
        setOfSpecies.update(reaction.listOfProducts)
    setOfSpecies.update(setOfSpecies_input)
    
    #setup for LP_solver
    candidate_reactions_dictionary={}
    reaction_compartment_dictionary={}
    species_compartment_dictionary={}
    for reaction in SOR:
        reaction_compartment_dictionary[reaction]=[]
    for species in setOfSpecies:
        species_compartment_dictionary[species]=[]
    for element in MCs:
        reactionlist=[]
        for reaction in SOR_r_class:
            if set(reaction.listOfReactants).issubset(element):
                reactionlist.extend(reaction.defined_name)
                reaction_compartment_dictionary[reaction.defined_name].append(element)
        for species in element:
            species_compartment_dictionary[species].append(element)
            
            

    
    
    iterator=set()
    compartment_shortform={}
    compartment_shortform_return={}
    i=0
    #setup of dictionaries used constraint construction
    for ding in MCs:
        #maps an index to each of the candidates in both directions
        compartment_shortform[i]= ding
        compartment_shortform_return[str(ding)]= i
        i+=1
        iterator.add(str(ding))
        
    # setup variables including compartments as instance of existence of a speciesset
    #also setup of reaction and species variable, which have to be present in at least one active compartment
        
        
    model = gp.Model()  
    reac = {rname: model.addVar(lb=0, name="r_{}".format(rname)) for rname in SOR}
    spec = {sname: model.addVar(lb=0, name="s_{}".format(sname)) for sname in setOfSpecies} 
    
    compartments = {comp: model.addVar(vtype=gp.GRB.BINARY, name="compartment_{}".format(comp)) for comp in compartment_shortform.keys()}
    
    #Reac=pulp.LpVariable.dicts("reaction_value",SOR,lowBound=0, cat="Integer")
    #species=pulp.LpVariable.dicts("species_value", frozenset(setOfSpecies),lowBound=0,cat="Integer")
    #compartments=pulp.LpVariable.dicts("compartment", compartment_shortform.keys(),lowBound=0, upBound=1,cat="Integer" )
    
    #define lp objective function for minimization of the number of compartments active
    #solve_min_comp=pulp.LpProblem("Solver Min Compartment",pulp.LpMinimize)
    model.setObjective(gp.quicksum(compartments.values()), gp.GRB.MINIMIZE)
    model.setParam( 'OutputFlag', False )
    model.params.TimeLimit = 300
    #objective_function=pulp.lpSum(compartments.values())  
    #solve_min_comp+=objective_function  
    
    #setup of reaction and species in relation to their compartments, in which the occur
    for reaction in SOR:
        model.addConstr(gp.quicksum(compartments[compartment_shortform_return[str(b)]] for b in reaction_compartment_dictionary[reaction])>=1)
    for single_species in setOfSpecies:
        model.addConstr(gp.quicksum(compartments[compartment_shortform_return[str(b)]] for b in species_compartment_dictionary[single_species])>=1)
    result= model.optimize()
    #print(result)
            
    if model.status == gp.GRB.OPTIMAL:
        pass
        #print("The model has an optimal solution with objective value", model.objVal)
    else:
        print("infeasible")
        print(model.status)
        return()     
    if model.status == gp.GRB.TIME_LIMIT:
        return()
    #variable extraction and output
    solutionlist=[]
    species_over_list=[]
    pattern1 = re.compile("compartment")
    for variable in model.getVars():
#        print("{} = {}".format(variable.name,variable.varValue))
        if pattern1.match(variable.VarName) and variable.X==1:
            solutionlist.append(compartment_shortform[int(variable.VarName[12:])])
    
    return(solutionlist)

"######################################################################################################################"
