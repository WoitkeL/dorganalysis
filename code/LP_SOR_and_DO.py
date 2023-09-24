# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:07:39 2023

@author: linus
"""

import time
import regex as re
from reactionnetwork import  ERC_meets_transitivity,generate_closure_for_reactions,generate_closure_for_species
from gurobipy import GRB
import gurobipy as gp
import sys




def basicLP(analyze,DO=False):
    
    ERC_dict=analyze.get_ERCs()
    lOR_solve=analyze.reaction_network.reactions
    
    setOfSpecies=analyze.reaction_network.species

    #generates the reaction names of the list of reactions
    listOfReactionsname = [i.defined_name for i in lOR_solve]
    #the exact stoichiometric parameters for reactions and the resulting sspecies quantity



    model = gp.Model()
    model.setParam(GRB.Param.OutputFlag, 0)
    
    reac = {rname: model.addVar(lb=0, name="r_{}".format(rname)) for rname in listOfReactionsname}
    
    # the boolean value for reactions
    Rbool = {rname: model.addVar(vtype=gp.GRB.BINARY, name="Rbool_{}".format(rname)) for rname in listOfReactionsname}
    """for DO"""
    if DO==True:
        species_exist = {sname: model.addVar( vtype=gp.GRB.BINARY, name="species_exist_{}".format(sname)) for sname in setOfSpecies}

    # defines optimization goal for LP problem as maximization problem
    """for DO"""
    if DO==True:
        model.setObjective(gp.quicksum(species_exist.values()), gp.GRB.MAXIMIZE)
    else:
        model.setObjective(gp.quicksum(Rbool.values()), gp.GRB.MAXIMIZE)
 
    speciesdict={}
    
    # Initialize empty lists for reactants and products of each species
    for species in setOfSpecies:
        speciesdict[species+'_r']=[]
        speciesdict[species+'_p']=[]  
    
    # Map each species to its reactions as either reactant or product
    for reaction in lOR_solve:
        for i in range(len(reaction.listOfReactants)):
            speciesdict[reaction.listOfReactants[i]+'_r'].append((reaction.defined_name,reaction.reac_stoich[i]))
        for i in range(len(reaction.listOfProducts)):
            speciesdict[reaction.listOfProducts[i]+'_p'].append((reaction.defined_name,reaction.prod_stoich[i]))    
    
    # Implement the stoichiometric matrix as equations for each species
    for ele in setOfSpecies:
        educt=[]
        product=[]
        for combi in speciesdict[ele+'_r']:
            educt.append(int(combi[1])*reac[combi[0]])
        for combi in speciesdict[ele+'_p']:
            product.append(int(combi[1])*reac[combi[0]])
            
        # Group the list of reactions to one expression using quicksum
        # Implement the equations for the changes of species
        model.addConstr(gp.quicksum(product) - gp.quicksum(educt) >= 0)
        
    #relate bool variable with reaction flux
    for key, element in reac.items():
        model.addConstr(Rbool[key]<=element)
        model.addConstr(element<=Rbool[key]*1000)
        
        
    if type(ERC_dict)==str:
        print("timeout ERC")
        return()
    
    for reaction1 in lOR_solve:
        for x in range (len(ERC_dict[reaction1.defined_name].reactions)-1):
            #solve_dist_org+= Rbool[reaction1.defined_name]<= Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]
            model.addConstr(Rbool[reaction1.defined_name]<= Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name])
            
        #a reaction which is not closed in regard to the species set has to be inactive
        if not reaction1.closed:
            model.addConstr( Rbool[reaction1.defined_name]==0)
        if reaction1.always:
            model.addConstr(Rbool[reaction1.defined_name]==1)
                
            
    """for DO"""
    if DO==True:
        spec_closure=generate_closure_for_species(analyze.reaction_network)
        for species in setOfSpecies:
            for x in range (len(spec_closure[species].reactions)):
                model.addConstr(species_exist[species]<= Rbool[spec_closure[species].reactions[x].defined_name])
                
# also sets dependance of species_exist to the occurence in reactions. Since if a reaction is active, the
#involved species are part of the DO
    """for DO"""
    if DO==True:
        for reaction2 in lOR_solve:
            for species in reaction2.listOfReactants:
                model.addConstr(species_exist[species]>= Rbool[reaction2.defined_name])
            for species in reaction2.listOfProducts:
                model.addConstr(species_exist[species]>= Rbool[reaction2.defined_name])





     
    sol_pool_var=50000 
    model.params.TimeLimit = 300
    model.params.PoolSearchMode = 2
    model.params.Threads = 8
    model.setParam(GRB.Param.PoolSearchMode, 2)
    model.params.PoolSolutions = int(sol_pool_var)
    model.optimize()
    # Get the time taken to solve the LP problem
    solve_time = model.Runtime
    
    analyze.metadata["Timer_LP"]=solve_time
    analyze.metadata["number_constr"]=len(model.getConstrs())
    model.write("model.lp")
    
    
    print("Time taken to solve LP:", solve_time, "seconds")
    if model.status == gp.GRB.TIME_LIMIT:
        print("Time limit reached")

    if model.status == gp.GRB.OPTIMAL:
        print("The model has an optimal solution with objective value", model.objVal)
    else:
        print("infeasible")
        print(model.status)
        return()
    sys.stdout.flush()
    pool_solutions = model.SolCount

    """for DO"""
    if DO==True:
        sort_para="species_exist"
        sort_para2="Rbool"
        pattern2 = re.compile(sort_para2)
    else:
        sort_para="Rbool"
    pattern1 = re.compile(sort_para)
    second_value=[]
    output=[]
    
    #iterate over each solution
    for i in range(pool_solutions):
        current_solution=[]
        second_value=[]
        #set the index of the solution to the currently available solution
        model.setParam(gp.GRB.Param.SolutionNumber, i)
        model.setParam(GRB.Param.SolutionNumber, i)
        for v in model.getVars():
            if pattern1.match(v.VarName) and int(v.Xn)==1:
                current_solution.append(v.VarName[(len(sort_para)+1):])
            if DO==True:
                if pattern2.match(v.VarName) and int(v.Xn)==1:
                    second_value.append(v.VarName[(len(sort_para2)+1):])
        if len(current_solution)>0:
            current_solution.sort()
        if len(second_value)>0:
            second_value.sort()
        if DO==True:
            output.append((tuple(current_solution), tuple(second_value)))
        else:
            output.append(tuple(current_solution))
   
    
    return(output)

def OP_LP(analyze, SOR):
    
    
    lOR_solve=[ele for ele in analyze.reaction_network.reactions if ele.defined_name in SOR]
    #print(lOR_solve)
    #lOR_solve=SOR
    setOfSpecies=analyze.reaction_network.species

    #generates the reaction names of the list of reactions
    #listOfReactionsname = SOR
    #the exact stoichiometric parameters for reactions and the resulting sspecies quantity

    model = gp.Model()
    model.setParam(GRB.Param.OutputFlag, 0)
    
    reac = {rname: model.addVar(lb=1, name="r_{}".format(rname)) for rname in SOR}
    
    # the boolean value for reactions
    #species_exist = {sname: model.addVar( vtype=gp.GRB.BINARY, name="species_exist_{}".format(sname)) for sname in setOfSpecies}
    species_ammount = {sname: model.addVar(lb=0,  name="s_ammount_{}".format(sname)) for sname in setOfSpecies}
    sbool={sbool: model.addVar(vtype=gp.GRB.BINARY, name="sbool_{}".format(sbool)) for sbool in setOfSpecies}
        
    
    # defines optimization goal for LP problem as maximization problem
    
    model.setObjective(gp.quicksum(sbool.values()), gp.GRB.MAXIMIZE)
    
    
    speciesdict={}
    
    # Initialize empty lists for reactants and products of each species
    for species in setOfSpecies:
        speciesdict[species+'_r']=[]
        speciesdict[species+'_p']=[]  
    
    # Map each species to its reactions as either reactant or product
    for reaction in lOR_solve:
        for i in range(len(reaction.listOfReactants)):
            speciesdict[reaction.listOfReactants[i]+'_r'].append((reaction.defined_name,reaction.reac_stoich[i]))
        for i in range(len(reaction.listOfProducts)):
            speciesdict[reaction.listOfProducts[i]+'_p'].append((reaction.defined_name,reaction.prod_stoich[i]))    
    
    # Implement the stoichiometric matrix as equations for each species
    for ele in setOfSpecies:
        educt=[]
        product=[]
        for combi in speciesdict[ele+'_r']:
            educt.append(int(combi[1])*reac[combi[0]])
        for combi in speciesdict[ele+'_p']:
            product.append(int(combi[1])*reac[combi[0]])
            
        # Group the list of reactions to one expression using quicksum
        # Implement the equations for the changes of species
#         print(ele)
#         print(product)
#         print(educt)
#         print(gp.quicksum(product))
#         print(gp.quicksum(educt))
        model.addConstr(species_ammount[ele]==gp.quicksum(product) - gp.quicksum(educt))

    #relate bool variable with reaction flux
    for key, element in species_ammount.items():
        model.addConstr(sbool[key]<=element)
        model.addConstr(element<=sbool[key]*1000)
           
        
# also sets dependance of species_exist to the occurence in reactions. Since if a reaction is active, the
#involved species are part of the DO
    """for DO"""
    for reaction1 in lOR_solve:
        model.addConstr(reac[reaction1.defined_name]>=1)
            
    
    model.params.TimeLimit = 30
    model.params.Threads = 8
    model.optimize()
    # Get the time taken to solve the LP problem
    solve_time = model.Runtime
   
    if model.status == gp.GRB.TIME_LIMIT:
        print("Time limit reached")

    if model.status == gp.GRB.OPTIMAL:
        print("The model has an optimal solution with objective value", model.objVal)
    else:
        print("infeasible")
        print(model.status)
        return()
    sys.stdout.flush()

    sort_para="sbool"
    pattern1 = re.compile(sort_para)
    output=[]
    
    #iterate over each solution
    
    current_solution=[]
    for v in model.getVars():
        if pattern1.match(v.VarName) and int(v.Xn)==1:
            current_solution.append(v.VarName[(len(sort_para)+1):])
   
    
    return(current_solution)


