# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 20:19:25 2023

@author: linus
"""

import os
import time
import csv
from Hasse import check_ORG_LOR,create_Hasse
from Analysis import Analysis
from reactionnetwork import create_closures, SBML_to_RN
from LP_compartments import get_min_compartments,get_MCs

def calculate_os(analyze):
    orgs=0
    for i in range(len(analyze.DOs)):
        if is_DO_O(analyze,i):
            orgs+=1
    return(orgs)

def is_DO_O(analyze,i):
    
    for rea in analyze.reaction_network.reactions:
        if rea.defined_name not in analyze.SOR_map[analyze.DOs[i]][-1]:
            if set(rea.listOfReactants).issubset(analyze.DOs[i]):
                return(False)
    return(True)
    
    
def iterate_over_database(path="C:/Users/linus/python/nice_BM2/", start_index=0,end_index=3000, exelname=False):
    
    header = ['ID',"name",'number_species', 'number_reactions', 'number_SORs',\
              'number_SOR_O', "number_SOR_DO_only","number_DO","number_Os",\
                  "max_MRC","comp2","comp3","comp4","comp5","comp6","timer_MRC","reaction_orders",\
                  "Timer_ERC","Timer_LP","Timer_LP_DO1","number_constr","time_complete"]
    more_timers=False
    if more_timers:
        header = ['ID',"name",'number_species', 'number_reactions', 'number_SORs',\
              'number_SOR_O', "number_SOR_DO_only","number_DO","number_Os",\
                  "max_MRC","comp2","comp3","comp4","comp5","comp6","timer_MRC","reaction_orders",\
                  "Timer_ERC","Timer_LP","Timer_LP_DO1","number_constr","timer_all1","timer_all2","timer_all3","timer_all4"]
    if exelname:
        exel_insert=str(exelname)
    else:
        exel_insert="output_table"
    if not os.path.exists("test-output/"):
        os.makedirs("test-output/")
        
    global safe_object
    with os.scandir(path) as sbml_folder, open('test-output/' + exel_insert + '.csv', 'w', newline='',encoding="utf-8") as exel:
        writer_csv = csv.writer(exel, delimiter=';')
        writer_csv2 = csv.DictWriter(exel, header,delimiter=';') 
        writer_csv.writerow(header)
        
        current_index=-1
        dead_entries=[560,595,909,908,1061]
        dead_entries=[]
        #only_entries=[496,497,255,560,595,909,908,1061]
        #only_entries=[496,497,255,1061]
        for entry in sbml_folder:
            #skips files before reaching start_index
            if start_index-1 >current_index:
                current_index+=1
                continue
            for i in dead_entries:
                if entry.name.endswith(str(i)+".xml"):
                    continue
            run=0
            #for i in only_entries:
            #    if entry.name.endswith(str(i)+".xml"):
            #        run=1
            #        break
            #if not run:
            #    continue
            
            if end_index <=current_index:
                break
            current_index+=1
            
            #check for file integrity
            if  ".orgml" in entry.name:
                continue
            
            
            if (entry.name.endswith(".sbml") or entry.name.endswith(".xml")) and entry.is_file():
                    #job_queueueue.append([entry.path,current_index])
                    result_line=calculate_RN([entry.path,current_index])
                    writer_csv2.writerow(result_line)
                    exel.flush()
    

def calculate_RN(input_list):
    more_timers=False
    entry,current_index=input_list
    print(current_index)
    RN=SBML_to_RN(str(entry))
    numbers_array=[0,0,0,0,0,0,0,0,0]
    for ele in RN.reactions:
        numbers=len(ele.listOfReactants)
        try:
            numbers_array[numbers]+=1
        except IndexError:
            numbers_array=numbers
            break
    
    RN.replace_inflow_by_selfreplication()
    current_SBML=Analysis(RN)
    #call of analysis 
    global database_dict
    if current_SBML.reaction_network is None:
        print(entry)
        return()
    
    
    #insert basic data    
    outputdict={}    
    outputdict["reaction_orders"]=numbers_array
    timer_all=time.time()
    
    outputdict['ID']=os.path.basename(entry)
    outputdict['name']=current_SBML.name
    
    #calculate SORs
    current_SBML.get_ERCs()
    current_SBML.all_SORs()
    SOR_o=0
    for ele in current_SBML.SORs:
        if check_ORG_LOR(ele,current_SBML.reaction_network):
            SOR_o+=1
            
    #insert data of SORs
    outputdict['number_SORs']=len(current_SBML.SORs)
    outputdict['number_SOR_O']=SOR_o
    outputdict['number_SOR_DO_only']=len(current_SBML.SORs)-SOR_o
    
    
    if hasattr(current_SBML,"abort_run"):
        outputdict["termination"]="ERC timeout"
        return(outputdict)
    
    if len(current_SBML.reaction_network.reactions)==0:
        return outputdict
    
    outputdict['number_species']=len(current_SBML.reaction_network.species)
    outputdict['number_reactions']=len(current_SBML.reaction_network.reactions)
    outputdict["number_constr"]=current_SBML.metadata["number_constr"]
    
    outputdict["Timer_ERC"]=str(current_SBML.metadata["Time_ERC"]).replace('.',',')
    outputdict["Timer_LP"]=str(current_SBML.metadata["Timer_LP"]).replace('.',',')
      
      
    
    start_MRC=time.time()
    cancel_mrc=False
    get_comp=True
    compartment_2=0
    compartment_3=0
    compartment_4=0
    compartment_5=0
    compartment_6=0
    maxnumbercomp=0
    timer1=str(time.time()-timer_all).replace('.',',')  
    
    if get_comp==True:
        for i in range(len(current_SBML.SORs)):
            if not check_ORG_LOR(current_SBML.SORs[i],current_SBML.reaction_network):
                if current_SBML.SORs[i][0]:
                    compartments=[]
                    try:
                        solution_gc=get_MCs(current_SBML.reaction_network,current_SBML.SORs[i])
                        maxnumbercomp=max(maxnumbercomp,len(solution_gc))
                        compartments=get_min_compartments(current_SBML.reaction_network, solution_gc, current_SBML.SORs[i])
                    except MemoryError:
                        cancel_mrc=True
                        break

                    #compartments=get_min_compartments(current_SBML.reaction_network, solution_gc, current_SBML.SORs[i])
                else: compartments=[]
                #if len(compartments)>2:
                 #   compartment_counter.append([i,len(compartments)])
                if len(compartments)==2:
                    compartment_2+=1
                if len(compartments)==3:
                    compartment_3+=1
                if len(compartments)==4:
                    compartment_4+=1
                if len(compartments)==5:
                    compartment_5+=1
                if len(compartments)>5:
                    compartment_6+=1
            
                if int(time.time()-start_MRC)>300:
                    cancel_mrc=True
                    print("cANCELA")
                    break
        if cancel_mrc:   
            outputdict["timer_MRC"]="timeout"
        else:
            outputdict["timer_MRC"]=str(time.time()-start_MRC).replace('.',',')



        outputdict["max_MRC"]=maxnumbercomp
        outputdict["comp2"]=compartment_2
        outputdict["comp3"]=compartment_3
        outputdict["comp4"]=compartment_4
        outputdict["comp5"]=compartment_5
        outputdict["comp6"]=compartment_6
    else:
    #outputdict["compartment_lengths"]=compartment_counter
        outputdict["max_MRC"]=""
        outputdict["comp2"]=""
        outputdict["comp3"]=""
        outputdict["comp4"]=""
        outputdict["comp5"]=""
        outputdict["comp6"]=""
        outputdict["timer_MRC"]=""
    timer2=str(time.time()-timer_all).replace('.',',')  
    #outputdict["Timer_LP_DO2"]=str(time.process_time()-start_DO2_time).replace('.',',')
    
    
    current_SBML.all_DOs()
    timer3=str(time.time()-timer_all).replace('.',',')
    outputdict["Timer_LP_DO1"]=str(current_SBML.metadata["Timer_LP"]).replace('.',',')
    
    Os=calculate_os(current_SBML)

    outputdict["number_Os"]=Os
    outputdict["number_DO"]=len(current_SBML.DOs)
      
    #if hasattr(current_SBML,"SORs_wrong"):
     #   outputdict['termination']=current_SBML.SORs_wrong
      #  timeout_counter=1
 
    if more_timers:
        
        outputdict["timer_all1"]=timer1
        
        outputdict["timer_all2"]=timer2
        
        outputdict["timer_all3"]=timer3
        outputdict["timer_all4"]=str(time.time()-timer_all).replace('.',',')
    else:
        outputdict["time_complete"]=str(time.time()-timer_all).replace('.',',')
        
    return outputdict
    


















def get_species_closure(loR,list_of_SORs):
    o_set=set()
    
    for SOR in list_of_SORs:
        if check_ORG_LOR(SOR[0],loR):
            species_set=set()
            for reaction in loR:
                if reaction.defined_name in SOR[0]:
                    species_set.update(reaction.listOfReactants)
                    species_set.update(reaction.listOfProducts)
            species_set_frozen=frozenset(species_set)
            o_set.add(species_set_frozen)
    return(o_set)

def get_species_closure2(loR,list_of_SORs):
    SOR_species=[]
    
    for SOR in list_of_SORs:
        species_set=set()
        for reaction in loR:
            if reaction.defined_name in SOR[0]:
                species_set.update(reaction.listOfReactants)
                species_set.update(reaction.listOfProducts)
        species_set_frozen=frozenset(species_set)
        SOR_species.append(species_set_frozen)
    return(SOR_species)

            
            
def calculate_number_DO(listofSpecies,SORs,LOR):
    species_min_SOR=set()
    for rea in SORs[0][0]:
        for reaction in LOR:
            if reaction.defined_name==rea:
                species_min_SOR.add(frozenset(set(reaction.listOfProducts).union(set(reaction.listOfReactants))))
    free_spec=set()
    for spec in listofSpecies:
        if spec not in species_min_SOR:
            free_spec.add(spec)
    #get ERC of species_existence
    ERC_dict1= create_closures(LOR, species=free_spec, Timer_in_sec=120)
    dict1,dict2=create_Hasse(SORs,LOR, analyze_only=True)
    #get SORs, where species can exist
    
    species_of_SORs=get_species_closure2(LOR,SORs)
            
    
    
    
    new_map={}
    combi_spec=0
    combi_spec=[]
    print(ERC_dict1)
    #für jede abh. spez
    for ele in ERC_dict1:
        if len(ERC_dict1[ele].reactions)!=0:
            print("stuff")
            print(ele)
            print(type(ERC_dict1[ele].reactions[0]))
            print(ERC_dict1[ele].reactions)
            #für jeden vertex
            for solution_index in range(len(SORs)):
                print("vertex")
                print(dict1[solution_index])
                #wenn rea1 grün
                if ERC_dict1[ele].reactions[0].defined_name in dict1[solution_index]:
                    if ele not in new_map:
                        new_map[ele]=[species_of_SORs[solution_index]]
                    else:
                        new_map[ele].append(species_of_SORs[solution_index])
                    print()
        else:
            combi_spec.append(ele)
    print(new_map)
    print(combi_spec)
    number_DOs=0
    number_DOs+=2**len(combi_spec)
    for ele in new_map:
        combi_sec_ele=combi_spec
        for opportunity in new_map[ele]:
            for species in opportunity:
                if species in combi_spec:
                    combi_sec_ele-=1
            number_DOs+=2**len(combi_spec)

