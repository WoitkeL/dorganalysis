# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 20:19:25 2023

@author: linus
"""

import os
import time
import csv
from Hasse import check_ORG_LOR,create_Hasse
from dorganalysis import Analysis
from reactionnetwork import create_closures, SBML_to_RN


def calculate_os(analyze):
    orgs=0
    print(analyze.SOR_map)
    for i in range(len(analyze.DOs)):
        print(analyze.DOs[i])
        print(analyze.SOR_map[analyze.DOs[i]])
        if is_DO_O(analyze.reaction_network,analyze.DOs[i],analyze.SOR_map[analyze.DOs[i]][0]):
            orgs+=1
    return(orgs)
def is_DO_O(RN,DO,SOR):
    SOR2=set()
    for rea in RN.reactions:
        if set(rea.listOfReactants).issubset(DO):
            SOR2.add(rea.defined_name)
    if SOR==SOR2:
        return(True)
    else:
        return(False)
    
def iterate_over_database(path="C:/Users/linus/python/nice_BM2/", start_index=0,end_index=3000, exelname=False):
    print(globals())
   
    header = ['ID',"name",'number_species', 'number_reactions', 'termination','number_SORs',\
              'number_SOR_O', "number_SOR_DO_only","number_DO","number_Os","number_reak_DO","numbero",\
                  "Timer_ERC","Timer_LP","Timer_LP_DO1","Timer_LP_DO2","number_constr","timer_all1","timer_all2","timer_all3"]
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
        dead_entries=[560,595,909,908]
        for entry in sbml_folder:
            #skips files before reaching start_index
            if start_index-1 >current_index:
                current_index+=1
                continue
            for i in dead_entries:
                if entry.name.endswith(str(i)+".xml"):
                    continue
            
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
    entry,current_index=input_list
    RN=SBML_to_RN(str(entry))
    current_SBML=Analysis(RN)
    timeout_counter=0
    #call of analysis 
    global database_dict
    if current_SBML.reaction_network is None:
        print(entry)
        return()
    
    timer_all=time.time()
    #insert basic data
    outputdict={}
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
      
    outputdict["timer_all1"]=str(time.time()-timer_all).replace('.',',')    
    
    dos_alter=set()
    start_DO2_time=time.process_time()
    try:
        
        for ele in current_SBML.SORs:
            dos_alter.update(tuple(current_SBML.get_DOs_of_SOR(ele)))
        outputdict["numbero"]=len(dos_alter)
    except:
        outputdict["numbero"]="error"
        #dos_alter.update(get_DOs_of_SOR(current_SBML.reaction_network, ele))
    
    outputdict["timer_all2"]=str(time.time()-timer_all).replace('.',',')    
    outputdict["Timer_LP_DO2"]=str(time.process_time()-start_DO2_time).replace('.',',')
    outputdict["numbero"]=len(dos_alter)
    current_SBML.all_DOs()
    
    outputdict["Timer_LP_DO1"]=str(current_SBML.metadata["Timer_LP"]).replace('.',',')
    Os=calculate_os(current_SBML)

    outputdict["number_Os"]=Os
    outputdict["number_DO"]=len(current_SBML.DOs)
    
    if hasattr(current_SBML,"SORs_wrong"):
        outputdict['termination']=current_SBML.SORs_wrong
        timeout_counter=1
 
    if timeout_counter!=0:
        try:
            if hasattr(current_SBML, "info_so_far"):
                outputdict["number_SOR_O"]=current_SBML.info_so_far[0]
                outputdict["number_SOR_DO_only"]=current_SBML.info_so_far[1]
        except KeyError:
            pass
    outputdict["timer_all3"]=str(time.time()-timer_all).replace('.',',')
    
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
            #combi_spec+=1
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

