# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 14:33:11 2023

@author: linus
"""


'''create sbml data with alg from veloz'''
def wat(elem):
    return(int(elem))
from pathlib import Path
import sys
import time
import os
import io
from contextlib import redirect_stdout
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from analyse_class import Analyse
from Hasse import analyze_DO,check_ORG_LOR,check_ORG_LOR_expanded
from reaction_ERC import Reaction
import time 
from LP_compartments import get_compartments, get_min_compartments
from iterate_over_DB_multi import get_species_closure
from pandas import DataFrame
#arguments
species=[10,12,14,16,18,20]
#species=[10,12]
# reactions=[100,110,120,130,140,30,40,50,60,70,80,90]
# reactions=[x for x in range(10,51,10)]
reactions=[30,40,50,60,70,80,90,100,110,120,130,140]
reactions=[i for i in range(10,51,5)]
#reactions=[10,15,20]
plsno=time.time()
i=0
markerlist=[]
markerlist2=[]
time_data_together={}
reaction_number_data=set()
datadict={}
ERC_data_together={}
SAR_hit_list_together={}
SAR_number_list_together={}
datadict_constr_together={}
#reaction_number_data=[10,20,40,60,80,100,200]

# folders= os.listdir("C:/Users/linus/python/final_generated_sbml/10_species")
path="C:/Users/linus/python/sbml_next"
#path="C:/Users/linus/python/sbml_new2"
folders2=os.listdir(path)


# speciesq1=[10,30,50]
# speciesq1=[30,50,70]
# reactionsq2=[x for x in range(30,71,10)]
q1=-1

for species_number in folders2:
    q1+=1
    q2=-1
    subfolders=os.listdir(path+"/"+species_number)
    print(species_number)
    ERC_avg_overall=[]
    time_data_overall=[]
    
    datadict_ERC=[]
    datadict_time=[]
    datadict_SAR=[]
    datadict_SAR_number=[]
    datadict_reaction_numbers=[]
    datadict_constr=[]
    SAR_hit_list=[]
    SAR_number_list=[]
    
    
    
    subfolders.sort(key=wat)
    for folder in subfolders:
        q2+=1
        run=False
        ERC_len_data=[]
        time_data=[]
        number_constr=0
        SAR_hits=0
        SAR_hits2=0
        SAR_numbers=0
        for subdir, dirs, files in os.walk(path+"/"+species_number +"/"+ folder):

            for file in files:
                time1= time.process_time()
                x=os.path.join(subdir, file)
                current_SBML=Analyse(x) 
                current_SBML.change_network(remove_species_from_reactions="null")
                ERC_lengths= current_SBML.print_ERC_len()

                current_SBML.SARs(see_constraints=True)
                number_constr+=current_SBML.number_constr
                if len(current_SBML.SAR_solution)>1:
                    SAR_hits+=1
                    if len(current_SBML.SAR_solution)>2:
                        SAR_hits2+=1
                    #SAR_numbers+=len(current_SBML.SAR_solution)-2
                    #if len(current_SBML.SAR_solution)>5:
                    #    current_SBML.draw_hasse()
                ERC_length_avg_file=sum(ERC_lengths)/len(ERC_lengths)
                ERC_len_data.append(ERC_length_avg_file)
                run=True
                time2= time.process_time()
                time_data.append(time2-time1)
                i+=1
    ##################################################################################################################
    #create after all files
        if run==False:
            continue
        current_SBML.get_species()    
        number_reactions=len(current_SBML.reaction_network)
        reaction_number_data.add(number_reactions)
#         if number_species==0:
#             number_species=len(current_SBML.species)
#             datadict_reaction_numbers[number_species]=[]
        
        datadict_ERC.append(sum(ERC_len_data)/len(ERC_len_data)/number_reactions)
        datadict_time.append(sum(time_data)/len(time_data))
        datadict_SAR.append(SAR_hits/len(time_data))
        datadict_SAR_number.append(SAR_hits2/len(time_data))
        #datadict_SAR_number.append(SAR_numbers/len(time_data))

        datadict_constr.append(number_constr/len(time_data))
    #number_species=len(current_SBML.species)  
    """quickfix"""
    #number_species=species_number[q1]
    print("here")
    print(datadict_ERC)
    time_data_together[species_number]=datadict_time
    ERC_data_together[species_number]=datadict_ERC
    SAR_hit_list_together[species_number]=datadict_SAR
    SAR_number_list_together[species_number]=datadict_SAR_number
    datadict_constr_together[species_number]=datadict_constr
print("test")
print(ERC_data_together)
print(reaction_number_data)

reaction_number_data=sorted(reaction_number_data)
results1=DataFrame(time_data_together, index=reactions)
results2=DataFrame(ERC_data_together, index=reactions)
results3=DataFrame(SAR_hit_list_together, index=reactions)
results4=DataFrame(SAR_number_list_together,index=reactions)
results5=DataFrame(datadict_constr_together, index=reactions)
print(results5)
from matplotlib import pyplot
pyplot.figure(figsize=(6, 4), 
           dpi = 600)
sns.heatmap(results1)
plt.show()
pyplot.figure(figsize=(6, 4), 
           dpi = 600)
sns.heatmap(results2)
plt.show()
pyplot.figure(figsize=(6, 4), 
           dpi = 600)
sns.heatmap(results3)
plt.show()
pyplot.figure(figsize=(6, 4), 
           dpi = 600)
sns.heatmap(results4)
plt.show()
pyplot.figure(figsize=(6, 4), 
           dpi = 600)
sns.heatmap(results5)
# Show the plot
plt.show()
# results1= DataFrame({str(folders2[0]):time_data_together[0],str(folders2[1]):time_data_together[1]},index=reaction_number_data)
# results2= DataFrame({str(folders2[0]):ERC_data_together[0],str(folders2[1]):ERC_data_together[1]},index=reaction_number_data)
# results3= DataFrame({str(folders2[0]):SAR_hit_list_together[0],str(folders2[1]):SAR_hit_list_together[1]},index=reaction_number_data)
# results4= DataFrame({str(folders2[0]):SAR_number_list_together[0],str(folders2[1]):SAR_number_list_together[1]},index=reaction_number_data)
# fig,(results1,results2)=plt.subplots(1,2)


fig,axes= plt.subplots(nrows=3,ncols=2)
results1.interpolate(method='linear').plot(ax=axes[0,0],title="time")
results2.interpolate(method='linear').plot(ax=axes[0,1],title="respective ERC-length")
results3.interpolate(method='linear').plot(ax=axes[1,0],title="# SBML with reactive SARs")
results4.interpolate(method='linear').plot(ax=axes[1,1],title="# reactive SARS")
results5.interpolate(method='linear').plot(ax=axes[2,1],title="# constraints")
fig.tight_layout(pad=1.0)
# print("plsno")
print(time.time()-plsno)
