# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:41:31 2023

@author: linus
"""
import os
from reactionnetwork import  create_closures,generate_closure_for_species,generate_closure_for_reactions,ERC_meets_transitivity2,ReactionNetwork
from setup_LP_DO_final import basicLP
from Hasse import create_Hasse_naive, analyze_DO, create_Hasse
from LP_compartments import get_min_compartments,get_MCs
import itertools
import time

class Analysis:
    """class to handle all properties of a reaction network analysis. 
    It saves any information in class attributes and pipes them into their advancing/ processing functions. 
    """ 
    def __init__(self, RN):
        self.reaction_network=ReactionNetwork(RN.reactions, RN.species, RN.name)
        self.name=RN.name
        self.metadata={}
        

        
        
    def get_ERCs(self):
        """Returns the transitive reduced dictionary of reaction-ERC pairs."""
        if not hasattr(self, "ERCs"):
            time_start=time.process_time()
            ERC_dict_complete= generate_closure_for_reactions(self.reaction_network)
            if ERC_dict_complete=="error":
                print("timeout")
                self.metadata["abort_run"]=True
                return()
            #self.ERCs_complete=ERC_dict.copy()
            ERC_dict_reduced=ERC_meets_transitivity2(ERC_dict_complete)
            self.metadata["Time_ERC"]=time.process_time()-time_start
            self.ERCs=ERC_dict_reduced
        return(self.ERCs)
    
    
    
    def all_SORs(self):
        """returns a list of the sets of organizational reactions (SORs). THe list is sorted from the smallest to the largest SOR.
        A single SOR is a set of reaction names."""
        self.SORs=basicLP(self)
        return(self.SORs)
        
    
    def all_DOs(self):
        """returns a sorted list of all DOs, where a DO is a set of species. also the function creates the attribute SOR_map,
        which maps a DO to a list of all possible SORs, the species set can perform in a DO."""
        output=basicLP(self,DO=True)
        only_DO=set()
        DO_SOR_map={}
        for ele in output:
            only_DO.add(ele[0])
        only_DO=[ele for ele in only_DO]
        only_DO.sort(key=lambda x: len(x), reverse=True)
        for DO in only_DO:
            DO_SOR_map[DO]=[]
        for ele in output:
            DO_SOR_map[ele[0]].append(ele[1])
        for ele in DO_SOR_map:
            DO_SOR_map[ele].sort(key=lambda x: len(x))
        self.DOs=only_DO
        self.SOR_map=DO_SOR_map
        return only_DO
    
    
    def largest_SOR(self):
        """returns the largest SOR, which is a set of reactions."""
        self.all_SORs()
        return(self.SORs[-1])
        
    def largest_DO(self):
        """returns the largest DO, which is a set of species. also the function creates the attribute SOR_map,
        which maps a DO to a list of all possible SORs, the species set can perform in a DO."""
        self.all_DOs()
        return(self.DOs[-1])
    
    def get_SORs_of_DO(self,DO):
        """returns all SORs that can be performed by the DO. input is a DO as a sorted tuple of species. 
        Alternative access via class attribute .SOR_map[DO]"""
        if not hasattr(self, "SOR_map"):
            self.all_DOs()
        return(self.SOR_map[DO])
    
    
    def get_DOs_of_SOR(self,SOR,limit=18):
        """returns all DOs that can be performed by the SOR. Input is a SOR as a sorted tuple of reactions. 
        Alternative access via class attribute .SOR_map[DO]"""
        species_min_SOR=set()
        RN=self.reaction_network
        # calculate minimal species set of DO by iterating over reactants and products of each reaction of the SOR.
        for rea in SOR:
            for reaction in RN.reactions:
                if reaction.defined_name==rea:
                    species_min_SOR.update(set(reaction.listOfProducts).union(set(reaction.listOfReactants)))
        
        # collect all species, which are not in the minimal species set and can therefore be added
        free_spec=set()
        for spec in RN.species:
            if spec not in species_min_SOR:
                free_spec.add(spec)
            
        species_closures= generate_closure_for_species(RN)
        
        #these species can only be added as non reactive species, if the the species do not trigger a reaction as a sole 
        #species in a compartment
        species_combis=set()
        for species in free_spec:
            if len(species_closures[species].reactions)==0:
                species_combis.add(species)
        combinations2 = []
        
        #a certain limit prevents computational overloading
        if len(species_combis)<limit:
            newl=list(species_min_SOR)
            newl.sort()
            combinations2.append(tuple(newl))
            #for each length of compartments
            for i in range(1,len(species_combis)+1):
                #for each combination of the given length
                for ele in itertools.combinations(species_combis, i):
                    new_ele=list(set(ele).union(species_min_SOR))
                    new_ele.sort()
                    combinations2.append(tuple(new_ele))
            combinations2.reverse()
                
        else:
            print("WARNING: computation canceled. limit of SORs reached.")
            print("current limit: "+ 2^limit)
            print("limit can be extended by chaning parameter\"limit\" to the logarithm to the base 2 of the desired size. default = 18 ")
            return()
        return(combinations2)

        
    def draw_DOs(self, solution=None, shortform=False, show_compartments="length", pdf=True,use_naive=False):
        
        #options={shortform:False, show_compartments:"length", pdf:True}
        options={"shortform":shortform,  "show_compartments":show_compartments, "pdf":pdf, "DOs":True}
        if solution==None:
            if hasattr(self,"DOs"):
                solution=self.DOs
            else:
                self.all_DOs()
                solution=self.DOs
                
        if use_naive:
            dot_object=create_Hasse_naive(solution, self,**options)
        else:
            dot_object=create_Hasse(solution, self, **options)
        return(dot_object)
    def draw_SORs(self, solution=None, shortform=False, show_species=False, show_compartments="length", pdf=True,use_naive=False):
        
        options={"shortform":shortform, "show_species":show_species, "show_compartments":show_compartments, "pdf":pdf}
        if solution==None:
            if hasattr(self,"SORs"):
                solution=self.SORs
            else:
                self.all_SORs()
                solution=self.SORs
                
        if use_naive:
            dot_object=create_Hasse_naive(solution, self,**options)
        else:
            dot_object=create_Hasse(solution, self, **options)
        return(dot_object)
    """  
    def evaluate_solution(self, solution_set="default", output_style="exel"):
        
        x=0
        try:
            bla=getattr(self,solution_set)
            x=analyze_DO(bla, self.reaction_network)
        except KeyError:
            print()
        return(x)"""
            
    def print_RN(self):
        line_List=[["reaction","reactants","products"]]
        padding1=len("reaction")
        padding2=len("reactants")
        padding3=len("products")
        for reaction in self.reaction_network.reactions:
            left_alignment=""
            right_alignment=""
            for ele in range(len(reaction.listOfReactants)):
                left_alignment+=str(reaction.reac_stoich[ele]) + " "+ reaction.listOfReactants[ele] + "  "
                if ele<len(reaction.listOfReactants)-1:
                    left_alignment+="+ "
            for ele in range(len(reaction.listOfProducts)):
                right_alignment+=str(reaction.prod_stoich[ele]) + " "+ reaction.listOfProducts[ele] + "  "
                if ele<len(reaction.listOfProducts)-1:
                    right_alignment+="+ "
                    
            padding1 = max(padding1, len(reaction.defined_name))
            padding2 = max(padding2, len(left_alignment))
            padding3 = max(padding3, len(right_alignment))
            line_List.append([reaction.defined_name, left_alignment, right_alignment])
        for a,b,c in line_List:
            print(f'{a:<{padding1}}:  {b:<{padding2}}->  {c:<{padding3}}')
            print()
            
    def print_ERC(self):
        """draw function to print Analysis.ERCs properly"""
        try:
            self.ERCs= self.ERCs
        except AttributeError:
            self.ERCs= create_closures(self.reaction_network, Timer_in_sec=120)
        if type(self.ERCs)==str:
                print("timeout ERC")
                return
        print("ERCs:", end=" ")
        [print(key,':',[ele.defined_name for ele in self.ERCs[key].reactions ]) for key in self.ERCs.keys()]
        print()
        
    def print_ERC_len(self):
        try:
            self.ERCs= self.ERCs
        except AttributeError:
            self.ERCs= create_closures(self.reaction_network, Timer_in_sec=120)
            if type(self.ERCs)==str:
                print("timeout ERC")
                return
            
        output_array=[len(self.ERCs[key].reactions)-1 for key in self.ERCs.keys()]
        return(output_array)


    def get_compartmentalization_of_SOR_DO_pair(self, SOR, species=""):
        """ function links function for MC calculation and calculation of a compartmentalization with minimal number of 
        compartments. Input is SOR and a DO. the minimal DO of the SOR is the default DO."""
        MCs=get_MCs(self.reaction_network, SOR, species)
        self.MCs=MCs
        x=get_min_compartments(self.reaction_network, MCs, SOR, species)
        return(x)
