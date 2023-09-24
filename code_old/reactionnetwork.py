# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:54:46 2023

@author: linus
"""

import time
import libsbml
import os
import copy
class ReactionNetwork:
    """ReactionNetwork class represents an instances of a reaction network. It can be printed and altered,
    before creating a analyse class object with it."""
    def __init__(self,list_of_reactions=[],species_set={},name="default_name"):
        """initialization with list of reaction class object as reaction-attribute and species set as
        species-attribute"""
        self.reactions=copy.copy(list_of_reactions)
        self.species=copy.copy(species_set)
        self.name=name
        
    def __str__(self):
        line_List=[["reaction","reactants","products"]]
        padding1=len("reaction")
        padding2=len("reactants")
        padding3=len("products")
        for reaction in self.reactions:
            left_alignment=""
            right_alignment=""
            #reactants are extacted
            for ele in range(len(reaction.listOfReactants)):
                left_alignment+=str(reaction.reac_stoich[ele]) + " "+ reaction.listOfReactants[ele] + "  "
                if ele<len(reaction.listOfReactants)-1:
                    left_alignment+="+ "
            #products are extacted
            for ele in range(len(reaction.listOfProducts)):
                right_alignment+=str(reaction.prod_stoich[ele]) + " "+ reaction.listOfProducts[ele] + "  "
                if ele<len(reaction.listOfProducts)-1:
                    right_alignment+="+ "
            #to ensure the output to be aligned, the print command uses the maximal length of the objects
            padding1 = max(padding1, len(reaction.defined_name))
            padding2 = max(padding2, len(left_alignment))
            padding3 = max(padding3, len(right_alignment))
            line_List.append([reaction.defined_name, left_alignment, right_alignment])
        #the terms are printed with a colon the reaction arrow
        output=""
        for a,b,c in line_List:
            output+=f'{a:<{padding1}}:  {b:<{padding2}}->  {c:<{padding3}}'+"\n"
        output+="\n"
        output+="speciesset: "+str(self.species)
        
        return(output)
    
    def copy_RN(self):
        """returns a deepcopy of the reaction network that is, also the list of species and list of reactions is cloned"""
        new_RN=ReactionNetwork(self.reactions, self.species, self.name)
        return(new_RN)
    
    
    
    def add_reaction(self,reaction):
        if isinstance(reaction, Reaction):
            self.reactions.append(reaction)
            self.species.update(reaction.listOfReactants)
            self.species.update(reaction.listOfProducts)
            
    def add_species(self,species):
        if type(species)==str:
            species={species}
        self.species.update(species)
        
    def replace_inflow_by_selfreplication(self):
        """replaces each inflow reaction of the form " -> x " by a self replication reaction " x -> 2 x" and marks this new
        reaction as mandatory that is, it must appear in at least one compartment."""
        for reaction in self.reactions:
            if reaction.always:
                reaction.listOfReactants.extend(reaction.listOfProducts)
                reaction.reac_stoich.extend([1]* len(reaction.listOfReactants))
                reaction.prod_stoich=[reaction.prod_stoich[i]+1 for i in range(len(reaction.prod_stoich))]
                
    def remove_reaction(self,reaction):
        if type(reaction)==str:
            for reaction2 in self.reactions:
                if reaction2.defined_name==reaction:
                    self.reactions.remove()
class Reaction:
    """Reaction class represents instances of reaction with all their used informations
    each reaction is initialized with the following parameters: name, reactants, products and their affiliated
    stoichiometric parameters"""
    def __init__(self, name, reactants, products, extract_stoich_rea=[], extract_stoich_prod=[]):
        """Creates an instance of a Reaction and appends this instance to the listOfReactions"""
        self.closed=True                         # properties of closure
        self.defined_name=name                   # name of raction
        self.reac_stoich=extract_stoich_rea      #stoichiometric factor of i-th reactant
        self.prod_stoich=extract_stoich_prod     #stoichiometric factor of i-th product
        self.listOfReactants=reactants           #species reference of i-th reactant
        self.listOfProducts=products             #species reference of i-th product
        self.always=(len(self.listOfReactants)==0) 
        self.reversible=False
    def __repr__(self):
        left_alignment=""
        right_alignment=""
        for ele in range(len(self.listOfReactants)):
                left_alignment+=str(self.reac_stoich[ele]) + " "+ self.listOfReactants[ele] + "  "
                if ele<len(self.listOfReactants)-1:
                    left_alignment+="+ "
            #products are extacted
        for ele in range(len(self.listOfProducts)):
            right_alignment+=str(self.prod_stoich[ele]) + " "+ self.listOfProducts[ele] + "  "
            if ele<len(self.listOfProducts)-1:
                right_alignment+="+ "
        output=self.defined_name + ":"+ left_alignment+" -> "+right_alignment
        return output
    
    def __str__(self):    
        return self.defined_name

#def generate_closure_for_species(ele,LOR):
 #   x=ERC(LOR,solospecies=ele)
 #   return(x)

#def generate_closure_for_reaction(ele,LOR):
 #   x=ERC(LOR,reaction=ele)
  #  return(x)

def generate_closure_for_species(RN):
    """returns all elementary species closures(ESC). The return value is a dictionary that maps a species to its elementary species closure,
    which represents the set of reactions that are active in the closure of the species."""
    dict_output={}
    seconds_to_timeout=180
    time_start=time.process_time()
    for specie in RN.species:
        #dict_output[specie]=generate_closure_for_species(specie,LOR)
        dict_output[specie]=elementary_species_closure(RN,solospecies=specie)
        if time_start:
            if time.process_time() > seconds_to_timeout + time_start:
                print("species closure timeout")
                return("error") 
                
    return(dict_output)

def generate_closure_for_reactions(RN):
    """returns all elementary reaction closures (ERCs). The return value is a dictionary that maps a reaction to its ERC, 
    which is the set of reactions that are active in the closure of the reactants and product species of the reaction."""
    dict_output={}
    seconds_to_timeout=180
    time_start=time.process_time()
    for ele in RN.reactions:
        dict_output[ele.defined_name]=ERC(RN,reaction=ele)
        if time_start:
            if time.process_time() > seconds_to_timeout + time_start:
                print("ERC timeout")
                return("error")
    return(dict_output)

    
def create_closures(RN, species=None, Timer_in_sec=False):
    time_start=False
    dict_output={}
    if type(Timer_in_sec)==int:
        time_start=time.process_time()
    
    if species==None:
        for ele in RN.reactions:
            #dict_output[ele.defined_name]=generate_closure_for_reaction(ele,LOR)
            dict_output[ele.defined_name]=ERC(RN,reaction=ele)
            if time_start:
                if time.process_time() > Timer_in_sec + time_start:
                    return("error")

    else:
        for specie in species:
            dict_output[specie]=elementary_species_closure(RN,solospecies=ele)
            if time_start:
                if time.process_time() > Timer_in_sec + time_start:
                    return("error")    
    return(dict_output)
    
class ERC:
    """Elementary Reaction Closures
    class object for smallest compartents to perform the reactions of the system."""

    def __init__(self, RN, reaction=[]):
        """initialization includes creating the list of species and the list of reactions of the ERC and adding the appropriate items of the starting reaction / solospecies"""
        if reaction !=[]:
            self.defined_name=reaction.defined_name
            self.reactions=[reaction]
            self.species=set()
            for element in reaction.listOfReactants:
                self.species.add(element)
            for element in reaction.listOfProducts:
                self.species.add(element)
            self.eRC_aufstellung(RN.reactions)

    def eRC_aufstellung(self, lOR_solve):
        """
        function iterates over all reactions and checks if a reaction, which is not yet part of the ERC, is supported
        by the species set
        if this is the case, the reaction and its products are added to the ERC, and the check of the reactions continues
        if the species set is checked for every reaction without adding species, the function terminates
        """
        itera= [(i, False) for i in range(len(lOR_solve))]
        checker={key: value for (key, value) in (itera)}

        last_match=(len(lOR_solve)-1)

        #infinite loop until broken by return command
        while True:
            #iterate over all reactions
            for checkreaction_index in range(len(lOR_solve)):
                #only check reactions, which are not a part already
                if not lOR_solve[checkreaction_index] in self.reactions:
                    #reaction supported
                    
                    if set(lOR_solve[checkreaction_index].listOfReactants).issubset(self.species):
                        #add reaction and species to ERC
                        
                        self.reactions.append(lOR_solve[checkreaction_index])
                        self.species.update(lOR_solve[checkreaction_index].listOfProducts)
                        
                        last_match=checkreaction_index
                        #last match is safed and loop is continued, so the loop only return after a whole loop over all reactions 
                        continue
                    
                    #also try for reversible reactions with swapped list of reactions and products
                    elif lOR_solve[checkreaction_index].reversible==True:
                        if set(lOR_solve[checkreaction_index].listOfProducts).issubset(self.species):
                            self.reactions.append(lOR_solve[checkreaction_index])
                            self.species.update(lOR_solve[checkreaction_index].listOfReactants)
                            last_match=checkreaction_index
                            continue
                    
                if checkreaction_index==last_match:
                    return
                
                
    def __str__(self):
        outp=[]
        for ele in self.reactions:
            outp.append(ele.defined_name)
        return str(outp)
        
    #def __repr__(self):
     #   for ele in self.reactions:
      #      ele.defined_name
        
class elementary_species_closure(ERC):
    def __init__(self, RN,solospecies):
        if solospecies==[]:
                self.species=set()
                self.defined_name="empty"
        else:
            self.defined_name=solospecies
            self.species={solospecies}
        self.reactions=[]
        self.eRC_aufstellung(RN.reactions)
#Es gibt 2 Möglichkeiten, die Reaktionen in das Programm einzuspeisen. Wird kein Dateipfad mit einer SBML-Datei gegeben, 
#greift es die im Programmtext manuell gegebenen Reaktionen ab.
#uses the 2 options to get the reaction network. If no path to an SBML-file is given, it takes the  manually alterable
#reaction network down below
def get_example(example=1):
    '''function to return default reaction network class:
        available networks:
            
        "generator"    
        network of 3 generators b1,b2,b3 producing the interacting species a1,a2,a3. most of the time,
        the generators can only exist in different compartments.
        
        "virus"
        network of virus infection model
        '''
    reverse_thing=[]
    listOfReactions = []
    if example==1 or example=="generator":
        listOfReactions.append(Reaction("r1",["b1"],["b1" , "a1"],[1],[1,1]))
        listOfReactions.append(Reaction("r2",["b2"],[ "b2" , "a2"],[1],[1,1]))
        listOfReactions.append(Reaction("r3",["b3"],[ "b3" , "a3"],[1],[1,1]))
        listOfReactions.append(Reaction("r4",["b2","b3"],["b3"],[1,1],[2]))
        listOfReactions.append(Reaction("r5",["b1","a2"],[],[1,1],[]))
        listOfReactions.append(Reaction("r6",["b1","b3"],["b2"],[1,1],[2]))
        listOfReactions.append(Reaction("r7",["a1","a2"],[],[1,1],[]))
        listOfReactions.append(Reaction("r8",["a1","a2","a3"],["d1"],[1,1,1],[1]))
        listOfReactions.append(Reaction("r9",["d1","b1"],["b1","a2"],[1,1],[2,1]))
        name="generator"
    if example==2 or example=="virus":
        listOfReactions.append(Reaction("r1",[],["h"],[],[1]))
        listOfReactions.append(Reaction("r2",["h","v"],["in","v"],[1,1],[1,1]))
        listOfReactions.append(Reaction("r3",["h","v"],["h","vb"],[1,1],[1,1]))
        listOfReactions.append(Reaction("r4",["vb"],["m"],[1],[1]))
        listOfReactions.append(Reaction("r5",["p"],["v"],[1],[1]))
        listOfReactions.append(Reaction("r6",["m"],["m" , "s"],[1],[1,1]))
        listOfReactions.append(Reaction("r7",["s"],[],[1],[]))
        listOfReactions.append(Reaction("r8",["in","s"],["s"],[1,1],[1]))
        listOfReactions.append(Reaction("r9",["v","s"],["s"],[1,1],[1]))
        name="virus"
    speciesSet=set()
    for reaction in listOfReactions:
        speciesSet.update(set(reaction.listOfReactants))
        speciesSet.update(set(reaction.listOfProducts))
    x=ReactionNetwork(listOfReactions,speciesSet,name)
    return(x)
def SBML_to_RN(path, consider_reverse=True, consider_constant=False, consider_init_ammount=False,alternative_reverse=False):
        
        listOfReactions = []
        #SBML-Reader from LibSBML-Package, is initialized
        reader = libsbml.SBMLReader()
        
        # Check for various errors
        if not os.path.isfile(path):
            print("no file in path found")
            return()
        if reader == None:
            print("no object created")
            
        # Read SBML model from file
        doc_extract = reader.readSBMLFromFile(path)
        
        # Check for errors in SBML file
        if doc_extract.getNumErrors() > 0:
            if doc_extract.getError(0).getErrorId() == libsbml.XMLFileUnreadable:
                print("XMLFileUnreadable")
            elif doc_extract.getError(0).getErrorId() == libsbml.XMLFileOperationError:
                print("XMLFileOperationError")
                
                

        # Extract model from SBML file
        model_extr=doc_extract.getModel()
        if model_extr is None:
            return(None, None)
        # Get model name
        model_name="no name"
        model_name= model_extr.getName()
        
        # reactionlist is extracted from model:
        listrea_extract=model_extr.getListOfReactions()
        #species list is extracted from model:
        reverse=[]
        
        # Loop through all reactions
        for i in range (len(listrea_extract)):
            # Get ID of reaction
            name=listrea_extract.get(i).getId()
            
            # Get lists of reactants and products         
            extract_rea_list=listrea_extract.get(i).getListOfReactants()
            extract_pro_list=listrea_extract.get(i).getListOfProducts()
            extract_rea_list_species=[]
            extract_pro_list_species=[]
            extract_rea_list_species_stoich=[]
            extract_pro_list_species_stoich=[]
            
            # Extract species and their stoichiometric factors from reactants list
            for j in range (len(extract_rea_list)):
                extract_rea_list_species.append(extract_rea_list[j].getSpecies())
                extract_rea_list_species_stoich.append(extract_rea_list[j].getStoichiometry())
            
            # Extract species and their stoichiometric factors from products list  
            for j in range (len(extract_pro_list)):
                extract_pro_list_species.append(extract_pro_list[j].getSpecies())
                extract_pro_list_species_stoich.append(extract_pro_list[j].getStoichiometry())

            # memorize reversible reactions       
            if listrea_extract.get(i).getReversible():
                reverse.append(name)
                
            # Initialize reaction with Reaction class        
            x= Reaction(name,extract_rea_list_species,extract_pro_list_species,extract_rea_list_species_stoich,extract_pro_list_species_stoich)
            listOfReactions.append(x)
            
        # Check for reversible reactions and add corresponding reverse reactions if necessary    
        if consider_reverse:
            # loop over all reactions
            for i in listOfReactions:
                if i.defined_name in reverse:
                    has_rea=0
                    
                    # check if reversed reaction already exists in reaction network
                    for j in listOfReactions:
                        if i==j:
                            continue
                        if set(j.listOfReactants)==set(i.listOfProducts):
                            if set(j.listOfProducts)==set(i.listOfReactants):
                                first_works=1
                                for k in range(len(i.listOfReactants)):
                                    first_works=0
                                    ind=j.listOfProducts.index(i.listOfReactants[k])
                                    if i.reac_stoich[k]==j.prod_stoich[ind]:
                                        first_works=1
                                    else: break
                                if first_works==1:
                                    if len(i.listOfProducts)==0:
                                        has_rea=1
                                    for k in range(len(i.listOfProducts)):
                                        ind=j.listOfReactants.index(i.listOfProducts[k])
                                        if i.prod_stoich[k]==j.reac_stoich[ind]:
                                            has_rea=1
                                            print("skip_reverse")
                                            
                                            
                    if not has_rea:
                        # add reversed reaction 
                        listOfReactions.append(Reaction(i.defined_name+str("_reverse"),i.listOfProducts,i.listOfReactants,i.prod_stoich,i.reac_stoich))
        
        
        # get the list of species in the model
        listspec_extract=model_extr.getListOfSpecies()
         
        species_set=[]
        for i in range(len(listspec_extract)):
            name_species=listspec_extract.get(i).getId()
            species_set.append(name_species)
            
            if consider_init_ammount and listspec_extract.get(i).getInitialAmount()>0:
                name_reaction=str(name_species)+str("init_ammount_inflow")
                listOfReactions.append(Reaction(name_reaction,[str(name_species)],[str(name_species)],[1],[2]))
                listOfReactions.append(Reaction(str(name_reaction+str(2)),[str(name_species)],[],[1],[0]))
                
        # return the ReactionNetwork object        
        model_name=model_extr.getName()   
        x=ReactionNetwork(listOfReactions, species_set,model_name)
        return(x)
    
    
def erc_to_matrix2(erc_dict):
    i=len(erc_dict)
    bool_matrix= [[False] * i for j in range(i)]
    j=0
    index_to_reac={}
    reac_to_index={}
    for key in erc_dict:
        index_to_reac[j]=key
        reac_to_index[key]=j
        j+=1
    for key in erc_dict:
        for i in erc_dict[key].reactions:
            bool_matrix[reac_to_index[key]][reac_to_index[i.defined_name]]=True
    return(bool_matrix, index_to_reac)

def sort_second(elem):
    return elem[1]

def ERC_meets_transitivity2(ERC_dict):
    bool_ma,index_to_reac= erc_to_matrix(ERC_dict)
    remover=[]
    for i in range(len(bool_ma)):
        for j in range(len(bool_ma)):
            if bool_ma[i][j]and i!=j:
                for k in range(len(bool_ma)):
                    if (bool_ma[j][k]and j!=k):
                        bool_ma[i][k] = False
                        if i!=k:
                            remover.append([i,k])
    remover.sort(key=sort_second)
    remover.reverse()
    for ele in remover: 
        oko=ERC_dict[index_to_reac[ele[0]]].reactions
        for i in oko:
            if i.defined_name==index_to_reac[ele[1]]:
                oko2=ERC_dict[index_to_reac[ele[0]]].reactions.index(i)
                del ERC_dict[index_to_reac[ele[0]]].reactions[oko2]
    
    return(ERC_dict)






def erc_to_matrix(erc_dict):
    '''transforms ERC into a boolean matrix, with a True value in [i,j] if the ERC of reaction i contains j '''
    bool_matrix = [[False] * len(erc_dict) for j in range(len(erc_dict))]
    #create both maps to encode reaction to integer, the decoding function is returned with the matrix
    index_to_reac = {i: key for i, key in enumerate(erc_dict)}
    reac_to_index = {key: i for i, key in enumerate(erc_dict)}
    
    #iteration through ERCs
    for key in erc_dict:
        for i in erc_dict[key].reactions:
            bool_matrix[reac_to_index[key]][reac_to_index[i.defined_name]] = True
    #return the matrix as well as the decoding function of indexes     
    return bool_matrix, index_to_reac

def ERC_meets_transitivity(ERC_dict):
    ''' takes a dictionary of EC numbers (ERC_dict) as input and returns the same dictionary with transitive relationships removed.'''
    #creating matrix as well as a dictionary, which maps the indexes of the matrix to the original reaction
    bool_ma, index_to_reac = erc_to_matrix(ERC_dict)
    
    #check for transitivity
    # for a set of relation a->b, b->c we need to check for a!=b,
    # b!=c and a!=c before declaring that a is not -> c 
    for i in range(len(bool_ma)):
        for j in range(len(bool_ma)):
            if bool_ma[i][j] and i != j:
                for k in range(len(bool_ma)):
                    if bool_ma[j][k] and j != k and i != k:
                        bool_ma[i][k] = False
    # transform matrix back to erc dict
    for i in range(len(bool_ma)):
        for j in range(len(bool_ma)):
            if bool_ma[i][j] and i != j:
                ERC_dict[index_to_reac[i]].reactions = [r for r in ERC_dict[index_to_reac[i]].reactions if r.defined_name != index_to_reac[j]]
    return ERC_dict
