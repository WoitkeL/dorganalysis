# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:14:02 2023

@author: linus
"""

"""HASSE"""
"""HASSE"""
"""HASSE"""

from LP_compartments import get_min_compartments,get_MCs
import regex as re
from graphviz import Digraph, Source
import copy

def check_ORG_LOR(SOR,RN):
    """function to check a subset of reactions for closure and therefore the attribute of being DO lvl 1 (being a O)"""
    speciesset=get_species_of_SOR(SOR,RN.reactions)
      
    #check for support of inactive reactions
    for reaction in RN.reactions:
        if reaction.defined_name not in SOR and set(reaction.listOfReactants).issubset(speciesset):
            return(False)
    return(True)
def get_species_of_SOR(knot_reactions,LOR):
    LOR_dict={}
    for reac in LOR:
        LOR_dict[reac.defined_name]=set(reac.listOfReactants+ reac.listOfProducts) 
    speciesset=set()
    for reaction in knot_reactions:
            speciesset.update(LOR_dict[reaction])
    return(speciesset)


def analyze_DO(x,y):
    niceDOs=[]
    niceOs=[]
    check_if_Org_dict={}
    outputdict={}
    #get dictionaries for new active reactions and new active overproduction
    #dict1 conatains all reactions, which did not appear in subsets of the SOR
    dict1,dict2= create_Hasse(x,y,analyze_only=True)
    
    new_reactions=0
    new_reactions_O=0
    equal_DO_DO=[0,0]
    equal_DO=[0,0]
    smaller_DO=[0,0]
    DO_large=0
    new_reactions_O=0
    check_if_spec_new={}
    check_if_Dist_Org_dict={}
    for i in range(len(x)):
        new_reac_bool=0
        
        if bool(dict1[i])==True:
            new_reac_bool=1
            
            
        #check if DO is also ORG
        #bool_result=check_ORG_LOR_expanded(x[i][0],y)
        bool_result=check_ORG_LOR(x[i][0],y)
        speciesout=get_species_of_SOR()
        speciesout=frozenset(speciesout)
        #for every org
        if i==0 and bool_result==False:
            DO_large=1
        if bool_result==True:
            
            niceOs.append(x[i][0])
            if new_reac_bool:
                new_reactions_O+=1
            #saves species of SOR belonging to org as key, and 
            #value determines if SOR of org has new rea
            if bool(dict1[i])==True:
                check_if_Org_dict[speciesout]=1
            else:
                check_if_Org_dict[speciesout]=0
        else:
            if new_reac_bool:
                new_reactions+=1
            #safe DOs which are not Os
            niceDOs.append(x[i][0])
            if speciesout in check_if_Dist_Org_dict:
                if check_if_Dist_Org_dict[speciesout]:
                    equal_DO_DO[1]+=1
                else:
                    equal_DO_DO[0]+=1
            #saves species of SOR belonging to org as key and 
            #value determines if SOR of org has new rea
            if bool(dict1[i])==True:
                check_if_Dist_Org_dict[speciesout]=1
            else:
                check_if_Dist_Org_dict[speciesout]=0
            
            #check for SOR of DO, if the species corresponds to SOR for DO or O respectively
            if speciesout in check_if_Org_dict.keys():
                if check_if_Org_dict[speciesout]:
                    equal_DO[1]+=1
                else:
                    equal_DO[0]+=1
            else:
                for species in check_if_Org_dict.keys():
                    if species in speciesout: 
                        if check_if_Org_dict[speciesout]:
                            smaller_DO[1]+=1
                        else:
                            smaller_DO[0]+=1
        
        
                    
        if bool(dict1[i])==True:
            if bool_result==True:
                check_if_spec_new[speciesout]="O"
            else:
                check_if_spec_new[speciesout]="DO"
    #    'number_SOR_O', "number_SOR_DO_only", "new_reactions","s(DO)=s(O)", "s(DO)=s(DO)","DO_is_largest","Os","unique_DOs",
     
    outputdict["number_SOR_O"]=len(niceOs)
    outputdict["number_SOR_DO_only"]=len(niceDOs)
    outputdict["new_reactions"]=new_reactions
    outputdict["new_reactions_O"]=new_reactions_O
    outputdict["s(DO)=s(O)"]=str(equal_DO[0])+"/"+str(equal_DO[1])
    outputdict["s(DO)=s(DO)"]=str(equal_DO_DO[0])+"/"+str(equal_DO_DO[1])
    outputdict["DO_is_largest"]=str(DO_large)
    outputdict["unique_DOs"]=niceDOs
    outputdict["Os"]=niceOs
    
    
    
    return(outputdict)

def sorted_list_subset_check(list1, list2):
    """simple function to check if first list is subset of secondlist, when bot are sorted"""
    x=0
    #for each element, we check if it matches an element of the other list
    for ele in list1:
        
        while True:
            try:
                if ele==list2[x]:
    #increase index in list2 if either match or no match
                    x+=1
    #if match break loop and therefore check for next element
                    break
                else:
                    x+=1
    #indexerror shows if list2 ends and therefor list1 isnt subset
            except IndexError:
                return(False)
    #if correct loop through all elements, return true
    return(True)

def get_highlighted_text(text_highlight, whole_text):
    """function to generate highlighted text by giving whole knot name and new/to be highlighted terms"""
    non_highlight_text = set(whole_text)-set(text_highlight)
    #goes over all cases of name1 or name 2 being empty
    #uses html as language in digraph file
    
    if text_highlight:
        #sort text_highlight
        text_highlight=list(text_highlight)
        text_highlight.sort()
        
        #use html to color text, default font color='green3'
        text_highlight="<font color='green3'>"+'%s'%', '.join(map(str, text_highlight))+"</font>"
    else: 
        text_highlight=""
        
    # add comma if both text exist
    if text_highlight and non_highlight_text:
        text_highlight=text_highlight+", "
        
    if non_highlight_text:
        #sort non_highlight_text
        non_highlight_text=list(non_highlight_text)
        non_highlight_text.sort()
        
        #add rest of label with html syntax
        non_highlight_text='%s'%', '.join(map(str, non_highlight_text))
    else:
        non_highlight_text=""
    
    #fuse both parts
    name_output= text_highlight+str(non_highlight_text)
    return(name_output)

def draw_hasse(label_highlight_reaction_dict,dot,list_of_solutions,analyze, DOs=False, shortform=False, show_species=False, only_show_new=False,show_compartments="length", second_value=False, pdf=True, display_bool=True):
    """function to handle visual parameters and transform data of hasse diagram into graphviz file."""
    RN=analyze.reaction_network
    comp=[]
    namedict={}
    LOR=RN.reactions
    
    if DOs==True:
        SORs=[analyze.SOR_map[do_ele][-1] for do_ele in list_of_solutions]
        
    else:
        SORs=list_of_solutions
        
    if shortform:
        short_form_knot={}    
        for i in range(len(LOR)):
            short_form_knot[LOR[i].defined_name]="r"+str(i+1)
    for i in range(len(list_of_solutions)):
        
        #initialize loop variables
        addinfo=""
        name="" 
        prod_name=""
        start=""
        
        #change names if shortform is active
        if shortform:  
            name1= [short_form_knot[element] for element in label_highlight_reaction_dict[i]]
            name2= [short_form_knot[element] for element in list_of_solutions[i]]
        else:
            name1= label_highlight_reaction_dict[i]
            name2= list_of_solutions[i]
        
        if only_show_new:
            name= ','.join(map(str,name1))
        else:
            start="<"
            name= get_highlighted_text(name1,name2)
            if name1==[] and name2==[]:
                name="none"
            name=start+name    
        if prod_name and second_value==True:
                name=name+" R: "+prod_name
        if start: name=name+"<br/>" 
        else: name=name+"\n"
        
        #handle species output
        if show_species:
            LOR_spec_dict={}
            for reac in LOR:
                LOR_spec_dict[reac.defined_name]=set(reac.listOfReactants+ reac.listOfProducts)
            speciesset=set()
            for reaction in list_of_solutions[i]:
                speciesset.update(set(LOR_spec_dict[reaction]))
            speciesset_ordered=list(speciesset)
            speciesset_ordered.sort()
            addinfo=addinfo+"species: "+'%s'%', '.join(map(str, speciesset_ordered))
            if start: addinfo=addinfo+"<br/>" 
            else: addinfo=addinfo+"\n"

####################################################################
#compartments

        if DOs==True and show_compartments and list_of_solutions[i]:
            MCs=get_MCs(RN,SORs[i], list_of_solutions[i])
            comp=get_min_compartments(RN, MCs,SORs[i], list_of_solutions[i])
            
        if show_compartments and list_of_solutions[i] and DOs==False:
            MCs=get_MCs(RN,SORs[i])
            comp=get_min_compartments(RN, MCs,list_of_solutions[i])
        
        if not list_of_solutions[i]:
            comp={}
        if show_compartments==True:
            comp=[list(comp_ele) for comp_ele in comp ]
            for comp_ele in comp:
                comp_ele.sort()
            comp1= ['%s'%', '.join(map(str, comp_ele))for comp_ele in comp]
            addinfo=addinfo+"{"+'%s'%' | '.join(map(str, comp1))+"}"
            #if start: addinfo=addinfo+"<br>" else: addinfo=addinfo+"\n"
        elif show_compartments=="length":
            addinfo=addinfo+"#"+ str(len(comp))
      
        name=name+addinfo
        if re.search("\<",name):
            name= name+ ">"
        namedict['node'+str(i)]=name
    
    for i in range(len(list_of_solutions)-1,-1,-1):

        is_o=True
        if DOs==True:
            true_rea_DO=[]
            knot_reac=SORs[1]
            for r in LOR:
                if set(r.listOfReactants).issubset(set(list_of_solutions[i])):
                    true_rea_DO.append(r)
                            
            for reac in true_rea_DO:
                if reac.defined_name not in knot_reac:
                    is_o=False
        else:
            knot_reac=list_of_solutions[i]
        
        if check_ORG_LOR(knot_reac,RN) and is_o:
            dot.node('node'+str(i), namedict['node'+str(i)], shape='box')
            
        else:
            dot.node('node'+str(i), namedict['node'+str(i)])
   
    line=dot.source
    line=[e+"]" for e in line.split("]")]
    line[-1]= re.sub("]","", line[-1] )
    for i in range(len(line)):
        if re.search("\<",line[i]):
            line[i]= re.sub("\"","", line[i] )
          
    if shortform:
        shortform_label_string="reaction: shortform \n"
        for ele in short_form_knot:
            shortform_label_string=shortform_label_string+ele+": "+short_form_knot[ele]+"\n"    
    addstring_merged = '\n'.join(map(str, line))
    
    dot = Source(addstring_merged)
    if pdf==True:
        dot.render('test-output/Hasse_Diagram.gv', view=True)  
    #only works in jupyter notebook
    display_bool=True
    if display_bool==True:
        pass
        #display(dot)

    return(dot)
    
def create_Hasse(list_of_solutions,analyze, DOs=False, shortform=False, show_species=False, only_show_new=False,show_compartments="length", second_value=False, analyze_only=False, pdf=True):
    RN=analyze.reaction_network
    LOR=RN.reactions
    max_knot=set(list_of_solutions[0])
    #creating inverse list of soulutions, to ensure existence of intersection of data. since orginial data is complete for union.
    inverse_list=[(max_knot.difference(set(list_of_solutions[i]))) for i in range(len(list_of_solutions))]
    #dictionary maps inverse elements to index of original vertex
    knot_to_index={}
    
    #extralists for special marking of new reactions/production
    label_highlight_reaction_dict={}
    
    #create elements for edge algorithm
    facedict={}
    borderlist=[]
    candidatelist=[]
    
    #create objects for hasse
    dot = Digraph(strict=True,comment='Hasse_Diagram')
    nodedict={}
    Os_dict={}
    species_closure={}
    
    #setup of objects, which are changed before looping over its knot
    for i in range(len(list_of_solutions)):
        
        knot_to_index[frozenset(inverse_list[i])]=i
        label_highlight_reaction_dict[i]=list_of_solutions[i][:]
        
        nodedict[frozenset(inverse_list[i])]=i
        
        if analyze_only:
            if check_ORG_LOR(list_of_solutions[i], RN):
                Os_dict[frozenset(list_of_solutions[i])]=True
            else: Os_dict[frozenset(list_of_solutions[i])]=False
        if DOs==False:
            species_closure[i]=get_species_of_SOR(list_of_solutions[i], LOR)


    #initialization of first knot
    borderlist.append(inverse_list[0])
    facedict[frozenset(inverse_list[0])]=[]

    for i in range(1,len(list_of_solutions)):
        
        #create candidates by intersection of border and current knot
        candidatelist=[]
        for borderelement in borderlist:
            candidatelist.append(set(inverse_list[i]) & (set(borderelement)))
            
    
        #create dictionary for memory of faces(all elements, which it is a subset of)
        facedict[frozenset(inverse_list[i])]=[]    
        
        #check if candidates have larger elements in face (elements smaller than knot and larger than candidate)
        for candidate in candidatelist:
            candidate_froz= frozenset(candidate)

            for faceobject in facedict[candidate_froz]:
                if faceobject.issubset(inverse_list[i]):
                    break
            else: 
                #create edge with specific handling if SORs are handled and 2 SORs have the same speciesset
                if DOs==False:
                    if species_closure[nodedict[candidate_froz]]==species_closure[i]:
                        dot.edge( str('node'+str(nodedict[candidate_froz])),str('node'+str(nodedict[frozenset(tuple(inverse_list[i]),)])), arrowhead = "none",color="red")
                    else:
                        dot.edge( str('node'+str(nodedict[candidate_froz])),str('node'+str(nodedict[frozenset(tuple(inverse_list[i]),)])), arrowhead = "none")
                else:
                    dot.edge( str('node'+str(nodedict[candidate_froz])),str('node'+str(nodedict[frozenset(tuple(inverse_list[i]),)])), arrowhead = "none")
                
                j=knot_to_index[candidate_froz]
                
                #save new reactions
                label_highlight_reaction_dict[j]= [e for e in label_highlight_reaction_dict[j] if e not in list_of_solutions[i]]
                
                #add knot to face of candidate
                add_element_to_facedict=inverse_list[i].copy()
                extract=facedict[candidate_froz]
                extract.append(frozenset(tuple(add_element_to_facedict)))
                facedict[candidate_froz]=extract
                #remove candidate from border
                #these candidates dont have to be in the border if there they are subsets of bigger knots in the border
                try:
                    borderlist.remove(candidate)
                except ValueError:
                    pass


        borderlist.append(inverse_list[i])
    if analyze_only:
        return(label_highlight_reaction_dict)
    else:
        x=draw_hasse(label_highlight_reaction_dict,dot,list_of_solutions,analyze,DOs, shortform, show_species, only_show_new, show_compartments, second_value, pdf)
        return(x)
def create_Hasse_naive(list_of_solutions,RN, DOs=False, shortform=False, show_species=False, only_show_new=False,show_compartments="length", second_value=False,analyze_only=False, pdf=True):
    dot = Digraph(strict=True,comment='Hasse_Diagram')
    #dot.attr( size='8,5')
    cover={}
    list_of_solutions.reverse()
    label_highlight_overprod_dict={}
    label_highlight_reaction_dict={}
    index_reactions_dict={}
    ################################################################################################################
    #naming
    list_of_solutions_empty=copy.deepcopy(list_of_solutions)
    if list_of_solutions[0][0]==[]:
        [s[0].append("empty_reaction") for s in list_of_solutions_empty]
    else:
        list_of_solutions_empty=list_of_solutions
    if shortform:
        short_form_knot={}    
        list_of_reactions=list_of_solutions[-1][0]
        for i in range(len(list_of_reactions)):
            short_form_knot[list_of_reactions[i]]="R"+str(i)
    ##########################################################################################################
    #setup first elements
    for i in range(len(list_of_solutions)):
        cover[i]=[]
        #umwandlung in frozenset um für dictionary vorzubereiten
        index_reactions_dict[i]=list_of_solutions_empty[i][0]
    cover[0]=[0]
    label_highlight_reaction_dict[0]=list_of_solutions[0][0]
    label_highlight_overprod_dict[0]=list_of_solutions[0][1]
##################################################################################################
#loop over all knots
    for i in range (1, len(list_of_solutions)):
        reactions_to_check=list_of_solutions[i][0][:]
        production_to_check=list_of_solutions[i][1][:]
        checkdic={}
        for j in range(i):
            checkdic[j]=""
        startcheck=i
        while True:
            for k in range(startcheck,-2,-1):
                if k in checkdic:
                    check_index=k
                    break
            
            startcheck=check_index
            if sorted_list_subset_check(list_of_solutions_empty[check_index][0],list_of_solutions_empty[i][0]):
                dot.edge("node"+str(i),"node"+str(check_index),  arrowhead = "none")
                reactions_to_check = [e for e in reactions_to_check if e not in list_of_solutions[check_index][0]]
                production_to_check = [e for e in production_to_check if e not in list_of_solutions[check_index][1]]
                
                #add all elements of linked knots to cover of current knot
                for ele in cover[check_index]:
                    if ele not in cover[i]:
                        cover[i]=cover[i]+[ele]
    
                for del_element in sorted(cover[i], reverse=True):
                    checkdic.pop(del_element, None)
            else:
                checkdic.pop(check_index, None)
            if len(checkdic)==0:
                cover[i]=cover[i]+[i]
                break
        label_highlight_reaction_dict[i]=reactions_to_check
        label_highlight_overprod_dict[i]=production_to_check
    if analyze_only:
        return(label_highlight_reaction_dict,label_highlight_overprod_dict)
    else:
   # draw_hasse(label_highlight_reaction_dict,label_highlight_overprod_dict,dot,list_of_solutions,LOR,DOs, shortform, show_species, show_new, show_compartments, overproduction)
        x=draw_hasse(label_highlight_reaction_dict,label_highlight_overprod_dict,dot,list_of_solutions,RN,DOs, shortform, show_species, only_show_new, show_compartments, second_value, pdf)
        return(x)