 # -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 20:48:14 2023

@author: linus
"""

from iterate_over_DB import iterate_over_database
if __name__ == '__main__':
    safe_object=[]
    database_dict={}
    iterate_over_database(path="C:/Users/linus/python/BioModels-files/",exelname="biomodels_paper14",start_index=500)
    #iterate_over_database(path="C:/Users/linus/python/BioModels-files/",use_comma=True,exelname="biomodels_paper6",start_index=200, end_index=400)
    #iterate_over_database(path="C:/Users/linus/python/BioModels-files/",use_comma=True,exelname="Kaleta_1_0_0_new",end_index=184)
    
    #iterate_over_database3(path="C:/Users/linus/python/BioModels-files/",use_comma=True,exelname="test_new3",end_index=200)
    print(safe_object)
