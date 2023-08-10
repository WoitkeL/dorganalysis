 # -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 20:48:14 2023

@author: linus
"""

from iterate_over_DB import iterate_over_database
if __name__ == '__main__':
    safe_object=[]
    database_dict={}
    iterate_over_database(path="C:/Users/linus/python/BioModels-files/",exelname="biomodels_paper")

    print(safe_object)
