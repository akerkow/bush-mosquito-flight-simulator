#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file provides the "Weather" class. After loading weather data from
the clothest weather stations to one of the specified model regions, the
"Weather" class returns the mean temperature for every day for a given period. 
"Model_Szenario.py" provides the start date and the number of years.
"""

#==============================================================================
# Load modules
#==============================================================================
import pandas as pd
import numpy as np

#==============================================================================
# Read the files from the German weather service and extract the average day 
# temperatures ("TMK") and dates
#==============================================================================

fields = ['MESS_DATUM', 'TMK']

def read_csv(filename):
    return pd.read_csv(filename,
                       sep = ';', 
                       delim_whitespace = False,
                       skipinitialspace=True,
                       usecols = fields, 
                       na_values = -999,
                       index_col = "MESS_DATUM")
    
produkt='produkt_klima_tag_'

# Stations for Region 1
St_294 = read_csv(produkt+'19920201_20191231_00294.txt')
St_294new = read_csv(produkt+'20190508_20201107_00294.txt')
St_5715 = read_csv(produkt+'19740101_20191231_05715.txt')
St_5715new = read_csv(produkt+'20190508_20201107_05715.txt')
St_769 = read_csv(produkt+'19720101_20191231_00769.txt')
St_769new = read_csv(produkt+'20190508_20201107_00769.txt')                       
# Stations for Region 2 
St_603 = read_csv(produkt+'19860801_20191231_00603.txt')
St_603new = read_csv(produkt+'20190508_20201107_00603.txt')
St_3490 = read_csv(produkt+'19480101_20191231_03490.txt')
St_3490new = read_csv(produkt+'20190508_20201107_03490.txt')
# Stations for Region 3
St_4177 = read_csv(produkt+'19480101_20191231_04177.txt')
St_4177new = read_csv(produkt+'20190508_20201107_04177.txt')

#==============================================================================
# Data preparation (Merge data from several stations for a region by date, 
# combine with recent data for 2020)
#==============================================================================
def clean(frame):
    if 20191231 in frame.index:
        frame.drop(frame[frame.index < 20200101].index, inplace=True)

clean(St_294new)
clean(St_5715new)
clean(St_769new)
clean(St_603new)
clean(St_3490new)
clean(St_4177new)

#==============================================================================
# Get the mean weather values from all weatherstations within the model region:
#==============================================================================
List_1 = [St_294, St_5715, St_769, St_294new, St_5715new, St_769new]
List_2 = [St_603,St_3490, St_603new, St_3490new]
List_3 = [St_4177, St_4177new]

def merge(StationsData_List):
    concat = pd.concat(StationsData_List)
    byIndex = concat.groupby(concat.index)
    return byIndex.mean().reset_index()

r1= merge(List_1)
r2= merge(List_2)
r3= merge(List_3)

#==============================================================================
# Check for missing data
#==============================================================================
#print("no-data-values for r1:", r1[r1['TMK'].isnull()]=False)
#print("no-data-values for r2:",r2[r2['TMK'].isnull()])
#print("no-data-values for r3:",r3[r3['TMK'].isnull()])

# proof if the year has 365 (or 366) days:
# len(r1[(r1.MESS_DATUM.astype('str')).str.startswith('2018')])

#==============================================================================
# Make the temperature data available by means of a weather class
#==============================================================================
class Weather():
    """ 
    This class returns the mean daily temperature at a given day 
    """
    def __init__(self,region,date):
        """ 
        Weather is initialised with a region (r1,r2,r3) and a starting date 
        """
        self.date=date
        self.region=region
        self.index=self.region[(
            self.region.MESS_DATUM.astype('str')
            ).str.startswith(self.date)].index[0]
        
    def next(self):
        res=self.region.iloc[self.index:self.index+1]
        s=0
        k=0
        for x in res['TMK']:
            if(pd.isna(x)):
                continue
            s+=x
            k+=1
        self.index+=1
        if k!=0:
            return s/k
        return 0