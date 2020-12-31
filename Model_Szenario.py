#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Model_Szenario.py" defines parameters and scenarios for certain model 
applications. The parameters are called from:
    - time_simu.py
    - visualize_time_simu.py
    - space_simu.py
    - main_simu.py
    - visualize_log.py
    - animation_NM.py
    - animation_IM.py.  
"""

species = "culex" # declare "japonicus" or "culex"
JapBirdBites = 0.165 # Japonicus: fraction of bird bites from all bites
CuBirdBites = 0.96 # Culex: fraction of bird bites from all bites
testregion ="R3" 
closed = True # if True, mosquitoes cannot emigrate 
startdate = "20180101"
years = 1
superagents = 1000
KB = 8226
highestKM = 19900000

R1factor = 0.35
R2factor = 0.66
R3factor = 0.79

if testregion == "R1":
    KM = int(highestKM*R1factor)
elif testregion == "R2":
    KM = int(highestKM*R2factor)
elif testregion == "R3":
    KM = int(highestKM*R3factor)
else:
    print('set testregion in Model_Szenario.py to "R1","R2" or "R3"')
print("KM = ", KM)