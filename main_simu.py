#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This main simulation combines the spatial component (space_simu.py) with the 
temporal component (time_simu.py). 
The components are called with each simulation step (every day).
"""

#==============================================================================
# Load modules
#==============================================================================

import Model_Szenario as MS
import space_simu 
import time_simu 
import numpy as np
np.set_printoptions(edgeitems=3)
np.core.arrayprint._line_width = 180
import read_temperature as RT
from datetime import datetime

#==============================================================================
# Load the ascii-grid map (mosquito world) with the help of GDAL
#
# The grid map covers the area of Germany, has a resolution of 100m x 100m and 
# shows the habitat suitability for the mosquito species Aedes j. japonicus as 
# values between 0 and 1. The map was developed using a model that integrates 
# climate data and land use data.
#
# There are two ways for reading the data:
# 1. Read the full grid map from the published ascii file and than slice the 
#    array to get the appropriate testregion. The map can be downloaded under 
#    the following link: https://doi.org/10.4228/zalf.dk.90 
#
# 2. Directly call dumped numpy arrays from the testregions (faster!)
#==============================================================================
start=datetime.now()

"""
Possibility 1: (uncommend for usage)
"""
#from osgeo import gdal   
#def read_ascii(filename):
#    g=gdal.Open(filename)
#    return g.ReadAsArray()
#landscape=read_ascii(path+"FuzzyResult_Version3_1981_2010.asc")
    
# northernmost, moderate mosquito density:
#R1=space_simu.region(landscape[2900:3150,2200:2450])

# central:
#R2=space_simu.region(landscape[4800:5050,700:950])   

# southernmost testregion, especially hot:
#R3=space_simu.region(landscape[6710:6960,1640:1890]) 

"""
Possibility 2: (prepare arrays and uncommend for usage)
"""   
R1=space_simu.region(np.load("./R1.npy"))
R2=space_simu.region(np.load("./R2.npy"))
R3=space_simu.region(np.load("./R3.npy"))

#==============================================================================
# Call details for a model szenario
#==============================================================================
testregion = MS.testregion
startdate = MS.startdate
years = MS.years
superagents = MS.superagents

#==============================================================================
# Definition of the main simulation class SIMU (controls the simulation)
# It handles the folling inputs and outputs:
# Inputs:
# Temperature, time steps
# Outputs: 
# LM,SM,EM,IM,SB,EB,IB,RB,DB; Locations of all mosquitoes and the infectious 
# mosquitoes
#==============================================================================
class SIMU():
    def __init__(self,
                 lat,   # latitude of the testregion
                 M      # space_simu.mosquitoes(Region)
                 ):    
        
        self.sir=time_simu.SIR(lat)
        self.space=M 
        self.sir.set_init_conditions(lm=0,# lm=50000
                                     sm=MS.KM*0.15, # =NMmin
                                     im=0,
                                     sb=MS.KB,
                                     ib=0
                                     )
        self.factor=superagents
        self.log_file=open('log_{}_super{}.csv'.format(testregion,
                                                       superagents),'w')
        s=("t mosquitoes LM SM EM IM NB SB EB IB RB DB " + 
           "migrate Emi_N Emi_NE Emi_NW Emi_S Emi_SE Emi_SW Emi_W Emi_E " +
           "EnerKillIM EnerKillNM\n")
        self.log_file.write(s)
        
        """
        Tag day with the highest number of mosquitoes (NM and IM)
        """
        self.highestm=1 # start with an asumption 
        self.highesti=1
        self.IM_killCounter=0
        self.SmEm_killCounter=0
        self.m=None
        self.inf=None
        
    def time_step(self, 
                  t,    # time step counter, starts with 0
                  tx,   # Temperature ('T' in time_simu)
                  yd    # day of the year ('dayOfYear' in time_simu)
                  ):
        """ 
        Run one step using the temperature tx and the yearday yd as inputs;
        save simulation results from the SEIR model (time_simu), corrected by 
        the spatial component (space_simu) in "log.csv"
        """
        self.sir.step(tx,yd) 
        NM = self.sir.SM + self.sir.EM + self.sir.IM
        NB = self.sir.SB + self.sir.EB + self.sir.IB + self.sir.RB
        s=("%d %f %f %f %f %f %f %f %f %f %f %f"+ 
           " %d %d %d %d %d %d %d %d %d %d %d\n") % (t,
                                           NM,
                                           self.sir.LM,
                                           self.sir.SM,
                                           self.sir.EM,
                                           self.sir.IM,
                                           NB,
                                           self.sir.SB,
                                           self.sir.EB,
                                           self.sir.IB,
                                           self.sir.RB,
                                           self.sir.DB,
                                           self.sir.migrate,
                                           self.sir.FirstEmi_N,
                                           self.sir.FirstEmi_NE,
                                           self.sir.FirstEmi_NW,
                                           self.sir.FirstEmi_S,
                                           self.sir.FirstEmi_SE,
                                           self.sir.FirstEmi_SW,
                                           self.sir.FirstEmi_W,
                                           self.sir.FirstEmi_E,
                                           self.IM_killCounter,
                                           self.SmEm_killCounter)
        self.log_file.write(s) 
        
    def space_step(self,t):
        """
        The space_step calls the number of mosquito imagos from the time_simu 
        and turns them into supermosquitoes. 
        
        Autumn-kill: If the space_simu has more mosquito imagos (NonIMs or IMs)
        than time_simu hass calculated for this day, the difference in IM will 
        be killed in the space_simu. 
        
        If the space-Simu has less IMs or non-IMs than time_simu, add the 
        difference to the spacial simulation 
        (Attention: individuals killed due to bad energy level would also 
        reappear. Therefore, their number is saved in space_simu and not 
        added back to the spacial simulation here)
        """
        
        # call mosquitoes from time_simu and make them to supermosquitoes:
        nonIMs = int((self.sir.SM + self.sir.EM)/self.factor)
        IM = int(self.sir.IM/self.factor)
                
        """
        Balance the IMs between the spacial and temporal component
        """
        diffI = IM - len(self.space.infected)
        
        if diffI<0: # if space_simu has more IMs than time_simu...
            self.space.kill_IM(-diffI) # kill diff in space_simu (= autum kill)
            
        if diffI>0: # if space_simu has less IMs than time_simu...
            nr=diffI
            """
            Infectious mosquitoes are not generated from existing susceptible 
            individuals but are newly born. We have implemented this in that 
            simplified way as we do not know exactly where the mosquitoes are 
            infected within the region because we do not consider the birds 
            locally
            """
            self.space.addD(nr,infect=True)
            # created mosquitoes to compensate for the difference to space_simu
        
        """
        Balance the SMs
        """
        EM = int(self.sir.EM/self.factor)
        space_nonIMs = len(self.space.m) - len(self.space.infected)
        
        """
        Calculate the difference between the non-infectious mosquitoes. 
        Do not include the EMs in the check, EMs will not be killed nor added
        on top of those already existent but hiding between the "space_nonIMs", 
        which consist of EMs and SMs
        """
        diff = nonIMs - space_nonIMs - EM # do not add or kill EMs!
        
        if diff<0: # if space_simu has more SMs than time_simu...
            self.space.kill_nonIM(-diff) # kill difference in space_simu
            print("autumn-kill")
        if diff>0: # if space_simu has less SMs than time_simu...
            nr=diff #(energykilled are not added)
            self.space.addD(nr,infect=False) # create difference in space_simu
            
        """
        Space step
        """
        self.space.step(name=testregion)
        
        """
        Check again if space_simu has fewer mosquitoes than time_simu due to
        energy-kill
        """
        
        spaceIM=len(self.space.EnergyKill_IM())     
        IM = int(self.sir.IM/self.factor) 
        
        if spaceIM < IM:
            self.sir.IM = spaceIM * self.factor 
            self.IM_killCounter += (IM-spaceIM) 
            print("Energykill IM: ", IM-spaceIM)
        
        spaceNonIM=(len(self.space.EnergyKill_NonIM())-
                    len(self.space.EnergyKill_IM()))
        nonIMs = int((self.sir.SM + self.sir.EM)/self.factor) 
        
        if spaceNonIM < nonIMs:
            self.sir.SM = spaceNonIM * self.factor 
            self.SmEm_killCounter += (nonIMs-spaceNonIM) 
            print("Energykill SM: ", nonIMs-spaceNonIM)
        
        """
        Update the number of susceptible single mosquitoes by means of the 
        space component. It is only updated when the number of mosquitoes in 
        the flight simulator is smaller than the super-mosquitofactor, because 
        otherwise the result would be smaller than 1.
        """
        if len(self.space.m)>self.factor:          
            self.sir.SM=(len(
                    self.space.m)-len(self.space.infected))*self.factor 
                
            """
            Now, highestm is compared with m to check if it is still the 
            highest number of mosquitoes
            """
            if self.highestm<len(self.space.m): # if highestm is not highest nr          
                self.m = self.space.get_matrixm() #get matrix for visualisation
                self.highestm=len(self.space.m) # update highestm 
                print('____________________________________',self.highestm) 
            if self.highesti<len(self.space.infected):
                self.inf = self.space.get_matrixi()
                self.highesti=len(self.space.infected)
                print('************************************',self.highesti)
            
        """
        Save all mosquitoes as well as the infectious within the region when 
        the simulation runs only for a short time period (one year)
        """
        if years <=1:
            np.save('DaylyM_{}_{}.npy'.format(testregion,t),
                    self.space.get_matrixm().astype(int))
            np.save('DaylyInf_{}_{}.npy'.format(testregion,t),
                    self.space.get_matrixi().astype(int))
        
        """
        Add number of mosquitoes that have tried to emigrate (also in specific
        directions) to the time_simu in order to save the information daily)
        """
        self.sir.migrate = self.space.migrate
        self.sir.FirstEmi_N = self.space.FirstEmi_N
        self.sir.FirstEmi_NE = self.space.FirstEmi_NE
        self.sir.FirstEmi_NW = self.space.FirstEmi_NW
        self.sir.FirstEmi_S = self.space.FirstEmi_S
        self.sir.FirstEmi_SE = self.space.FirstEmi_SE
        self.sir.FirstEmi_SW = self.space.FirstEmi_SW
        self.sir.FirstEmi_W = self.space.FirstEmi_W
        self.sir.FirstEmi_E = self.space.FirstEmi_E
                
    def save_matrix_mostInfections(self):
        """
        Save matrix at the moment with the highest number of infections
        """
        np.save('matrixm_HI_{}_super{}.npy'.format(testregion,
               superagents),self.m)
        np.save('matrixi_HI_{}_super{}.npy'.format(testregion,
                superagents),self.inf)

#==============================================================================
# Execute the simulation
#==============================================================================        
if testregion=="R1":
    M=space_simu.mosquitoes(R1)
    rt=RT.Weather(RT.r1,startdate)
    print("Region 1")
elif testregion=="R2":
    M=space_simu.mosquitoes(R2)
    rt=RT.Weather(RT.r2,startdate)
    print("Region 2")
elif testregion=="R3":
    M=space_simu.mosquitoes(R3) 
    rt=RT.Weather(RT.r3,startdate)
    print("Region 3")
else:
    print("Set 'testregion' correctly as R1, R2 or R3")
    
simu=SIMU(50.6, # latitude
          M) 
tx=np.arange(years*365)

""" 
Iterate over one year (the step function contains the "dayOfYear") for every 
day + year of the whole simulation period (rt) and add infectious birds at 
specific days:
    
60 (Begin of March): Birds start to build their nests
210 (End of Juli): End of main chick hatching period 
280 (Begin of October): Autumn migration   
"""
T0=rt.next() # no smoothing of the first data value
for t in tx:
    if t%365==60: 
        simu.sir.IB+=2
    if t%365==210:
        simu.sir.IB+=2
    if t%365==280:
        simu.sir.IB+=4
    T=(rt.next())*0.1+T0*0.9 # weighted average filter
    T0=T
    simu.time_step(t,T,t%365)
    simu.space_step(t) # Coupling with the space model

simu.save_matrix_mostInfections()
simu.log_file.close()

print("Time: ",datetime.now()-start)
