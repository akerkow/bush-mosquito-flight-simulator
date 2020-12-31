#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 
This simulation is derived from the SEIR model from Laperriere et. al. 2010 
(doi:10.1098/rspb.2003.2608). It implements the development of host birds, 
vector mosquitoes and virus transmissions dependent on local temperatures and 
daylenghts. 

We adapted the model to a probable European sentinel bird species (Pica pica) 
and another vector species (Aedes japonicus japonicus). Additionally, we 
adopted the model in a way that it takes into account mosquito blood meals on 
organisms other than the host birds. 
"""

import numpy as np
from scipy.stats import gamma
import Model_Szenario as MS


class SIR:
    def __init__(self,lat):
        """
        Regional parameter
        """
        self.lat = lat # latitdude of the region
        
        """ 
        Mosquito init parameter 
        """
        self.species = MS.species
        self.KM = MS.KM   # carrying capacity of mosquito larvae
        self.NMmin = self.KM*0.15 # minimum number of mosquito imagos
        
        """ 
        Bird init parameter
        """
        self.KB = MS.KB # carrying capacity of magpies (equal in all regions)        
        self.mB = 0.001404   # Eurasian magpie (WNV)
        
        # Magpies interaction with West-Nile-virus:         
        self.alphaB = 0.29   # removal rate of magpies 
        self.gammaB = 0.333  # incubation rate among magpies
        self.nuB = 0.43      # portion of dead magpies due to WNV-2 infection
        
        """ 
        Probability of Transmission from IM to SB 
        """
        self.pM = 1.0 # every mosquito-bird contact leads to an infection
        
        """ 
        Probability of Transmission from IB to SM
        """
        if self.species == "japonicus":
            self.pB = 0.17 
            self.birdBites= MS.JapBirdBites
            print("Aedes japonicus")
            
        elif self.species == "culex":
            self.pB= 0.06
            self.birdBites= MS.CuBirdBites
            print("Culex pipiens ssp.")
                            
        else:
            print('Set "species" (init condition) as "culex" or "japonicus"')
                
    """ 
    Mosquito parameter dependent on temperature or daylength
    """    
    def biting_rate(self,T): # T: temperature in °C
        return (0.344/(1+1.231*np.exp(-0.184*(T-20))))
    
    def birdBites(self):
        if self.species=="japonicus":
            return MS.JapBirdBites
        elif self.species =="culex":
            return MS.CuBirdBites
        else:
            print('Set "species" (init condition) as "culex" or "japonicus"')
    
    def bL(self,T): # Birth rate of mosquito larvae
        return 2.325 * self.biting_rate(T)
    
    def bM(self,T): # Birth rate of mosquito imagos
        return self.bL(T)*0.1 # from Rubel
    
    def mL(self,T): # Mortality of the larvae
        return 0.0025*T**2-0.094*T+1.0257 
    
    def mM(self,T): # Mortality rate of mosquito imagos
        return 0.1*self.mL(T)
    
    def betaM(self,T): # Transmission rate
        return self.biting_rate(T)*self.pM*self.birdBites*0.75
    
    def daylength(self,dayOfYear):
        """
        Function computes the length of the day as time between sunrise and 
        sunset, given the day of the year and latitude of the location 
        according to the Brock model
        -----------------------------------------------------------------------
        Inputs:
        -----------------------------------------------------------------------
        dayOfYear : int
        The day of the year. 1 corresponds to 1st of January and 365 to 31st of
        December (on a non-leap year).
        -----------------------------------------------------------------------
        lat : float
        Latitude of the location in degrees. Positive values for the north and 
        negative for the south.
        -----------------------------------------------------------------------
        d : float
        Daylength in hours.
        -----------------------------------------------------------------------
        """
        latInRad = np.deg2rad(self.lat)
        declinationOfEarth = 23.45*np.sin(np.deg2rad(360.0*(283.0+dayOfYear
                                                            )/365.0))
        if -np.tan(latInRad)*np.tan(np.deg2rad(declinationOfEarth)) <= -1.0:
            return 24.0
        elif -np.tan(latInRad)*np.tan(np.deg2rad(declinationOfEarth)) >= 1.0:
            return 0.0
        else:
            hourAngle = np.rad2deg(np.arccos(
                    -np.tan(latInRad)*np.tan(np.deg2rad(declinationOfEarth))))
        return 2.0*hourAngle/15.0
    
    def deltaM(self,dayOfYear): # Fraction of active, not diapausing mosquitoes
        if self.species=="japonicus":
            return 1.0-1.0/(1.0+40000*np.exp(1.559*(
                    self.daylength(dayOfYear-30)-18.177)))
        elif self.species =="culex":
            return 1.0-1.0/(1.0+1775.7*np.exp(1.559*(
                    self.daylength(dayOfYear)-18.177)))
        else:
            print('Set "species" (init condition) as "culex" or "japonicus"')
    
    
    def gammaM(self,T): # rate infected-infectious
        if(T<=15): 
            return 0.0
        return 0.0093*T-0.1352
    
    """ 
    Bird parameter dependent on temperature, daylength an day of the year
    """   
    def bB(self,dayOfYear): # Bird birth rate
        x = dayOfYear-120 # transformed Julian calender day (120=loc)
        if(x<=0):
            return 0.0
        return 0.614 * gamma.pdf(x, a=4.43, scale=7.67,loc=0)
        
    def betaB(self,T): # virus transmission rate to susceptible mosquito
        return self.biting_rate(T)*self.pB*self.birdBites*0.75
        # 0.75 = fraction of host bird bites from all bird bites
    
    
    """ 
    Simulation part 
    """
    def set_init_conditions(self,
                            lm=1,
                            sm=MS.KM*0.15, #=NMmin
                            im=0,
                            sb=MS.KB,
                            ib=0,
                            migrate = 0,
                            FirstEmi_N = 0,
                            FirstEmi_NE = 0,
                            FirstEmi_NW = 0,
                            FirstEmi_S = 0,
                            FirstEmi_SE = 0,
                            FirstEmi_SW = 0,
                            FirstEmi_W = 0,
                            FirstEmi_E = 0):
        
        """ 
        Init the simulation by setting the starting parameters 
        """
        self.SM = sm     # susceptible mosquitoes
        self.LM = lm     # mosquito larvae
        self.EM = 0.0    # exposed mosquitoes
        self.IM = im     # infectious mosquitoes
        
        self.SB = sb     # susceptible birds
        self.EB = 0.0    # exposed birds
        self.IB = ib     # infectious birds
        self.RB = 0.0    # recovered and immune birds
        self.DB = 0.0    # birds succumbed to infection
        
        """"
        Derivations from the space_simu
        """
        self.migrate = 0 # counts mosquitoes that tried to leave R
        self.FirstEmi_N = 0 #counts attempts to escape to the north...
        self.FirstEmi_NE = 0
        self.FirstEmi_NW = 0
        self.FirstEmi_S = 0
        self.FirstEmi_SE = 0
        self.FirstEmi_SW = 0
        self.FirstEmi_W = 0
        self.FirstEmi_E = 0
        
    def lambdaMB(self,T, dayOfYear): # Virus transfer from mosquitoes to birds 
        phi = self.KM/self.KB
        return self.deltaM(dayOfYear)*self.betaM(T)*phi *self.IM/self.KM
        
    def lambdaBM(self,T, dayOfYear): # Virus transfer from birds to mosquitoes
        return self.deltaM(dayOfYear)*self.betaB(T)*self.IB/self.KB 
            
    def step(self,T,dayOfYear):
        """ 
        One time step with a simple euler method dh=1,
        T= temp in degree, dayOfYear
        """
        SM = self.SM # save the old state for reconstraction
        EM = self.EM
        IM = self.IM
        LM = self.LM

        SB = self.SB
        EB = self.EB
        IB = self.IB
        RB = self.RB

        dayOfYear+=1
        
        """
        Bird population loop
        """ 
        bB = self.bB(dayOfYear)
        NB = SB+EB+IB+RB
        self.SB+= ((bB-(bB-self.mB)*NB/self.KB)*NB-self.lambdaMB(T,dayOfYear)*
                   SB-self.mB*SB)
        self.EB+= self.lambdaMB(T,dayOfYear)*SB-self.gammaB*EB-self.mB*EB
        self.IB+= self.gammaB*EB-self.alphaB*IB-self.mB*IB 
        self.RB+= (1-self.nuB)*self.alphaB*IB-self.mB*RB
        self.DB+= self.nuB*self.alphaB*IB
        
        """
        Mosquito population loop
        """
        NM = SM+EM+IM 
        self.LM+= (self.bL(T)*self.deltaM(dayOfYear)*NM-self.mL(T)*LM)*(
                1-LM/self.KM)-self.bM(T)*self.LM            
        self.SM+= self.bM(T)*LM-self.mM(T)*SM-self.lambdaBM(T,dayOfYear)*SM
        self.EM+= self.lambdaBM(T,dayOfYear)*SM-self.gammaM(T)*EM-self.mM(T)*EM
        self.IM+= self.gammaM(T)*EM-self.mM(T)*IM
        
        """
        Check if valid
        """
        if self.SM<0 :
            self.SM=100
        NM = self.SM+self.EM+self.IM
        if NM<self.NMmin : 
            self.SM = self.NMmin
        if(self.EM<0):
            self.EM = 0.1
        if(self.IM<0):
            self.IM = 0.1
        if(self.LM<0):
            self.LM = 0.1
        if(self.SB<0):
            self.SB = 0.1