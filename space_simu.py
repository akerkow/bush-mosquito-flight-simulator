#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
space_simu is the spatial extension for the SEIR model and is called by the 
main simulation. The spatial simulation models the movements of the 
mosquitoes and is therefore also called "flight simulator". It is based on the 
agent-based modeling technique. 
    
Space_simu imports the number of non-diapausing mosquito imagos with their 
state of infection from the SEIR model. 

Returning to the main simulation, it provides information on how many 
mosquitoes died due to staying in unsuitable habitats. It also returns the 
distribution of the mosquitoes within the defined study area.

"Model_Szenario.py" controlls if the model runs in an open or closed mode. If 
the model is in closed mode, the mosquitoes cannot leave the study region.
"""

#==============================================================================
# Load modules
#==============================================================================
import numpy as np
import random
from math import sqrt
import matplotlib.pylab as plt
import pickle
import Model_Szenario as MS

#==============================================================================
# The class "region" provides the region object. One of 3 testregions can be
# selected in "main_simu.py" and is referred to as "self.R". 
# The testregion is a numpy array. Different functions allow e.g. to easily 
# get the size of the testregion or a specific position within the array.
#==============================================================================
class region:

    def __init__(self, nparray):
        self.R=nparray
        
    def get_shape(self):
        """ 
        Get the shape (number of columns and raws) of the region 
        """
        return self.R.shape[0], self.R.shape[1]
  
    def get(self,i,j):
        """ 
        Get the position within the array (row, column) 
        """
        return self.R[i,j]
    
    def show(self):
        """ 
        Show the habitat characteristics within the region 
        """
        plt.figure(figsize=(12,12))
        plt.imshow(self.R)
        plt.colorbar()
        plt.show()
        
    def save_region(self,filename):
        """ 
        Save the region as numpy array 
        """
        np.save(filename,self.R)
                
    def load_region(self,filename):
        """ 
        Reload region from a dumped numpy array 
        """
        try:
            self.R=np.load(filename)
        except OSError:
            print('cannot open:', filename)
            return False
        return True
        
#==============================================================================
# The class "mosquito" contains the mosquito objects (every individual mosquito 
# or supermosquito) that is not diapausing). They have individual proberties 
# which are e.g. the location within the study region, the study region they 
# fly on, the infection status and the amount of energy they have.
#        
# Energy: Depending on the habitat quality of its location, the mosquito has a 
# motivation to leave that patch within a defined, but rondomly controlled
# range within a closer or wider moore neighbourhood. The mosquito can gain 
# energy from beeing located in a good patch and looses energy constantly as 
# well as during the flight. The energy loss during the flight depends on the 
# flight distance.
#==============================================================================
class mosquito:
    def __init__(self, 
                 i1, j1, 
                 id, 
                 R, 
                 energy=100.0, 
                 infect=False,
                 dirFirstEmi=None, 
                 pulledBack=False):
        
        self.i=i1          # row of the location in R
        self.j=j1          # column of the location in R
        self.id=id         # name of the mosquito
        self.infect=infect # tag if the mosquito is infected 
        self.energy=energy # energy level at the begin of the simulation
        self.R=R           # use the defined region
        self.dirFirstEmi=dirFirstEmi # direction of first emigration attempt
        self.pulledBack=pulledBack # tag if mosquito tried to emigrate before
        
    def jump(self,motivation,Region):
        """ 
        Mosquito random flight 
        """
        def jumpConsequences(i,j):
            height,width=self.R.get_shape() # call extention of study region
            
            """
            Check in which direction the mosquito wanted to leave the area and
            don't execute the jump. The mosquito stays at its place, but that
            does not mean that it gets killed yet.
            """
            
            # check if it did not try to leave R before:
            if self.dirFirstEmi is None: # no attempt yet!
                I = self.i+i 
                J = self.j+j
                #array indices: 0...249, width and height = 250
                
                # north
                if (I<0 and J>=0 and J<width):
                    self.dirFirstEmi="N"
                    return False # no execution of jumpConsequences(i,j)
                # northeast
                elif (I<0 and J>(width-1)):
                    self.dirFirstEmi="NE"
                    return False
                # northwest
                elif (I<0 and J<0):
                    self.dirFirstEmi="NW"
                    return False
                # south
                elif (I>(height-1) and J>=0 and J<width):
                    self.dirFirstEmi="S"
                    return False
                # southeast
                elif (I>(height-1) and J>(width-1)):
                    self.dirFirstEmi="SE"
                    return False
                # southwest
                elif (I>(height-1) and J<0):
                    self.dirFirstEmi="SW"
                    return False
                # west
                elif (I>=0 and I<height and J<0):
                    self.dirFirstEmi="W"
                    return False
                # east
                elif (I>=0 and I<height and J>(width-1)):
                    self.dirFirstEmi="E"
                    return False
                
            else:
                self.pulledBack = True # tried to emigrate before!
                # also no execution when it tries to emigrate again...
                if((self.i+i)<0 or (self.i+i)>=height or (self.j+j)<0 or 
                   (self.j+j)>=width):
                    return False
            
            self.i += i # my new row position
            self.j += j # my new column position
            
            """
            Energy lost during the flight (using Hypotenuse for calculation
            of the flight distance and a scaling factor of 50)
            """
            SteadyLoss = 4
            self.energy -= (8*(sqrt(np.abs(i)**2 + np.abs(j)**2)))-SteadyLoss
            
            """
            Energy gain from R
            """
            self.energy += 50*self.R.get(self.i,self.j)
            if(self.energy > 100.0):
                self.energy = 100.0
                return True 
        
        if motivation == 0:
            """
            Low flight motivation
            """
            FlightDist=round(np.random.normal(loc=0, scale=0.6))
            if FlightDist !=0:
                """ Begin flight in vertical direction (N or S) """
                iDist = random.randint(0,abs(FlightDist))
                if FlightDist > 0:
                    i = iDist 
                elif FlightDist < 0:
                    i = -iDist
                """ Continue the Flight in horizontal direction (W or E) """
                jDist = abs(FlightDist)-iDist
                j = [-1,1][random.randrange(2)]*jDist
            else:
                i=0
                j=0
            jumpConsequences(i,j)

        elif motivation == 1: 
            """
            Average flight motivation
            """
            FlightDist = round(np.random.normal(loc=0, scale=2.4))
            if FlightDist !=0:
                """ Begin flight in vertical direction (N or S) """
                iDist = random.randint(0,abs(FlightDist))
                if FlightDist > 0:
                    i = iDist 
                elif FlightDist < 0:
                    i = -iDist
                """ Continue the Flight in horizontal direction (W or E) """
                jDist = abs(FlightDist)-iDist
                j = [-1,1][random.randrange(2)]*jDist # random choice W or E
            else:
                i=0
                j=0
            jumpConsequences(i,j)
            
        elif motivation == 2:
            """
            Extremly high flight motivation
            -> mosquito flights are corrected by the local mean wind
               conditions and thus slightly different for every region
            """
            if Region == "R1":
                i=np.random.randint(-3,4) # North, South
                j=np.random.randint(-4,3) # West, East
                jumpConsequences(i,j)
            elif Region == "R2":
                i=np.random.randint(-3,4) # North, South
                j=np.random.randint(-4,3) # West, East
                jumpConsequences(i,j)
            elif Region == "R3":
                i=np.random.randint(-3,4) # North, South
                j=np.random.randint(-4,3) # West, East
                jumpConsequences(i,j)
            else:
                print("Region not defined")
        else:
            print("Jump motivation unvalid")
            
#==============================================================================
# The class "mosquitoes" controls all mosquitoes of the simulation.
#==============================================================================           
class mosquitoes:
    
    def __init__(self,R):
        self.m=[] # list for all mosquitoes (infected and non-infected)
        self.infected=[] # list for infectious mosquitoes
        self.steps=0 # count simulation steps (days)
        self.R=R     # use the defined region
        self.hist2D=self.calc_hist2D()  # array(50,50) with habitat quality 
                                        # sums for every 25ha sized patch
        self.matrix=np.zeros(self.R.get_shape())#store mosquitoes in this array
        self.energy=np.zeros(self.R.get_shape()) # store energy in this array
        self.distribution = None # define the "distribution" object
        self.migrate = 0 # counts mosquitoes that tried to leave R
        self.FirstEmi_N = 0 # count a try to emigrate to the north...
        self.FirstEmi_NE = 0
        self.FirstEmi_NW = 0
        self.FirstEmi_S = 0
        self.FirstEmi_SE = 0
        self.FirstEmi_SW = 0
        self.FirstEmi_W = 0
        self.FirstEmi_E = 0
           
    def calc_hist2D(self):
        """ 
        Calculates the spacial hist2D. This is the region of the simulation
        that provides the habitats for the mosquitoes, but in a lower spatial 
        resolution (the occurrence probabilities are summarised in 5x5 unit 
        subarrays which is 500x500m2 or 25ha)
        """
        hight, wide = self.R.get_shape()
        hight //= 5 #250/5= 50 grid units= 5km (floor devision returns integer)
        wide //= 5 
        hist2D=np.zeros((hight, wide)) # empty copy from 1/5 testregion
        for i in range(hight):  # fill the 2D historgram
            for j in range(wide):
                hist2D[i,j]=np.sum(self.R.R[5*i:5*i+5,5*j:5*j+5])
                #print(np.sum(hist2D), np.shape(hist2D))
        # array with 25ha sized patches that sum up the habitat qualities:  
        hist2D = hist2D/hist2D.sum() #get probability values between 0 and 1 
        return hist2D
    
          
    def set_one_mosquito(self,
                         sumOfhist2D,   # sum of the broad hist2D habitat array
                         infect=False   # status will get overset by addD()
                         ):
        """ 
        Implementation of the mosquito birth process: 
        """
        # thresh is an individual random habitat quality threshold for the 
        # mosquito to be placed in a 25ha patch
        
        #new:
        thresh=np.random.choice(self.hist2D.ravel(),1,
                                p=self.hist2D.ravel(), # weighted choice
                                replace=False)
        possiblePositions=list(zip(*np.where(self.hist2D>=thresh)))
        FinalPosition=np.random.permutation(possiblePositions)[0]
        self.distribution[FinalPosition[0], FinalPosition[1]]+=1 
        
        """
        Define random position within a 25ha sized patch
        """
        k1=np.random.randint(0,5) #includs 0,1,2,3 or 4
        k2=np.random.randint(0,5)
                    
        """
        Set the mosquito properties:
        i,j = position in the hist2d array 
        i1, j1 = position in the original array
        k1, k2 = range of hist2D within original array
                    """
        m=mosquito(i1 = FinalPosition[0]*5+k1,
                   j1 = FinalPosition[1]*5+k2,
                   id = 0,
                   R = self.R, 
                   infect=infect, 
                   dirFirstEmi=None, 
                   pulledBack=False) 
                    
        """
        Add mosquito to the list with all mosquitoes
        """
        self.m.append(m) # get born!
        
        """
        Eventually also add it to the list of infected mosquitoes
        """
        if(infect==True): 
            self.infected.append(m)
        return True
                
    def addD(self,
             nr, # number of mosquitoes to add
             infect # infection status
             ): 
        """
        When the main simulation detects that the space_simu has eather less 
        IMs or SMs than time simu, than this function gets called. 
        "addD()" calculates a distribution for the missing mosquitoes of the 
        respective infection status and adds them to the 2D-histogram. 
        The infection status is controlled in main_simu.
        
        The 2D-histogram is a matrix with (shape[0]/5,shape[1]/5) cells, which 
        is (50:50) and every cell covers 25ha. 
        The result is the "distribution" matrix. This matrix has the same size 
        like hist2D and gets fiiled with a part of the number of mosquitoes to 
        be reborn (nr). The reborn-process is handled by 
        "self.set_one_mosquito"
        """
        
        hight, wide = self.R.get_shape() 
        hight //= 5 # get size of the hist2D-matrix
        wide //= 5
        
        self.distribution=np.zeros((hight, wide))  # result matrix (empty)
        sumOfhist2D=np.sum(self.hist2D) # self.hist2D=self.calc_hist2D()
        
        for k in range(nr):
            """
            add to m[] and infected[] if infect=True, otherwise just to m[]
            """
            self.set_one_mosquito(sumOfhist2D,infect)
        return True 
    
    def step(self,name):
        self.steps+=1 # step one day for all stored mosquitoes
        if MS.closed==False:
            ni = 0 # counter for infectious mosquitoes with unvalid jump_tag
            nm = 0 # counter for susceptible mosquitoes with unvalid jump_tag

        self.matrix[:,:]=0.0  # reset the mosquito matrix (original,250x250)
        allm=self.m+self.infected
        for m in allm: # fill matrix with the mosquitoes
            self.matrix[m.i,m.j]+=1
            
            """
            Get jump motivation for the mosquitoes depending on the habitat 
            quality at their location in the main matrix.
            Also check which region is handeled in the moment, this information
            is important in order to control the wind effect in the open
            landscapes with bad habitat quality
            """
            if (self.R.get(m.i,m.j)>=0.8): # test for low motivation
                jump_tag=m.jump(motivation=0, Region=name)
            elif (self.R.get(m.i,m.j)>= 0.3): # test for average motivation
                jump_tag=m.jump(motivation=1, Region=name)
            else:
                jump_tag=m.jump(motivation=2,Region=name) # low motivation
            
            if MS.closed==False:
                if(jump_tag is None):
                    if(m.infect in self.infected):
                        ni+=1 # count infectious mosquitoes
                        self.infected.remove(m)
                        self.m.remove(m)
                    elif(m.infect==False):
                        nm+=1 # count susceptible mosquitoes
                        self.m.remove(m)
            
            """
            Count mosquitoes that tried to leave the region for the first time,
            all of them, and separately for every possible direction of the 
            Moore neighbourhood
            """
            if(jump_tag is None and m.pulledBack==False):
                self.migrate += 1 # all
                if m.dirFirstEmi == "N": 
                    self.FirstEmi_N+=1
                elif m.dirFirstEmi == "NE":
                    self.FirstEmi_NE+=1
                elif m.dirFirstEmi == "NW":
                    self.FirstEmi_NW+=1
                elif m.dirFirstEmi == "S":
                    self.FirstEmi_S+=1
                elif m.dirFirstEmi == "SE":
                    self.FirstEmi_SE+=1
                elif m.dirFirstEmi == "SW":
                    self.FirstEmi_SW+=1
                elif m.dirFirstEmi == "W":
                    self.FirstEmi_W+=1
                elif m.dirFirstEmi == "E":
                    self.FirstEmi_E+=1
        
            """
            Out-of-energy-kill:
            Mosquitoes with bad energy level get removed from their lists
            """
            if(m.energy<1):
                if(m.infect==False):
                    self.m.remove(m)
                elif(m.infect==True and m in self.infected):
                    self.infected.remove(m)
                    self.m.remove(m)
        
        if MS.closed==False:
            """
            Add mosquitoes with "jump_tag is None" to the grid again. These 
            mosquitoes tried to leave the region.
            """                
            self.addD(nr=nm,infect=False) # add susceptable mosquitoes
            self.addD(nr=ni,infect=True) # add infectious mosquitoes
    
    def kill_IM(self,n):
        """
        IM random kill for autumn:
        Removes infectious mosquitoes at random locations when space_simu has 
        more than time_simu has calculated for this day. This situation happens 
        after the mosquito population reached its maximum size in summer.
        """
        IMs = [m for m in self.m if m.infect == True]
        ms=np.random.permutation(IMs)
        for k in range(n):
            self.m.remove(ms[k])
            self.infected.remove(ms[k])
        return True
    
    def kill_nonIM(self,n):
        """ 
        SM + EM random kill for autumn:
        """
        nonIMs = [m for m in self.m if m.infect == False]
        ms=np.random.permutation(nonIMs)
        for k in range(n):
            self.m.remove(ms[k])
        return True
    
    '''
    Out-of-energy-kill:
    Mosquitoes with bad energy level get first tagged and counted in order to
    keep them removed even after the harmonisation of the spatial and temporal
    component. Than, they get removed from their lists.
    '''
    def IM_noEnergy(self):
        killCount = [m for m in self.m if (m.energy<1 and m.infect==True and 
                                           m in self.infected)]
        nr = len(killCount)
        return nr # counter for killing them also in main/time_simu
    
    def NonIM_noEnergy(self):
        killCount = [m for m in self.m if (m.energy<1 and m.infect==False)]
        nr = len(killCount)
        return nr
                
    def EnergyKill_IM(self):
        kill_from_infected = [m for m in self.infected if (m.energy<1)]  
        
        for k in kill_from_infected:
            print("energy-kill: {} infected mosquitoes".format(k))
            del k
        return self.infected 
    
    def EnergyKill_NonIM(self):
        kill_from_m = [m for m in self.m if m.energy<1]
    
        for k in kill_from_m:
            del k
        return self.m
    
    """
    Visualising the individual mosquitoes at their respective locations within 
    the matrix and functions for saving and loading the matrix:
    """  
    def show(self):
        """ 
        Show the non-infectious mosquitoes
        """
        self.matrix[:,:]= 0.0  # reset the mosquito matrix
        for m in self.m:      # fill it with the mosquitoes
            self.matrix[m.i,m.j]+=1
        plt.imshow(self.matrix)
        plt.colorbar()
        plt.show()
        
    def show_energy(self):
        """ 
        Show the energy level of the mosquitoes 
        """
        self.energy[:,:]= 0.0
        for m in self.m:
            self.energy[m.i,m.j]+=m.energy
        plt.imshow(self.energy)
        plt.colorbar()
        plt.show()
        
    def get_number(self):
        """ 
        Retruns the number of mosquitoes wich are alive and in R 
        """
        return len(self.m)
    
    def save_mosquitoes(self,filename):
        """ 
        Save the list that contains all mosquitoes 
        """
        f=open(filename,'wb')
        pickle.dump(self.m,f)
        f.close()
        
    def load_mosquitoes(self,filename):
        """ 
        Reloads the mosquitoes from disk 
        """
        try:
            f=open(filename,'rb')
        except OSError:
            print('can not open:', filename)
            return False
        self.m=pickle.load(f)
        f.close()
        return True
    
    def get_matrixm(self): 
        """ 
        Matrix with non-infectious mosquitoes
        """
        matrix=np.zeros(self.R.get_shape())  # reset the mosquito matrix
        for m in self.m:      # fill it with the non-infectious mosquitoes
            matrix[m.i,m.j]+=1.0
        return matrix
    
    def get_matrixi(self): 
        """ 
        Matrix with mosquitoes tagged as infectious
        """
        matrix=np.zeros(self.R.get_shape())  # reset the mosquito matrix
        for m in self.infected:      # fill it with the infectious mosquitoes
            matrix[m.i,m.j]+=1.0 
        return matrix
        

        
        
    
