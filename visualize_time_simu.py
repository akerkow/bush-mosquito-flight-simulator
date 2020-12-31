# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 23:15:47 2020

@author: antje
"""
import time_simu
import numpy as np
import pylab as plt
import read_temperature as rt 
import Model_Szenario as MS

#==============================================================================
# Provide a weather plot for the period specified in "Model_Szenario.py"
#==============================================================================        
years = MS.years
date = MS.startdate
    
def makeList(w):
    TX=[]
    T0=0
    for i in range (365*years):
        rawData=(w.next())
        T1 = rawData *0.1 + T0*0.9
        T0 = T1
        TX.append(T1)
    return(TX)        

TX_r1 = makeList(rt.Weather(region = rt.r1, date = date))[0:365*years+1]
TX_r2 = makeList(rt.Weather(region = rt.r2, date = date))[0:365*years+1]
TX_r3 = makeList(rt.Weather(region = rt.r3, date = date))[0:365*years+1]


fig, (ax1, ax2, ax3) = plt.subplots(3,sharex=True,
                                    sharey=True,
                                    figsize=(12,6),
                                    dpi=300) 

ax1.plot(TX_r1, label="Region 1",color="black")
ax1.minorticks_on()
ax1.grid(b=True, which='minor', color='0.9', linestyle='--')
ax1.grid(b=True, which='major', color='0.7', linestyle='-')
ax1.tick_params(axis='y', which='major', labelsize=12)
ax1.legend()

ax2.plot(TX_r2, label="Region 2",color="black")
ax2.minorticks_on()
ax2.grid(b=True, which='minor', color='0.9', linestyle='--')
ax2.grid(b=True, which='major', color='0.7', linestyle='-')
ax2.tick_params(axis='y', which='major', labelsize=12)
ax2.legend()

ax3.plot(TX_r3, label="Region 3",color="black")
ax3.minorticks_on()
ax3.grid(b=True, which='minor', color='0.9', linestyle='--')
ax3.grid(b=True, which='major', color='0.7', linestyle='-')
ax3.tick_params(axis='y', which='major', labelsize=12)
ax3.legend()

Pos=(np.arange(0,365*(years+1),365))
dates=['01/01/2018', '01/01/2019', '01/01/2020', '01/01/2021']
plt.xticks(Pos, dates, size="large")
plt.tight_layout()
plt.show()

#==============================================================================
# Show plots for the temperature and season dependant parameters
#==============================================================================
def show_tests():
    lat=50.6
    sir=time_simu.SIR(lat)
    
    months=["1\n(Jan)","60\n(Mar)","121\n(May)","182\n(July)","244\n(Sept)",
            "305\n(Nov)","366\n(Jan)"]
    dates=[1,60,121,182,244,305,366]
    
    """ 
    biting_rate, bL und bM sowie mL and mM 
    """
    
    T=np.linspace(-5.0,45.0,1000)    
    biting_rate=sir.biting_rate(T)
    bL=sir.bL(T)
    bM=sir.bM(T)    
    mL=[]
    mM=[]
    for TX in T:
        mL.append(sir.mL(TX))
        mM.append(sir.mM(TX))
        
    plt.figure(figsize=(12,4),dpi=300)
    plt.subplot(121)
    plt.plot(T,biting_rate,linestyle='dotted',color="black",
             label="biting rate")
    plt.plot(T,bL,linestyle='dashed',color="black",
             label="birth rate of larvae")
    plt.plot(T,bM,linestyle='solid',color="black",
             label="birth rate of imagos")
    plt.grid(linestyle="dotted",c="0.8")
    plt.legend()
    plt.xlabel('Temperature (°C)')
    plt.ylabel('Daily rate')
    plt.title('a)' , loc='left',fontweight='bold')
    
    plt.subplot(122)
    plt.plot(T,mM, linestyle='solid', color="black", label="larvae")
    plt.plot(T,mL, linestyle='dashed', color="black", label="imagos")
    plt.grid(linestyle="dotted",c="0.8")
    plt.legend()
    plt.xlabel('Temperature (°C)')
    plt.ylabel('Daily mortality rate')
    plt.title('b)' , loc='left',fontweight='bold')
    plt.tight_layout(w_pad=3)
    #plt.savefig('PopulationControl_new.pdf',bbox_inches='tight')
    plt.show()
    
    """ 
    gammaM and EIP
    """
    gammaM=[]
    EIP=[]
    for TX in T:
        gammaM.append(sir.gammaM(TX))
    for TX in T:
        if(sir.gammaM(TX)<=0):
            EIP.append(170)
        else:
            EIP.append(1/sir.gammaM(TX))
        
    plt.figure(figsize=(12,4),dpi=300)
    
    plt.subplot(122)
    plt.plot(T,gammaM, color="black")
    plt.grid(linestyle="dotted",c="0.8")
    plt.xlim([12,35])
    plt.ylim([0.0,0.2])
    plt.xlabel('Temperature (°C)')
    plt.ylabel('$\gamma_{M}$',labelpad=10)
    plt.title('b)' , loc='left',fontweight='bold')
    
    plt.subplot(121)
    plt.plot(T,EIP, color="black")
    plt.grid(linestyle="dotted",c="0.8")
    plt.ylim([0,100])
    plt.xlim([12,35])
    plt.xlabel('Temperature (°C)')
    plt.ylabel('EIP in days',labelpad=10)
    plt.title('a)' , loc='left',fontweight='bold') 
    plt.tight_layout(w_pad=3)
    #plt.savefig('gammaM_EIP.pdf',bbox_inches='tight')
    plt.show()
    
    """ 
    deltaM 
    """
    t=np.linspace(1,365,365)
    deltaM=[]
    for tx in t:
        deltaM.append(sir.deltaM(tx))
    plt.figure(dpi=300)
    plt.plot(t,deltaM, color="black")
    plt.grid(linestyle="dotted",c="0.8")
    plt.xlabel('time in days')
    plt.xticks(dates,months,rotation = 45)
    plt.ylabel('deltaM')
    plt.show()
      
    """ 
    k(T)*p and transmission rate (betaB) from bird to mosquito 
    """
    betaB_j = sir.biting_rate(T)*0.17*0.165*0.75
    kp_j = sir.biting_rate(T)*0.17
    betaB_c1 = sir.biting_rate(T)*0.06 *0.33 *0.75
    betaB_c2 = sir.biting_rate(T)*0.06 *0.66 *0.75
    betaB_c3 = sir.biting_rate(T)*0.06 *0.96 *0.75
    kp_c = sir.biting_rate(T)*0.06
    
    plt.figure(figsize=(12,4),dpi=300)
    plt.subplot(121)
    plt.plot(T,kp_j,color="black",linestyle ="--",
             label="$Ae.$ $j.$ $japonicus$")
    plt.plot(T,kp_c,color="black",linestyle ="-",
             label="$Culex$ $agg.$")
    plt.grid(linestyle="dotted",c="0.8")
    plt.legend()
    plt.xlabel('Temperature (°C)')
    plt.ylabel('$k(T)p_{B}$',labelpad=10)
    plt.title('a)' , loc='left',fontweight='bold')
        
    plt.subplot(122)
    plt.plot(T,betaB_j, color="black",linestyle ="--",
             label="$Ae.$ $j.$ $japonicus$, $P_{J} = 0.165$")
    plt.plot(T,betaB_c1,color="0.75",linestyle ="-",
             label="$Culex$ $agg.$, $P$ = 0.33")
    plt.plot(T,betaB_c2,color="0.5 ",linestyle ="-",
             label="$Culex$ $agg.$, $P$ = 0.66")
    plt.plot(T,betaB_c3,color="0.0",linestyle ="-",
             label="$Culex$ $agg.$, $P$ = 0.96")
    plt.grid(linestyle="dotted",c="0.8")
    plt.legend()
    plt.xlabel('Temperature (°C)')
    plt.ylabel(r'$\beta_{B}$',labelpad=10)
    plt.title('b)' , loc='left',fontweight='bold') 
    plt.tight_layout(w_pad=3)
    #plt.savefig('betaB_kp.pdf',bbox_inches='tight')
    plt.show()
        
    """ 
    Bird birth as function of the year
    """
    t=np.linspace(1,365,365)
    bB=[]
    Pre_bB=[]
    for tx in t:
        bB.append(sir.bB(tx))
        Pre_bB.append(sir.bB(tx)/0.614)
    X = [125, 135, 145, 156, 166, 176, 186, 196, 206]
    Y = np.array([1.52, 15.27, 33.27, 25.27, 10.27, 5.77, 4.27, 2.77, 1.60])
    
    plt.figure(figsize=(6,4),dpi=300)
    plt.bar(X, Y/1000, width=11,color="0.75", label="observations")
    plt.plot(t, Pre_bB, color="0.5", label = "fit")
    plt.plot(t, bB, color="black", label = "fit * 0.614")
    
    plt.grid(linestyle="dotted",c="0.8")
    plt.legend()
    plt.xticks(dates,months,rotation=0)
    plt.ylabel('$b_{B}$')
    plt.subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=0.15)
    #plt.savefig('MagpieBirth_New.pdf',bbox_inches='tight')
    plt.show()
    
    """ 
    Daylength 
    """
    daylength=[]
    for tx in t:
        daylength.append(sir.daylength(tx))
    plt.figure(figsize=(6,4),dpi=300)
    plt.plot(t,daylength, color="black")
    plt.grid(linestyle="dotted",c="0.8")
    plt.xticks(dates,months,rotation=0)
    plt.ylabel('D',labelpad=10)
    plt.subplots_adjust(bottom=0.15)
    #plt.savefig('Daylength.pdf',bbox_inches='tight')
    plt.show()

    return sir

sir=show_tests()                                