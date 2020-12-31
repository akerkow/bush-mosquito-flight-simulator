#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Visualization of the simulation results saved in log.csv and the locations of 
all mosquitoes and the infectious mosquitoes at the day with the highest number
of infectious mosquitoes
"""

#==============================================================================
# Load modules
#==============================================================================
import Model_Szenario as MS
import pandas as pd
pd.set_option('display.expand_frame_repr', False)
import numpy as np
import pylab as plt
from matplotlib.colors import PowerNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams.update({'font.size': 12})

#==============================================================================
# Get testregion from the last modelrun
#==============================================================================
testregion = MS.testregion
superagents = MS.superagents

#==============================================================================
# Read the data from the log-file that contains the number of birds and mosqui-
# toes in their respective stage of infection for every day; 
# Save the data in a pandas data frame and show a description
#==============================================================================
log=pd.read_csv('log_{}_super{}.csv'.format(testregion,superagents), sep=' ')
print(log.describe())

#==============================================================================
# Visualize all mosquitoes, all birds, infected mosquitoes and infected 
# birds in 4 subplots
#==============================================================================
NM=log.mosquitoes.values
LM=log.LM.values
EM=log.EM.values
SM=log.SM.values
IM=log.IM.values

NB=log.NB.values
SB=log.SB.values
RB=log.RB.values
IB=log.IB.values

years=MS.years
T=[0,365*years]

plt.figure(figsize=(14,5),dpi=300)
ax1=plt.subplot(221)
plt.plot(LM, label="$L_{M}$", linestyle='dashed', color="black")
plt.plot(NM, label="$N_{M}$", linestyle='solid', color="black")
plt.ylabel('mosquitoes')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.grid(linestyle="dotted",c="0.8")
plt.title('a)' , loc='left',fontweight='bold', y=1.08)
plt.legend(loc='center left')

ax2=plt.subplot(222)
plt.plot(IM, label="$I_{M}$", linestyle='solid', color="black")
plt.ylabel('mosquitoes')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.grid(linestyle="dotted",c="0.8")
plt.legend(loc='center left')
plt.title('b)' , loc='left',fontweight='bold', y=1.08)

ax3=plt.subplot(223,sharex=ax1)
plt.plot(NB, label="$N_{B}$", linestyle='solid', color="black")
plt.plot(SB, label="$S_{B}$", linestyle='dashed', color="black")
plt.plot(RB, label="$R_{B}$", linestyle='dotted', color="black")
plt.ylabel('birds')
plt.xlabel('day')
plt.xticks(np.arange(0,365*years,60))
plt.grid(linestyle="dotted",c="0.8")
plt.title('c)' , loc='left',fontweight='bold')
plt.legend(loc="center left")

ax4=plt.subplot(224,sharex=ax2)
plt.plot(IB, label="$I_{B}$", linestyle='solid', color="black")
plt.ylabel('birds')
plt.xlabel('day')
plt.xticks(np.arange(0,365*years,60))
plt.grid(linestyle="dotted",c="0.8")
plt.title('d)' , loc='left',fontweight='bold')
plt.legend(loc="center left")    

plt.tight_layout(pad=2.5,h_pad=2)
plt.savefig('SEIR_results_{}_super{}.pdf'.format(testregion,superagents), 
            bbox_inches='tight')
plt.show()

#==============================================================================
# Show how many mosquitoes have emigrated to the respective cardinal directions
# at each simulation step
#==============================================================================
migrate=log.migrate.values
Emi_N=log.Emi_N.values
Emi_NE=log.Emi_NE.values
Emi_NW=log.Emi_NW.values
Emi_S=log.Emi_S.values
Emi_SE=log.Emi_SE.values
Emi_SW=log.Emi_SW.values
Emi_W=log.Emi_W.values
Emi_E=log.Emi_E.values
others=Emi_NE+ Emi_SE+ Emi_SW+ Emi_NW

sumEmigrations=migrate.max()
sumEmigrators=Emi_N.max()+Emi_S.max()+Emi_W.max()+Emi_E.max()+others.max()
print("A total of {} supermosquitoes tried to emigrate".format(sumEmigrators)+
      " {} times alltogether (some several times)!".format(sumEmigrations))

legendlist=["North","East","South","West","Others"]

plt.figure(figsize=(6,3.5),dpi=300)
plt.plot(Emi_N, linestyle="solid", c="0.8")
plt.plot(Emi_E, linestyle="solid", c="0.35")

plt.plot(Emi_S, linestyle="dotted", c="0.0")
plt.plot(Emi_W, linestyle="dotted", c="0.6")
plt.plot(others, linestyle="dashed", c="0.0")
plt.legend(legendlist)
plt.xlabel("Day")
plt.ylabel("Supermosquitoes")
plt.minorticks_on()
plt.grid(b=True, which='minor', color='0.95', linestyle='-')
plt.grid(b=True, which='major', color='0.85', linestyle='-')
plt.savefig('Emigrations_{}_super{}.pdf'.format(testregion,superagents),
            bbox_inches='tight')
plt.show()

#==============================================================================
# Show the spatial distribution of (i) all mosquitoes and (ii) infectious mos-
# quitoes at the day with the highest number of infectious mosquitoes
#==============================================================================
m=np.load('matrixm_HI_{}_super{}.npy'.format(testregion,superagents))
i=np.load('matrixi_HI_{}_super{}.npy'.format(testregion,superagents))
inf=np.int16(i)

norm = PowerNorm(0.4, vmin=0, vmax=np.max(m)+1)
color = plt.cm.inferno
monochrom = "gray_r"

fig, (ax1, ax2) = plt.subplots(figsize=(20,8),ncols=2, dpi=300)

NM = ax1.imshow(m,
                cmap=monochrom,
                norm = norm)

divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="5%", pad=0.2)
fig.colorbar(NM, cax=cax1)
plt.xticks()
plt.yticks()
ax1.set_title('a)' , loc='left',fontweight='bold', y=1.08)
ax1.set_title('All active mosquitoes ($N_{M}$)', y=1.08)


IM = ax2.imshow(inf, cmap=color, vmin=0,vmax=np.max(inf))
divider = make_axes_locatable(ax2)
cax2 = divider.append_axes("right", size="5%", pad=0.2)
fig.colorbar(IM, cax=cax2, ticks=range(np.max(inf)+1))


plt.xticks()
plt.yticks()
ax2.set_title('b)' , loc='left',fontweight='bold', y=1.08)
ax2.set_title('Infectious mosquitoes ($I_{M}$)', y=1.08)
plt.tight_layout(pad=4, w_pad=8)

plt.savefig('MosquitoDistributions_{}_super{}.pdf'.format(testregion,
            superagents), bbox_inches='tight')
plt.show()