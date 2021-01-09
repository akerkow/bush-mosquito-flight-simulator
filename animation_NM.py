#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The movements of all (the non-infectious and infectious) super-mosquitoes 
(SM+EM+IM) are shown here in a video. One super-mosquito in the view bundles 
10, 100 or 1000 mosquito individuals as specified in "Model_Szenario.py". 
"""

#==============================================================================
# Load modules
#==============================================================================
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import PowerNorm
import Model_Szenario as MS

#==============================================================================
# Get testregion from the last modelrun
#==============================================================================
testregion = MS.testregion
supermosquitoes = MS.superagents

#==============================================================================
# Read the spatial outputs (numpy arrays) from the last modelrun,
# read them in the order like they were dumped by using the time stamp of the
# file
#==============================================================================
DaylyMs = []
for array in sorted(glob.glob('DaylyM_{}_*.npy'.format(testregion)),
                    key=os.path.getmtime):
    #print(array)#check order of reading the files
    DaylyMs.append(np.load(array))
    
#==============================================================================
# Animate all mosquitoes
# -> Uncommend "ani.save..." to save the animation (this may take a minute)
#==============================================================================
cmap="inferno"
maxM = np.max(DaylyMs)
norm = PowerNorm(0.4, vmin=0, vmax=maxM)

fig, ax = plt.subplots(figsize=(10,10),dpi=150)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05) 
im = ax.imshow(DaylyMs[0], cmap=cmap, norm = norm) # axesimage object

def updatefig(i): 
    im.set_array(DaylyMs[i])# set data in the axesimage object
    ax.set_title("Day {}".format(i+1))
    return [im] # return artists set    

fig.colorbar(im,cax=cax, orientation='vertical')
    
ani = animation.FuncAnimation(fig, updatefig, frames = len(DaylyMs), 
                              interval=300, repeat = False, 
                              blit= False, 
                              cache_frame_data = True)
ani.save('DaylyMs_{}_super{}.mp4'.format(testregion, supermosquitoes),
         fps=2.5, bitrate= -1, dpi=300)
plt.show()
