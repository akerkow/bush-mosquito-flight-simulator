#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The movements of the infectious "super mosquitoes" (IM) are shown here in a 
video. 
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
DaylyInfs_365 = []
for array in sorted(glob.glob('DaylyInf_{}_*.npy'.format(testregion)),
                    key=os.path.getmtime):
    print(array)#check order of reading the files
    DaylyInfs_365.append(np.load(array))

#==============================================================================
# Remove days before the infectious period starts to make the video shorter
#==============================================================================
DaylyInfs = DaylyInfs_365.copy()
daycount = 0
while np.sum(DaylyInfs[0])==0:
    DaylyInfs.pop(0)
    daycount += 1
    #print(daycount, " removed")
print("Start date of the infectious period: ",daycount)

#==============================================================================
# Remove days after the infectious period ended to make the video shorter:
#==============================================================================
while np.sum(DaylyInfs[-1])==0:
    DaylyInfs.pop(-1)
    #print("removed last item")
print("Period with infectious mosquitoes: ",len(DaylyInfs))

#==============================================================================
# Animate the infectious mosquitoes
# -> Uncommend "ani.save..." to save the animation (this may take a minute)
#==============================================================================
maxI = np.max(DaylyInfs)
cmap = "inferno"

fig, ax = plt.subplots(figsize=(10,10),dpi=150)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)

im = ax.imshow(DaylyInfs[0],cmap=cmap,vmin=0,vmax=maxI) # axesimage object

def updatefig(i): 
    im.set_array(DaylyInfs[i]) # set data in the axesimage object
    ax.set_title("Day {}".format(daycount+i+1))
    return [im] # return artists set

plt.colorbar(im, cax=cax, orientation='vertical', ticks=range(maxI+1))
  
ani = animation.FuncAnimation(fig, updatefig, frames = len(DaylyInfs), 
                              interval=300, repeat = False, 
                              blit= False, 
                              cache_frame_data = True)
#ani.save('Infectious_{}_super{}.mp4'.format(testregion, supermosquitoes),
#         fps=2.5, bitrate= -1, dpi=300)
plt.show()
