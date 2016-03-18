# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 18:26:25 2016
richtingen van de stimuli
@author: alban
"""
import matplotlib.pyplot as plt
ax = plt.axes()
import math
pi = math.pi

hratio = 0.75
DisplaySizeX = 640
DisplaySizeY = 480
ms = [DisplaySizeX / 2, DisplaySizeY / 2]
radius = hratio * (ms[1] + 100)
triallength = 24
trialnumber = [1,2,3,4,5,6,7,8]
colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'brown']
print radius

for i in range (0,len(trialnumber)):
    direction = trialnumber[i]
    angle = direction * (2 * pi / 8)
    beginpoint = ms[0] + math.cos(angle) * radius, ms[1] + math.sin(angle) * radius
    endpoint = ms[0] - math.cos(angle) * radius, ms[1] - math.sin(angle) * radius
    chosencolor = colors[i]
    print 'current direction is', direction 
    print 'current color is', chosencolor
    print 'beginpoint is ', beginpoint
    print 'endpoint is ', endpoint    
    if i <= 4:
        plt.plot([beginpoint[0],endpoint[0]],[beginpoint[1],endpoint[1]], chosencolor, 
                 ls = '--',  lw=2)
    else:
        plt.plot([beginpoint[0],endpoint[0]],[beginpoint[1],endpoint[1]], chosencolor, 
                 ls = ':',  lw=2)
    
plt.axis([0,600, 0, 500])
plt.show()