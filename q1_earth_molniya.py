#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 17:59:01 2023

@author: bhaskar
"""
import numpy as np
import matplotlib.pyplot as plt
import math

R = 6370 #km, earth's radius
J2 = 0.00108 # Earth's J2
mu = 398600 # Earth's GM

def omegaDot(meanMotion,semimajor,inclination,eccentricity):
    return 0.75*meanMotion*J2*((R/semimajor)**2)*\
        ((5*math.cos(inclination)**2 -1)/(1-eccentricity**2)**2)
        
def OmegaDot(meanMotion,semimajor,inclination,eccentricity):
    return -1.5*meanMotion*J2*((R/semimajor)**2)*\
        ((math.cos(inclination))/(1-eccentricity**2)**2)


T = 24*3600/3 # seconds, time period is 1/3rd of a day

a = (T*mu**0.5/(2*math.pi))**(2/3) # km, semi major axis

n = 2*math.pi/T # mean motion

rp = np.arange(600+ R,a,100)# range of perigee radius(from center of earth)
ra = 2*a - rp #corresponding range of apogees
e = (ra-rp)/(ra+rp) #list of eccentricities for the above data

#for omegaDot to be zero, we must have 5 cos^2 i -1 = 0
i=np.empty(2)
i[0] = math.acos(1/5**0.5) #~63.45 degrees
i[1] = math.acos(-1/5**0.5)
print('i[0] = ', i[0]*180/math.pi)
print('i[1] = ', i[1]*180/math.pi)

OmegaDotValues_positive_cosi = OmegaDot(n,a,i[0],e)
OmegaDotValues_negative_cosi = OmegaDot(n,a,i[1],e)

plt.figure(1)
plt.plot(rp-R,OmegaDotValues_positive_cosi,label='$cos\ i = 1/\sqrt{5}$')
plt.plot(rp-R,OmegaDotValues_negative_cosi,label='$cos\ i = -1/\sqrt{5}$')
plt.legend()
plt.xlabel('Perigee altitude (km)')
plt.ylabel(r'$\dot{\bar{\Omega}}$')
plt.title(r'Plot of $\dot{\bar{\Omega}}$ for a range of perigee altitudes')
# Show the plot
plt.savefig('1_OMEGAdot.png',bbox_inches='tight')
plt.show()
#plt.yscale("log")

#put the perigee of choice and obtain rest of the values
choice_rp = a#3000+R
choice_ra = 2*a - choice_rp
choice_e=(choice_ra-choice_rp)/(choice_ra+choice_rp)
choice_OMEGAdot = OmegaDot(n,a,i[0],choice_e)