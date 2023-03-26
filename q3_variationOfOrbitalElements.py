#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 19:32:58 2023

@author: bhaskar
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import ode
from scipy.integrate import LSODA
from scipy.integrate import solve_ivp

R = 6370 #km, earth's radius
J2 = 0.00108 # Earth's J2
mu = 398600 # Earth's GM

def trueanomalyFROMmeananomaly(M,e):
    #First get eccentric anomaly using iteration
    #set mean anomaly as first guess
    E = M
    tolerance =1e-8
    for i in range(1000):
        f = E -e*math.sin(E)-M
        fprime = 1-e*math.cos(E)
        dE = -f/fprime
        E = E+dE
        
        if abs(dE) < tolerance:
            print("Obtained True anomaly from Mean anomaly")
            break
    #now we have eccentric anomaly
    
    true_anomaly = 2*math.atan2( (1+e)**0.5 * math.sin(E/2),\
                                (1-e)**0.5 * math.cos(E/2))
    return true_anomaly

def dhdt(h,e,r,i,omega,theta):
    return -(1.5*J2*mu*R**2/r**3)*math.sin(i)**2 * math.sin(2*(omega+theta))


def dedt(h,e,r,i,omega,theta):
    return (1.5*J2*mu*R**2/(h*r**3))*( (h**2/(mu*r))*math.sin(theta)*\
           (3*math.sin(i)**2*math.sin(omega+theta)**2 -1) \
           -math.sin(2*(omega+theta))*math.sin(i)**2 *( (2+e*math.cos(theta))\
                                *math.cos(theta)+e))

def dthetadt(h,e,r,i,omega,theta):
    return (h/r**2)+(1.5*J2*mu*R**2/(e*h*r**3))*((h**2/(mu*r))*math.cos(theta)*\
           (3*math.sin(i)**2*math.sin(omega+theta)**2 -1) \
           -math.sin(2*(omega+theta))*math.sin(i)**2 *(2+e*math.cos(theta))\
           *math.sin(theta))

def dOMEGAdt(h,e,r,i,omega,theta):
    return (-3*J2*mu*R**2/(h*r**3))*math.sin(omega+theta)**2 *math.cos(i)


def didt(h,e,r,i,omega,theta):
    return (-0.75*J2*mu*R**2/(h*r**3))*math.sin(2*(omega+theta))*math.sin(2*i)


def domegadt(h,e,r,i,omega,theta):
    return (1.5*J2*mu*R**2/(e*h*r**3))*((h**2/(mu*r))*math.cos(theta)\
           *(1-3*math.sin(i)**2 * math.sin(omega+theta)**2)-(2+e*math.cos(theta))\
           *math.sin(2*(omega+theta))*math.sin(i)**2*math.sin(theta)+2*e\
               *math.cos(i)**2 *math.sin(omega+theta)**2)

def f(t, y):
    #y = [h,e,i,omega,OMEGA,theta]
    #this function collects all the derivatives and pass it to the solver as an
    # array
    dydt = np.zeros(6)
    h=y[0]
    e=y[1]
    r = (y[0]**2/mu)/(1+y[1]*math.cos(y[5]))
    i=y[2]
    omega=y[3]
    theta=y[5]
    dydt[0] = dhdt(h, e, r, i, omega, theta)
    dydt[1] = dedt(h, e, r, i, omega, theta)
    dydt[2] = didt(h, e, r, i, omega, theta)
    dydt[3] = domegadt(h, e, r, i, omega, theta)
    dydt[4] = dOMEGAdt(h, e, r, i, omega, theta)
    dydt[5] = dthetadt(h, e, r, i, omega, theta)
    return dydt

case = 'actual' #test,actual      
#put in 'test', and adjust te total time below to try the test case                                      

if (case=='test'):
    a0 = 8309 #km
    e0 = 0.19629
    i0 = math.radians(28) #radians
    omega0 = math.radians(30)#radians
    OMEGA0 =math.radians(45)#radians
    M0 =math.radians(40)#radians, fraction of the orbit elapsed

if (case=='actual'):
    a0 = 26600 #km
    e0 = 0.74
    i0 = 1.10654 #radians
    omega0 = 5   * math.pi/180 #radians
    OMEGA0 =90   * math.pi/180 #radians
    M0 =10      * math.pi/180 #radians, fraction of the orbit elapsed



totalTime=100*24*3600 #seconds

a=[a0]
e=[e0]
h=[(a[0]*mu*(1-e[0]**2))**0.5] #initial angular momentum
i=[i0]
M=[M0]
omega=[omega0]
OMEGA=[OMEGA0]
theta=[trueanomalyFROMmeananomaly(M[0], e[0])] # initial true anomaly
r=[(h[0]**2/mu)/(1+e[0]*math.cos(theta[0]))] #distance of satellite from Earth
time=[0]

t=0

#solve ivp
# Create the `ode` object with `lsoda` solver
integrator = ode(f).set_integrator('lsoda', rtol=1e-6, atol=1e-9, max_step=100.0)

# Set the initial conditions
integrator.set_initial_value([h[0], e0, i0, omega0, OMEGA0, theta[0]], 0)

#integrate the system
while integrator.successful() and integrator.t < totalTime:
    integrator.integrate(integrator.t+1)
    h.append(integrator.y[0])
    e.append(integrator.y[1])
    i.append(integrator.y[2])
    omega.append(integrator.y[3])
    OMEGA.append(integrator.y[4])
    theta.append(integrator.y[5])
    r.append((h[t]**2/mu)/(1+e[t]*math.cos(theta[t])))
    a.append((h[t]**2/mu)/(1-e[t]**2))
    time.append(integrator.t)
    t+=1


# while time[t]<totalTime:
    
#     de_by_dt = dedt(h[t],e[t],r[t],i[t],omega[t],theta[t])
#     dt = 1.e-8/abs(de_by_dt)
    
#     dh_by_dt = dhdt(h[t],e[t],r[t],i[t],omega[t],theta[t])
#     de_by_dt = dedt(h[t],e[t],r[t],i[t],omega[t],theta[t])
#     dtheta_by_dt = dthetadt(h[t],e[t],r[t],i[t],omega[t],theta[t])
#     dOMEGA_by_dt = dOMEGAdt(h[t],e[t],r[t],i[t],omega[t],theta[t])
#     di_by_dt = didt(h[t],e[t],r[t],i[t],omega[t],theta[t])
#     domega_by_dt = domegadt(h[t],e[t],r[t],i[t],omega[t],theta[t])
    
#     dt = np.min([(h[t]/e[t])*1.e-8/abs(dh_by_dt),1.e-8/abs(de_by_dt),\
#                (i[t]/e[t])*1.e-8/abs(di_by_dt),\
#              (OMEGA[t]/e[t])*1.e-8/abs(dOMEGA_by_dt),\
#                  (omega[t]/e[t])*1.e-8/abs(domega_by_dt)])
#     #print('dt= ',dt)
#     dh = dh_by_dt*dt
#     de= de_by_dt*dt
#     dtheta = dtheta_by_dt*dt
#     dOMEGA = dOMEGA_by_dt*dt
#     di = di_by_dt*dt
#     domega = domega_by_dt*dt
    
#     time.append(time[t]+dt/3600)
#     h.append(h[t]+dh)
#     e.append(e[t]+de)
#     OMEGA.append(OMEGA[t]+dOMEGA)
#     i.append(i[t]+di)
#     omega.append(omega[t]+domega)
    
#     a.append((h[t+1]**2/mu)/(1-e[t+1]**2))
#     theta.append(theta[t]+dtheta)
#     r.append((h[t+1]**2/mu)/(1+e[t+1]*math.cos(theta[t+1])))
    
#     t=t+1
#     if (math.floor(time[t])/time[t]>0.99):
#         print('dt= ',dt)
#         print('time= ',time[t])

time =np.array(time)/(3600*24) #change the time array to show days
plt.figure(1,figsize=(10,4))
plt.plot(time,a)
plt.xlabel('Time (days)')
plt.ylabel('Semi-major axis, a (km)')
plt.savefig('3_a.png',bbox_inches='tight')
plt.figure(2,figsize=(10,4))
plt.plot(time,np.array(i)*180/math.pi)
plt.xlabel('Time (days)')
plt.ylabel('Inclination, i (degrees)')
plt.savefig('3_i.png',bbox_inches='tight')
plt.figure(3,figsize=(10,4))
plt.plot(time,e)
plt.xlabel('Time (days)')
plt.ylabel('Eccentricity, e')
plt.savefig('3_e.png',bbox_inches='tight')
plt.figure(4,figsize=(10,4))
plt.plot(time,np.array(omega)*180/math.pi)
plt.xlabel('Time (days)')
plt.ylabel('Argument of periapsis, $\omega$ (degrees)')
plt.savefig('3_omega.png',bbox_inches='tight')
plt.figure(5,figsize=(10,4))
plt.plot(time,np.array(OMEGA)*180/math.pi)
plt.xlabel('Time (days)')
plt.ylabel('Longitude of ascending node, $\Omega$ (degrees)')
plt.savefig('3_oomega.png',bbox_inches='tight')
plt.show()