# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 18:40:55 2020

@author: Vikram
"""


# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 17:41:25 2020

@author: Vikram
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


N_steps=2000
step=0.01
t = np.arange(0.0, N_steps*step,step)
N=14000000
st=[0.99*N,0,0*N,0*N,0.01*N,0*N,0*N,0*N,0*N,0*N]
time_factor=1
k0=1000 * time_factor
b=0.06
p=1/5  *time_factor
m=1/3*time_factor
v=1/25*time_factor
w=1/3*time_factor
q=0.1*time_factor
vacc=0.0000
lmbda=2*time_factor
r=1/10000*time_factor
alpha=0.0
pt=0.00*time_factor
imm=0 


def model(state,t):
    s,sq,e,eq,iu,iq,ii,r,d,f =state
#    modelling step and impuse functions
# if t > 80 and t <80.2:
#        imm=10
#        alpha=0
#        k0=30
#        print("h")
#    elif t > 80:
#        imm=0
#        alpha=0.5
#        k0=20
#        print("h")
#    else:
#        imm=0
#        alpha=0.7
#        k0=10

#        print(t)
    
    
    fear=(1-f/N)
    govt_action=1-alpha
    k=k0*govt_action*fear
    
    dsdt = -k*b*s*iu/N  - q*k*(1-b)*s*iu/N -   vacc*s +  r*sq
    dsqdt =  q*k*(1-b)*s*iu/N  -  r*sq

    dedt= + (1-q)*k*b*s*iu/N  - p*e  -pt*e
    deqdt= + q*k*b*s*iu/N  -p*eq
    
    diudt = p*e  -(w+v+m)*iu + imm -pt *iu
    diqdt = p*eq  -(w+v+m)*iq
    diidt = w*(iu+iq)  -(v+m)*ii  + pt*(e+iu)
    
    drdt = v * (iu+ii+iq)+vacc*s
    dddt=   m * (iu+ii+iq)
    dfdt=   d - lmbda * f
    
    return dsdt,dsqdt,dedt,deqdt,diudt,diqdt,diidt,drdt,dddt,dfdt

y=odeint(model,st,t)
print(y[-1,:])
plt.figure()
plt.plot(t,y/N)
plt.legend(["Susceptible","Susceptible Q","Exposed","Exposed Q","Infectious U","Infectious Q","Infectious Iso","Removed","Dead","Fear"])

plt.show()
