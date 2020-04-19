# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 17:41:25 2020

@author: Vikram
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


N_steps=500
step=1
t = np.arange(0.0, N_steps*step,step)
N=14000000

beta=1
sigma=1/3
gamma=1/5
print(beta/gamma)
def model(state,t):
    s,e,i,r =state
    
    
    #dxdt = 2*x+y
    #dydt =  x+2*y#-0.1*np.sin(0.5*t)
    
    dsdt = -beta*s*i/N
    dedt= + beta*s*i/N-sigma*e
    didt = sigma*e-gamma*i#-0.1*np.sin(0.5*t)
    drdt = gamma*i #dxdt = 0.01*(-x**3-mu*x+lambdaa)-kappa*y
    #dydt = 0.01*(x-y)/tau#-0.1*np.sin(0.5*t)
    return dsdt,dedt,didt,drdt
st=[(1)*N,0*N,1,(0)*N]
y=odeint(model,st,t)/N
fig,ax=plt.subplots(ncols=2,nrows=1)
ax[0].plot(t[0:-29+90],y[29:90,:])
ax[0].legend(["Susceptible","Exposed","Infectious","Removed"])
ax[0].set_xlabel("Days")
ax[0].set_ylabel("Fraction")
ax[0].set_xlim([0,60])
ax[0].set_ylim([0,1])
ax[1].set_aspect("equal","box")

for i in range(11):
    init_s=0.1*i
    st=[(init_s)*N,0*N,1,(1-init_s)*N]
    y=odeint(model,st,t)/N
    ax[1].plot(y[:,0],y[:,2],color="k",linewidth=0.5)
    
    
#    plt.tightplot()

ax[1].set_xlabel("Susceptible Fraction")
ax[1].set_ylabel("Infective Fraction")
ax[1].plot([1,0],[0,1],color="k",linewidth=0.5)
ax[1].set_aspect("equal","box")
ax[1].set_xlim([0,1])
ax[1].set_ylim([0,1])
#plt
fig.tight_layout()
plt.savefig("phase.svg")

#plt.subplot(2,2,2)
#plt.plot(y[:,0],y[:,1])
#
#plt.subplot(2,2,3)
#plt.plot(y[:,0],y[:,2])
#
#plt.subplot(2,2,4)
#plt.plot(y[:,0],y[:,1])

plt.show()
