import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


N_steps=3500
step=0.01
t = np.arange(0.0, N_steps*step,step)
N=14000000
#st=[0.99*N,0*N,0*N,0.01*N,0*N,0*N,0*N,0*N]
k=10
b=.5
p=1/5
m=1/10
v=1/5
w=1/10
q=0.1

def model(state,t):
#    s,e,eq,iu,iq,ii,r,d =state
    N_districts=3
    district_closeness=np.random.rand(N_districts,N_districts)
    s,e,eq,iu,iq,ii,r,d=np.zeros([8,N_districts])
    dsdt,dedt,deqdt,diudt,diqdt,diidt,drdt,dddt=np.zeros([8,N_districts])
    for s_unit in range(N_districts):
        s[s_unit],e[s_unit],eq[s_unit],iu[s_unit],iq[s_unit],ii[s_unit],r[s_unit],d[s_unit] =state[s_unit*8:(s_unit+1)*8]
#    print(s/N)
    for s_unit in range(N_districts):
        fear=1+0*((1-d[s_unit]/N)**10)

        dsdt[s_unit] = -k*b*s[s_unit]*iu[s_unit]/N*fear
        dedt[s_unit]= + k*b*s[s_unit]*iu[s_unit]/N*fear-p*e[s_unit]
        deqdt[s_unit]= + q*k*b*s[s_unit]*iu[s_unit]/N*fear-p*eq[s_unit]
        imm_fac=1
        emm_fac=0
#        for dis in range(N_districts):
#            if s_unit!=dis:
#                imm_fac=imm_fac+iu[dis]*district_closeness[dis,s_unit]*0.00001
#                emm_fac=emm_fac+iu[s_unit]*district_closeness[s_unit,dis]*0.0001
                
        diudt[s_unit] = p*e[s_unit]-(w+v+m)*iu[s_unit]-emm_fac+imm_fac
        diqdt[s_unit] = p*eq[s_unit]-(w+v+m)*iq[s_unit]
        diidt[s_unit] = w*(iu[s_unit]+iq[s_unit])-v*ii[s_unit]-m*ii[s_unit]
        
        drdt[s_unit] = v*(iu[s_unit]+ii[s_unit]+iq[s_unit])
        dddt[s_unit]=   m*(iu[s_unit]+ii[s_unit]+iq[s_unit])
    ret=[]
    for s_unit in range(N_districts):
        ret.append(dsdt[s_unit])
        ret.append(dedt[s_unit])
        ret.append(deqdt[s_unit])
        ret.append(diudt[s_unit])
        ret.append(diqdt[s_unit])
        ret.append(diidt[s_unit])
        ret.append(drdt[s_unit])
        ret.append(dddt[s_unit])
    return ret

st=[0.9*N,0*N,0*N,0.1*N,0*N,0*N,0*N,0*N   ,
    1*N,0*N,0*N,0*N,0*N,0*N,0*N,0*N,  
    1*N,0*N,0*N,0*N,0*N,0*N,0*N,0*N]

y=odeint(model,st,t)/N
#print(y[-1,:])
plt.subplot(3,1,1)
plt.plot(t,y[:,0:8])

plt.subplot(3,1,2)
plt.plot(t,y[:,8:16])

plt.subplot(3,1,3)
plt.plot(t,y[:,16:24])
plt.legend(["Susceptible","Exposed","Exposed Q","Infectious U","Infectious Q","Infectious Iso","Removed"])

plt.show()
