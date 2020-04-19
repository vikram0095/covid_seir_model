import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
fig, ax = plt.subplots()
N=51
N_steps=50
d=7
no_of_con=5
decay_rate=3
midN=math.floor((N-1)/2)
M=np.zeros([N_steps,N,N])
M[0,midN,midN]=10

for iter in range(1,N_steps):
    for i in range(d,N-d):
        for j in range(d,N-d):
            M[iter,i,j]=max(0,M[iter-1,i,j]-decay_rate);
            for  con in range(no_of_con):
                i_inf=math.floor(np.random.rand(1)*(2*d+1))-d
                j_inf=math.floor(np.random.rand(1)*(2*d+1))-d
                if np.random.rand(1)>0.5:
                    M[iter,i,j]=(M[iter,i,j]+M[iter-1,i+i_inf,j+j_inf])
#                if M[iter,i,j]>100:
#                    M[iter,i,j]=100
#                print(M[i-1:i+2,j-1:j+2])
mat=ax.matshow(M[1,d:-d,d:-d],vmax=100,vmin=0,cmap=plt.get_cmap('hot_r'))
cbar=plt.colorbar(mat)
                
def animate(frame):
    mat.set_array(M[frame,d:-d,d:-d])

ani = animation.FuncAnimation(
    fig, animate, frames=N_steps, interval=50)

plt.show()

