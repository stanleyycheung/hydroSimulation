import numpy as np
import os
import matplotlib.pyplot as plt
import sys

kerg=1.380626e-16
mu=1.3
mH=1.67e-24
pc=3.08e8

data=np.loadtxt("../inputfile")

au=1.5e3

i=1
lastT=data[0,2]
lastn=data[0,1]/mu/mH
LHS=[]
LHSarray=[]
RtoPlot=[]
movingAverage=0.0
counter=0
while i < len(data[:,0]):
    counter=counter+1
    n = data[i,1]/mu/mH
    T=data[i,2]
    v=data[i,4]*1.e5
    dT=T-lastT
    dn=n-lastn
    RtoPlot.append(data[i,0]/au)
    LHS.append((v**2*mu*mH/kerg) - T - n*(dT/dn))    
    movingAverage=movingAverage+(((v**2*mu*mH/kerg) - T - n*(dT/dn))/counter)
    if(LHS[-1] > 0.0): 
        print "critical point on grid", data[i,0]/au, movingAverage
    #    if(data[i,0]/au < 0.95*data[len(data[:,0],0]):
    #       if((v**2*mu*mH/kerg) - T - n*(dT/dn)>0.0):

    lastn=n
    lastT=T
    i=i+1


fig, axes= plt.subplots(figsize=(11,7))
axes.set_xlabel('R, AU', fontsize=18)
axes.set_ylabel('LHS of equation 9', fontsize=18)

axes.tick_params(labelsize=17)


axes.set_ylim(-400,400)
axes.axhline(y=0.0, linewidth=3, color="black")
axes.plot(RtoPlot[:], LHS[:], linewidth=3, color="blue")

#plt.semilogy(RtoPlot[:], LHS[:])

plt.savefig("criticalRad.pdf")
plt.show()
