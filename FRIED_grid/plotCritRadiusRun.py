import numpy as np
import os
import matplotlib.pyplot as plt
import sys

kerg=1.380626e-16
mu=1.3
mH=1.67e-24
pc=3.08e8

data=np.loadtxt("../radial0018.dat")

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
    if(movingAverage >0.0):
        print "critical point on grid", data[i,0]/au, movingAverage
    #    if(data[i,0]/au < 0.95*data[len(data[:,0],0]):
    #       if((v**2*mu*mH/kerg) - T - n*(dT/dn)>0.0):

    lastn=n
    lastT=T
    i=i+1



plt.xlabel("R, AU")
plt.ylabel("LHS")
#plt.xlim(0,0.4)
plt.ylim(-400,400)
plt.axhline(y=0.0, color='r', linestyle='-')
plt.plot(RtoPlot[:], LHS[:])

#plt.semilogy(RtoPlot[:], LHS[:])

plt.show()
