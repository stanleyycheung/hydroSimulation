import numpy as np


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pylab


dataFin50=np.loadtxt("1Msol_50AU_1000G0_10pc_finish.dat")
dataAbund50=np.loadtxt("1Msol_50AU_1000G0_10pc_abundances.dat")

dataFin100=np.loadtxt("1Msol_100AU_1000G0_20pc_finish.dat")
dataFin100B=np.loadtxt("1Msol_100AU_1000G0_20pc_finishB.dat")
dataFin100C=np.loadtxt("1Msol_100AU_1000G0_20pc_finishC.dat")
dataAbund100=np.loadtxt("1Msol_100AU_1000G0_20pc_abundances.dat")



dataFin100_1p3=np.loadtxt("1p3Msol_100AU_1000G0_20pc_finish.dat")
dataAbund100_1p3=np.loadtxt("1p3Msol_100AU_1000G0_20pc_abundances.dat")

dataFin200=np.loadtxt("1Msol_200AU_1000G0_10pc_finish.dat")
dataAbund200=np.loadtxt("1Msol_200AU_1000G0_10pc_abundances.dat")

au=1.5e3
mu=1.3
mH=1.67e-24

#plt.tick_params(labelsize=24)
#plt.axis(fontsize=28)
fig, axes= plt.subplots(figsize=(11,7))
axes.set_xlabel('R, AU', fontsize=18)
axes.set_ylabel('log$_{10}$(n), gcm$^{-3}$', fontsize=18)


axes2=axes.twinx()
axes2.set_ylabel('log$_{10}$($\dotM$), M$_\odot$yr$^{-1}$', fontsize=18)
axes.text(110, 14., "Disc outer edge", fontsize=18)

axes.tick_params(labelsize=16)
axes2.tick_params(labelsize=16)



#axes.plot(dataFin50[:,0]/au, np.log10(dataFin50[:,2]/mu/mH))
#axes.plot(dataFin100[:,0]/au, np.log10(dataFin100[:,2]/mu/mH))
#axes.plot(dataFin100_1p3[:,0]/au, np.log10(dataFin100_1p3[:,2]/mu/mH), color="blue")
#axes2.plot(dataFin100_1p3[:,0]/au,np.log10(dataFin100_1p3[:,15]), color="red")
#axes.plot(dataFin200[:,0]/au, np.log10(dataFin200[:,2]/mu/mH))

axes2.set_ylim(-6.4, -6.0)

#dataFin100[:,2]*

lns1=axes.plot(dataFin100[:,0]/au, np.log10(dataFin100[:,2]/mu/mH), linewidth=3, color="blue", label="Number density")
lns2=axes2.plot(dataFin100[:,0]/au,np.log10(dataFin100[:,15]), linewidth=3, color="red", label="Mass loss rate")

#axes.plot(dataFin100B[:,0]/au, np.log10(dataFin100B[:,2]/mu/mH),  color="blue")
axes2.plot(dataFin100B[:,0]/au,np.log10(dataFin100B[:,15]), linewidth=3, color="red")

#axes.plot(dataFin100C[:,0]/au, np.log10(dataFin100C[:,2]/mu/mH), linewidth=3, color="blue")
axes2.plot(dataFin100C[:,0]/au,np.log10(dataFin100C[:,15]), linewidth=3, color="red")

lns=lns1+lns2
labs=[l.get_label() for l in lns]
axes.legend(lns, labs, loc=0)

plt.axvline(x=100., linewidth=3, color='black')

#axes.legend(loc=2)
#axes2.legend(loc=1)

plt.savefig("Example_n_Mdot.pdf")
plt.show()
