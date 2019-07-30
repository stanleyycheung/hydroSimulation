import numpy as np


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pylab

#dataFin100_1p3=np.loadtxt("1p3Msol_100AU_1000G0_20pc_finish.dat")
#dataAbund100_1p3=np.loadtxt("1p3Msol_100AU_1000G0_20pc_abundances.dat")
dataAbund100=np.loadtxt("1Msol_100AU_1000G0_20pc_abundances.dat")

au=1.5e3
mu=1.3
mH=1.67e-24

plt.axis(fontsize=24)
fig, axes= plt.subplots(figsize=(11,7))
axes.set_xlabel('R, AU', fontsize=18)
axes.set_ylabel('log$_{10}$ species abundance', fontsize=18)

axes.tick_params(labelsize=17)


#axes2=axes.twinx()
#axes2.set_ylabel('Temperature', fontsize=18)

axes.set_ylim(-15, 0)

#axes.plot(dataFin50[:,0]/au, np.log10(dataFin50[:,2]/mu/mH))
#axes.plot(dataFin100[:,0]/au, np.log10(dataFin100[:,2]/mu/mH))
axes.plot(dataAbund100[:,0]/au, np.log10(dataAbund100[:,4]), linewidth=4, label="H$_2$")
axes.plot(dataAbund100[:,0]/au, np.log10(dataAbund100[:,5]), linewidth=4, label="H")
axes.plot(dataAbund100[:,0]/au, np.log10(dataAbund100[:,7]), linewidth=4, label="O")
axes.plot(dataAbund100[:,0]/au, np.log10(dataAbund100[:,11]), linewidth=4, label="C")
axes.plot(dataAbund100[:,0]/au, np.log10(dataAbund100[:,10]), linewidth=4, label="C$^+$")
axes.plot(dataAbund100[:,0]/au, np.log10(dataAbund100[:,13]), linewidth=4, label="CO")
axes.plot(dataAbund100[:,0]/au, np.log10(dataAbund100[:,20]), linewidth=4, label="HCO$^+$")
axes.legend()


#axes2.plot(dataFin100[:,0]/au,(dataFin100[:,3]), color="red")
#axes.plot(dataFin200[:,0]/au, np.log10(dataFin200[:,2]/mu/mH))
plt.savefig("Abundances.pdf")
plt.show()
