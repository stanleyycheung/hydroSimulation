import numpy as np


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pylab


# dataFin50 = np.loadtxt("1Msol_50AU_1000G0_10pc_finish.dat")
# dataAbund50 = np.loadtxt("1Msol_50AU_1000G0_10pc_abundances.dat")

# dataFin100 = np.loadtxt("1Msol_100AU_1000G0_10pc_finish.dat")
# dataAbund100 = np.loadtxt("1Msol_100AU_1000G0_10pc_abundances.dat")

dataFin100_1p3 = np.loadtxt("FRIED_grid/1p3Msol_100AU_1000G0_20pc_finish.dat")
# dataAbund100_1p3 = np.loadtxt("1p3Msol_100AU_1000G0_20pc_abundances.dat")

# dataFin200 = np.loadtxt("1Msol_200AU_1000G0_10pc_finish.dat")
# dataAbund200 = np.loadtxt("1Msol_200AU_1000G0_10pc_abundances.dat")

au = 1.5e3
mu = 1.3
mH = 1.67e-24

plt.axis(fontsize=24)
fig, axes = plt.subplots(figsize=(10, 7))
axes.set_xlabel('R, AU', fontsize=18)
axes.set_ylabel('log$_{10}$ number density', fontsize=18)

axes2 = axes.twinx()
axes2.set_ylabel('Temperature', fontsize=18)

# axes.plot(dataFin50[:,0]/au, np.log10(dataFin50[:,2]/mu/mH))
# axes.plot(dataFin100[:,0]/au, np.log10(dataFin100[:,2]/mu/mH))
axes.plot(dataFin100_1p3[:, 0]/au, np.log10(dataFin100_1p3[:, 2]/mu/mH), color="blue")
axes2.plot(dataFin100_1p3[:, 0]/au, (dataFin100_1p3[:, 3]), color="red")
# axes.plot(dataFin200[:,0]/au, np.log10(dataFin200[:,2]/mu/mH))

plt.show()
