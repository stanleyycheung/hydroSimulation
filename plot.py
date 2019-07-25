import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt('outputsod/output01248.dat')
x = data[:, 0]
rho = data[:, 1]

plt.plot(x, rho, '.-')
plt.xlabel('x')
plt.ylabel(r'$\rho$')
plt.grid()
plt.show()
