# ffmpeg -r 60 -f image2 -s 1920x1080 -i moviedisc/img%05d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p moviedisc.mp4
import numpy as np
import matplotlib.pyplot as plt
import glob


path = 'outputdisc/output*'
files = glob.glob(path)
counter = 1
gamma = 7/5

plt.rcParams["figure.figsize"] = [24, 8]

interval = 1
files = files[::interval]


for file in sorted(files):
    data = np.loadtxt(file)
    x = data[:, 0]
    rho = data[:, 1]
    rhou = data[:, 2]
    rhoe = data[:, 3]
    rhov = data[:, 4]

    u = rhou/rho
    v = rhov/rho
    e_th = rhoe/rho - u**2/2
    p = (gamma-1)*rho*e_th

    plt.subplot(1, 3, 1)
    plt.plot(x, rho, '.-')
    plt.grid()
    plt.xlabel(r'R (m)')
    plt.ylabel(r'$\rho$ $(kgm^{-3})$')

    plt.subplot(1, 3, 2)
    plt.plot(x, p, '.-')
    plt.grid()
    plt.xlabel(r'R (m)')
    plt.ylabel(r'$p$ $(Pa)$')

    plt.subplot(1, 3, 3)
    plt.plot(x, v, '.-')
    plt.grid()
    plt.xlabel(r'R (m)')
    plt.ylabel(r'$v$ $(ms^{-1})$')

    plt.tight_layout()
    plt.savefig('moviedisc/img{:05d}.png'.format(counter))
    plt.clf()

    print('creating img{:05d}.png'.format(counter))
    counter += 1
