# ffmpeg -r 120 -f image2 -s 1920x1080 -i moviesod/img%05d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p moviesod2.mp4

import numpy as np
import matplotlib.pyplot as plt
import glob

path = 'outputsod/output*'
files = glob.glob(path)
counter = 1
gamma = 7/5

plt.rcParams["figure.figsize"] = [8, 8]


for file in sorted(files):
    data = np.loadtxt(file)
    x = data[:, 0]
    rho = data[:, 1]
    rhou = data[:, 2]
    rhoue = data[:, 3]

    u = rhou/rho
    e_th = rhoue/rho - u**2/2
    p = (gamma-1)*rho*e_th

    plt.subplot(2, 2, 1)
    plt.plot(x, rho, '.-')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel(r'$\rho$')
    plt.ylim(0, 1.1)

    plt.subplot(2, 2, 2)
    plt.plot(x, u, '.')
    plt.grid()
    plt.ylim(0, 0.85)
    plt.xlabel('x')
    plt.ylabel(r'$u$')

    plt.subplot(2, 2, 3)
    plt.plot(x, p, '.-')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel(r'$p$')

    plt.subplot(2, 2, 4)
    plt.plot(x, e_th, '.')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel(r'$e_{th}$')
    plt.ylim(1.8, 2.6)

    plt.tight_layout()
    plt.savefig('moviesod/img{:05d}.png'.format(counter))
    plt.clf()

    print('creating img{:05d}.png'.format(counter))
    counter += 1
