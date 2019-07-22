# ffmpeg -r 60 -f image2 -s 1920x1080 -i movie/img%05d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p movie.mp4

import numpy as np
import matplotlib.pyplot as plt
import glob


path = 'output/output*'
files = glob.glob(path)
counter = 1


for file in sorted(files):
    data = np.loadtxt(file)
    x = data[:, 0]
    rho = data[:, 1]
    rhou = data[:, 2]

    plt.plot(x, rho)
    plt.grid()
    plt.xlabel('x')
    plt.ylabel(r'$\rho$')
    plt.ylim(0, 4)
    plt.savefig('movie/img{:05d}.png'.format(counter))
    plt.clf()
    print('creating img{:05d}.png'.format(counter))
    counter += 1
