# ffmpeg -r 60 -f image2 -s 1920x1080 -i moviedisc_iso_grav/img%05d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p moviedisc_iso_grav.mp4
import numpy as np
import matplotlib.pyplot as plt
import glob


path = 'outputdisc_iso_grav/output*'
files = glob.glob(path)
counter = 1
gamma = 7/5

plt.rcParams["figure.figsize"] = [12, 8]

# physical constants - mass in solar masses, length in m
G = 6.67e-11
AU = 1.496e11
k = 1.38e-23
M_solar = 1.9891e30
M_star = M_solar
starDist = 100 * AU
m_h = 1.673e-27
mu = 1.3  # following Facchini et al. 2016
N = 10**9.2
T = 20
gamma = 7. / 5.


for file in sorted(files):
    data = np.loadtxt(file)
    # x = data[:, 0]/AU + 100
    x = data[:, 0]/1000
    rho = data[:, 1]
    rhou = data[:, 2]

    u = rhou/rho
    p = rho * k * T / (mu * m_h)

    '''
    plt.plot(x, rho, '.-')
    plt.grid()
    plt.xlabel(r'R (AU)')
    plt.ylabel(r'$\rho$ $(kgm^{-3})$')
    '''

    plt.plot(x, p, '.-')
    plt.grid()
    plt.xlabel(r'R (km)')
    plt.ylabel(r'$P$ $(Pa)$')
    # plt.xlim(100, 200)

    plt.tight_layout()
    plt.savefig('moviedisc_iso_grav/img{:05d}.png'.format(counter))
    plt.clf()

    print('creating img{:05d}.png'.format(counter))
    counter += 1
