# ffmpeg -r 60 -f image2 -s 1920x1080 -i moviedisc_iso_grav/img%05d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p moviedisc_iso_grav.mp4
import numpy as np
import matplotlib.pyplot as plt
import glob


path = 'outputdisc_iso_grav/output*'
files = glob.glob(path)
counter = 1
gamma = 7/5

plt.rcParams["figure.figsize"] = [24, 8]

# physical constants - mass in solar masses, length in m
G = 6.67e-11
AU = 1.496e11
k = 1.38e-23
M_solar = 1.9891e30
M_star = M_solar
starDist = 100 * AU
# m_h = 1.673e-27
# mu = 1.3  # following Facchini et al. 2016
m_air = 28.97/1000/(6.02e23)
N = 10**9.2
T = 288.16
gamma = 7. / 5.


def solution(z):
    # plot the analytic solution
    return 101325*np.exp(-z / (7.2))


'''
data = np.loadtxt(files[0])
x = data[:, 0]/1000 + 6.371e3
x_plot = data[:, 0]/1000
plt.plot(x_plot, solution(x_plot))
plt.show()

'''
for file in sorted(files):
    data = np.loadtxt(file)
    # x = data[:, 0]/AU + 100
    x = data[:, 0]/1000 + 6.371e6
    rho = data[:, 1]
    rhou = data[:, 2]

    u = rhou/rho
    p = rho * k * T / (m_air) / 1000

    # plt.plot(x, solution(x))
    plt.subplot(1, 2, 1)
    plt.plot(x, p, '.-')
    plt.grid()
    plt.xlabel(r'R (km)')
    plt.ylabel(r'$P$ $(kPa)$')
    # plt.xlim(100, 200)

    plt.subplot(1, 2, 2)
    plt.plot(x, u, '.-')
    plt.grid()
    plt.xlabel(r'R (km)')
    plt.ylabel(r'$u$ $(ms^{-1})$')

    plt.tight_layout()
    plt.savefig('moviedisc_iso_grav/img{:05d}.png'.format(counter))
    plt.clf()

    print('creating img{:05d}.png'.format(counter))
    counter += 1
