#!/usr/bin/env python3
import numpy as np
import pynbody as pyn
from scipy import ndimage
from matplotlib import rc, ticker
import matplotlib.pyplot as plt

def annotate_lines(loc_labels, slope, color):
    angle = np.arctan(slope)*180./np.pi
    for label, loc in loc_labels.items():
        trans_angle = plt.gca().transData.transform_angles(np.array((angle,)),
                np.array((10,10)).reshape((1,2)))[0]
        plt.text(loc[0], loc[1], label, rotation=trans_angle, color=color, rotation_mode='anchor')

if __name__ == "__main__":
    font = {'size': 14}
    rc('font', **font)
    rc('axes', grid=False)
    rc('axes', facecolor='white')
    # Plot lines of constant cooling time
    coolrates = np.genfromtxt('coolrates.dat', skip_header=2, names=True)
    Edot_per_erg = coolrates['Edot']/(coolrates['nH']/0.74*1.672e-24)
    tcool = -coolrates['E']/(Edot_per_erg*3.155e7)
    tcool[np.where(tcool < 0)] = np.nan #To keep contours from joining up
    plt.contour(coolrates['nH'][::61], coolrates['T'][:61], np.reshape(tcool, (81,61)).T, 
            (1e4, 1e6, 1e8), linewidth=50, colors='Blue')
    plt.gcf().set_size_inches(11,8.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-7,1e4)
    plt.ylim(1e2,1e9)
    plt.ylabel('T $(K)$')
    plt.xlabel(r'$\rho (m_p /cm^3)$')
    plt.text(2e-1, 2e8, "$10^8yr$", color="Blue", rotation='60')
    plt.text(2e1, 2e8, "$10^6yr$", color="Blue", rotation='60')
    plt.text(1e3, 2e8, "$t_{cool}=10^4yr$", color="Blue", rotation='60')

    # Plot simulation data
    pdata = np.load('g1536_gas.npz')
    plt.hexbin(pdata['rho'], pdata['temp'], xscale='log', yscale='log', bins='log', extent=(-7,4,2,9),
            mincnt=2, gridsize=200, cmap='binary')
    # Add interesting density values
    plt.axvline(1,color='r',linestyle=':')
    plt.text(1.2,5e8,'Mean Milky Way ISM', color='r', rotation='270')
    plt.axvline(5e-6,color='r',linestyle=':')
    plt.text(6e-6,5e8,r'$\rho_{crit}\; (\Omega = 1)$', color='r', rotation='270')
    plt.axvline(0.308*5e-6,color='r',linestyle=':')
    plt.text(0.308*6e-6,5e8,r'$\Omega_m = 0.308$', color='r', rotation='270')
    plt.axvline(0.045*5e-6,color='r',linestyle=':')
    plt.text(0.045*6e-6,5e8,r'$\Omega_B = 0.045$', color='r', rotation='270')
    for i in range(-9,8,3):
        plt.plot(np.power(10.,np.arange(-10,10)), np.power(10.,i+np.arange(15,-5,-1)), '--', color='grey')
    plt.text(2e-5,5e3,r"Isobar", rotation=-47, color='grey')
    for i in range(-16,4,3):
        plt.plot(np.logspace(-8,5), np.power(np.logspace(i-8,i+5), -0.66), linestyle='dashdot', color='cyan')
    plt.text(2e-6,3e6,r"$t_{Courant}=$ const", rotation=-36, color='cyan')
    for i in range(3,19,3):
        plt.plot(np.logspace(-8,5), np.power(np.logspace(i-8,i+5), 2./3), 'k:')
    plt.text(2e-4,5e7,r"Adiabat", rotation=42, color='k')
    # Lines of Constant Jeans Mass
    plt.plot(np.logspace(-6,5), np.power(np.square(1e5/3e4)*np.logspace(-3,8), 1./3), 'm', linestyle='-')
    plt.plot(np.logspace(-6,5), np.power(np.square(1e6/3e4)*np.logspace(-3,8), 1./3), 'm', linestyle='-')
    plt.plot(np.logspace(-6,5), np.power(np.square(1e7/3e4)*np.logspace(-3,8), 1./3), 'm', linestyle='-')
    plt.plot(np.logspace(-6,5), np.power(np.square(1e8/3e4)*np.logspace(-3,8), 1./3), 'm', linestyle='-')
    anns = {r"$10^5 M_\odot$":(1.5e3,3e2), 
            r"$M_J = 10^6 M_\odot$":(3e2,9e2), 
            r"$M_J = 10^7 M_\odot$":(3e2,4.5e3), 
            r"$M_J = 10^8 M_\odot$":(3e2,2.2e4)}
    annotate_lines(anns, 1./3, 'm')
    # Lines of Constant Jeans Length
    plt.plot(np.logspace(-10,10), np.square(1e1/12.3)*np.logspace(-10,10), '-', color='green')
    plt.plot(np.logspace(-10,10), np.square(1e2/12.3)*np.logspace(-10,10), '-', color='green')
    plt.plot(np.logspace(-10,10), np.square(1e3/12.3)*np.logspace(-10,10), '-', color='green')
    plt.plot(np.logspace(-10,10), np.square(1e4/12.3)*np.logspace(-10,10), '-', color='green')
    plt.text(9e-5,5e2,r"$\lambda_{J}=10^4pc$", rotation=48, color='green')
    plt.text(9e-3,5e2,r"$\lambda_{J}=10^3pc$", rotation=48, color='green')
    plt.text(9e-1,5e2,r"$\lambda_{J}=100pc$", rotation=48, color='green')
    plt.text(1e2,5e2,r"$\lambda_{J}=10pc$", rotation=48, color='green')
    # Show sound speed as second y axis
    psound = plt.twinx()
    psound.set_ylabel('$c_s (km/s)$')
    psound.set_yticks([0.078, 0.194, 0.300, 0.416, 0.522, 0.638, 0.744, 0.861])
    psound.set_yticklabels(['3','10', '30', '100', '300', '1,000', '3,000', '10,000'])
    plt.figtext(0.71,0.04,"BW Keller's Reference Phase Diagram", fontsize=8)
    plt.figtext(0.825,0.02,"project2501.ca", fontsize=8)
    plt.savefig('phase_guide.pdf', orientation='landscape',bbox_inches='tight')
