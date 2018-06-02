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

# Return the temperatures for an array of densities (in H/cc) at constant Jeans mass (in Msol)
def temp_jeans_mass(rho, mass):
    return np.power(np.square(mass/22.1)*rho, 1./3)

# Return the temperatures for an array of densities (in H/cc) at constant Jeans length (in pc)
def temp_jeans_length(rho, length):
    return np.square(length/9.628)*rho

# Return the temperatures for an array of densities (in H/cc) at a constant entropy (in keV cm^2)
def temp_adiabat(rho, K):
    return K/(8.617e-8*np.power(rho, -2./3))

# Return the temperatures for an array of densities at a constant pressure (in K/cc/kB)
def temp_pressure(rho, P):
    return P/rho

# Return the temperature for an array of densities given a constant courant time (in years) NOTE: 
# There is a missing prefactor here that is code dependent, n^-1/3 is NOT exactly delta-x
def temp_courant_time(rho, tc):
    return np.power(rho, -2./3)*1.67e-27/1.38e-23*3./5./np.square(tc/3.154e7)

# Return the temperature for an array of sound speeds (in km/s)
def temp_sound_speed(c_s):
    return np.square(1e3*c_s)*(3./5*1.67e-27/1.38e-23)

# Return the densities for an array of free-fall times (in Myr)
def rho_freefall(t_ff):
    return np.square(51.5/t_ff)


if __name__ == "__main__":
    font = {'size': 14}
    rc('font', **font)
    rc('axes', grid=False)
    rc('axes', facecolor='white')
    rc('figure', figsize=(11,8.5))
    rho_range = (-7,4)
    temp_range = (2,9)
    rho = np.logspace(*rho_range) # we will use this to calculate isocontours

    # Plot lines of constant cooling time
    coolrates = np.genfromtxt('coolrates.dat', skip_header=2, names=True)
    Edot_per_erg = coolrates['Edot']/(coolrates['nH']/0.74*1.672e-24)
    tcool = -coolrates['E']/(Edot_per_erg*3.155e7)
    tcool[np.where(tcool < 0)] = np.nan #To keep contours from joining up
    plt.contour(coolrates['nH'][::61], coolrates['T'][:61], np.reshape(tcool, (81,61)).T, 
            (1e4, 1e6, 1e8), linewidth=50, colors='Blue')
    plt.text(2e-1, 2e8, "$10^8yr$", color="Blue", rotation='60')
    plt.text(2e1, 2e8, "$10^6yr$", color="Blue", rotation='60')
    plt.text(1e3, 2e8, "$t_{cool}=10^4yr$", color="Blue", rotation='60')

    #Set up axes properties
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(np.power(10., rho_range))
    plt.ylim(np.power(10., temp_range))
    plt.ylabel('T $(K)$')
    plt.xlabel(r'$\rho (m_p /cm^3)$')

    # Plot simulation data
    pdata = np.load('g1536_gas.npz')
    plt.hexbin(pdata['rho'], pdata['temp'], xscale='log', yscale='log', bins='log', extent=(-7,4,2,9),
            mincnt=2, gridsize=200, cmap='binary')

    # Add interesting density values
    plt.axvline(1,color='Salmon',linestyle=':')
    plt.text(1.2,5e8,'Mean MW ISM', color='Salmon', rotation='270')
    plt.axvline(5e-6,color='Salmon',linestyle=':')
    plt.text(6e-6,9e6,r'$\rho_{crit}\; (\Omega = 1)$', color='Salmon', rotation='270')
    plt.axvline(0.308*5e-6,color='Salmon',linestyle=':')
    plt.text(0.308*6e-6,9e6,r'$\Omega_m = 0.308$', color='Salmon', rotation='270')
    plt.axvline(0.045*5e-6,color='Salmon',linestyle=':')
    plt.text(0.045*6e-6,9e6,r'$\Omega_B = 0.045$', color='Salmon', rotation='270')

    # Annotate lines of constant courant time
    plt.plot(rho, temp_courant_time(rho, 1e2), color='cyan', linestyle='-.')
    plt.plot(rho, temp_courant_time(rho, 1e3), color='cyan', linestyle='-.')
    plt.plot(rho, temp_courant_time(rho, 1e4), color='cyan', linestyle='-.')
    annotate_lines({r"$t_{courant} = const$":(3e1,1e6)}, -2./3, 'cyan')

    # Lines of constant pressure
    plt.plot(rho, temp_pressure(rho, 1e2), color='grey', linestyle='--')
    plt.plot(rho, temp_pressure(rho, 1e3), color='grey', linestyle='--')
    plt.plot(rho, temp_pressure(rho, 1e4), color='grey', linestyle='--')
    plt.plot(rho, temp_pressure(rho, 1e5), color='grey', linestyle='--')
    anns = {"(isobar)":(1.1e-7,4e8),
            r"$P/k = 100\; K\; cm^{-3}$":(1.5e-7,8e8),
            r"$P/k = 10^3\; K\; cm^{-3}$":(1.5e-6,8e8),
            r"$P/k = 10^4\; K\; cm^{-3}$":(1.5e-5,8e8),
            r"$P/k = 10^5\; K\; cm^{-3}$":(1.5e-4,8e8)}
    annotate_lines(anns, -1, 'grey')

    # Lines of constant entropy
    plt.plot(rho, temp_adiabat(rho, 0.1), color='k', linestyle=':')
    plt.plot(rho, temp_adiabat(rho, 1), color='k', linestyle=':')
    plt.plot(rho, temp_adiabat(rho, 10), color='k', linestyle=':')
    plt.plot(rho, temp_adiabat(rho, 100), color='k', linestyle=':')
    anns = {r"$K = 100\; eV\; cm^2$ (adiabat)":(6e-7, 1e2),
            r"$K = 1\; keV\; cm^2$":(2e-7, 5e2),
            r"$K = 10\; keV\; cm^2$":(2e-7, 5e3),
            r"$K = 100\; keV\; cm^2$":(2e-7, 5e4)}
    annotate_lines(anns, 2./3, 'k')

    # Lines of constant Jeans Mass
    plt.plot(rho, temp_jeans_mass(rho, 1e3), color='m', linestyle='-')
    plt.plot(rho, temp_jeans_mass(rho, 1e4), color='m', linestyle='-')
    plt.plot(rho, temp_jeans_mass(rho, 1e5), color='m', linestyle='-')
    plt.plot(rho, temp_jeans_mass(rho, 1e6), color='m', linestyle='-')
    anns = {r"$10^3 M_\odot$":(1.5e3,1.8e2), 
            r"$10^4 M_\odot$":(1.5e3,8e2), 
            r"$M_J = 10^5 M_\odot$":(3e2,2.5e3), 
            r"$M_J = 10^6 M_\odot$":(3e2,1.2e4)}
    annotate_lines(anns, 1./3, 'm')

    # Lines of constant Jeans Length
    plt.plot(rho, temp_jeans_length(rho, 1e1), color='green', linestyle='-')
    plt.plot(rho, temp_jeans_length(rho, 1e2), color='green', linestyle='-')
    plt.plot(rho, temp_jeans_length(rho, 1e3), color='green', linestyle='-')
    plt.plot(rho, temp_jeans_length(rho, 1e4), color='green', linestyle='-')
    anns = {r"$\lambda_J = 10 pc$":(7e1,1e2),
            r"$\lambda_J = 100 pc$":(7e-1,1e2),
            r"$\lambda_J = 10^3 pc$":(7e-3,1e2),
            r"$\lambda_J = 10^4 pc$":(7e-5,1e2)}
    annotate_lines(anns, 1, 'g')

    # Show sound speed as second y axis
    psound = plt.twinx()
    psound.set_ylabel('$c_s (km/s)$')
    temp_dyn_range = temp_range[1] - temp_range[0]
    temp_cs_vals = (np.log10(temp_sound_speed(np.array([3e0,1e1,3e1,1e2,3e2,1e3,3e3])))-temp_range[0])/temp_dyn_range
    psound.set_yticks(temp_cs_vals)
    psound.set_yticklabels(['3','10', '30', '100', '300', '1,000', '3,000'])

    # Show free-fall time as second x axis
    pfall = plt.twiny()
    pfall.set_xlabel('$t_{ff} (Myr)$')
    rho_dyn_range = rho_range[1] - rho_range[0]
    rho_ff_vals = (np.log10(rho_freefall(np.array([1e0, 1e1, 1e2, 1e3, 1e4, 1e5])))-rho_range[0])/rho_dyn_range
    pfall.set_xticks(rho_ff_vals)
    pfall.set_xticklabels([r'$1$', r'$10$', r'$100$', r'$10^3$', r'$10^4$', r'$10^5$'])
    plt.figtext(0.71,0.04,"BW Keller's Reference Phase Diagram", fontsize=8)
    plt.figtext(0.825,0.02,"project2501.ca", fontsize=8)
    plt.savefig('phase_guide.pdf', orientation='landscape',bbox_inches='tight')
