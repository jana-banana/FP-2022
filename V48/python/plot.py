import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

import os

if os.path.exists("../build") == False:
    os.mkdir("../build")



t_1, T_1, I_1 = np.genfromtxt('data/m1.txt', unpack=True)

T_1 += 273.15 # in kelvin
I_1 *= 10 #pico ampere

t_2, T_2, I_2 = np.genfromtxt('data/m2.txt', unpack=True)

T_2 += 273.15 # in kelvin
I_2 *= 10 #pico ampere


plt.figure()
plt.plot(T_1, I_1, '.k', label='Messung 1')
plt.plot(T_2, I_2, '.r', label='Messung 2')

# in matplotlibrc leider (noch) nicht möglich
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot.pdf')
