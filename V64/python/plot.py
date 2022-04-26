import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

import os

if os.path.exists("build") == False:
    os.mkdir("build")


#Kontrast

phi, U_min, U_max = np.genfromtxt('data/kontrastmessung.txt', unpack=True)

K = (U_max - U_min)/(U_max + U_min)
print("Kontrast: ", K)

def line(Phi, A):
    return A*(np.abs(np.cos(Phi)*np.sin(Phi)))

phi *= (2*np.pi)/360

popt, pcov = curve_fit(line, phi, K)
errors = np.sqrt(np.diag(pcov))
print("A =", popt[0])
print("A fehler =",errors[0])

plt.plot(phi, K, '.', label = 'Messwerte')
plt.xlabel(r"$\theta \, / \, \mathrm{rad}$")
plt.ylabel('K')

plt.ylim(0,1.2)
pi= np.pi
lol = np.arange(-2 * pi, 2 * pi+pi/2, step=(pi / 2))
plt.xticks(lol, ['-2π', '-3π/2', 'π', '-π/2', '0', 'π/2', 'π', '3π/2', '2π'])
x = np.linspace(0, np.pi, 1000)
plt.plot(x, line(x, popt[0]),'r-', label='Ausgleichsfunktion')
plt.grid()
plt.tight_layout()
plt.legend()
plt.savefig('build/kontrast_ausgleich.pdf')

plt.clf()

#Brechungsindex Glas
i, M = np.genfromtxt('data/Brechungsindex_Glas_Messung.txt', unpack=True)

#korrektur: winkel in rad, nicht grad

delta   = 10*(2*pi/360)              #Platten sind jeweils um 10° geneigt        
theta   = 10*(2*pi/360)              #Platten wurden von 0 bis -10° verschoben
D       = 1                 #Dicke der Platten in mm
lamda   = 632.990 *10**(-6) #Wellenlänge Laser in mm

n_Glas = (2*delta*theta*D)/(2*delta*theta*D - M*lamda)
print("Brechungsindex Glas: ", n_Glas)
print("Mittelwert n_Glas: ", np.mean(n_Glas))
abw = (np.mean(n_Glas) - 1.45)/1.45
print("Abweichung", abw)

#Brechungsindex Luft -- Korrektur

p, M_1, M_2, M_3, M_4 = np.genfromtxt('data/Brechungsindex_Luft_Messung.txt', unpack=True)
#M_1 und M_2 ohne Haube, M_3 und M4 mit Haube

L = 100 

n_Luft_1 = (M_1*lamda)/L +1
print("n_Luft erster Durchgang ohne Haube: ")
for i in range(20):
    print(f'{n_Luft_1[i]:.8f}')

n_Luft_2 = (M_2*lamda)/L +1
print("n_Luft zweiter Durchgang ohne Haube: ")
for i in range(20):
    print(f'{n_Luft_2[i]:.8f}')

n_Luft_3 = (M_3*lamda)/L +1
print("n_Luft dritter Durchgang mit Haube: ")
for i in range(20):
    print(f'{n_Luft_3[i]:.8f}')

n_Luft_4 = (M_4*lamda)/L +1
print("n_Luft vierter Durchgang mit Haube: ")
for i in range(20):
    print(f"{n_Luft_4[i]:.8f}")


# Lorentz - Lorenz Gesetz

T = 294.25 #Temperatur an dem Tag in Kelvin
R = const.R #universelle Gaskonstante
p_plot = np.linspace(0, 1000, 10000)

def line(p, a, b):
    return a*p /(T*R) + b 

params_1, covariance_matrix_1 = curve_fit(line, p, n_Luft_1)
# params_1, covariance_matrix_1 = np.polyfit(p, n_Luft_1, deg=1, cov=True)


errors_1 = np.sqrt(np.diag(covariance_matrix_1))

for name, value, error in zip('ab', params_1, errors_1):
    print(f'{name} = {value:.8f} ± {error:.8f}')

plt.plot(p, n_Luft_1, 'r+', label = 'Messwerte 1')
plt.plot(p_plot, line(p_plot,params_1[0], params_1[1]),'r-', label='Fit 1')
# plt.plot(p_plot, p_plot*params_1[0] + params_1[1],'r-', label='Fit 1')



params_2, covariance_matrix_2 = curve_fit(line, p, n_Luft_2)
# params_2, covariance_matrix_2 = np.polyfit(p, n_Luft_2, deg=1, cov=True)


errors_2 = np.sqrt(np.diag(covariance_matrix_2))

for name, value, error in zip('ab', params_2, errors_2):
    print(f'{name} = {value:.8f} ± {error:.8f}')

plt.plot(p, n_Luft_2, 'b+', label = 'Messwerte 2')
plt.plot(p_plot, line(p_plot, params_2[0], params_2[1]),'b-', label='Fit 2')
# plt.plot(p_plot, p_plot*params_2[0] + params_2[1],'b-', label='Fit 2')



params_3, covariance_matrix_3 = curve_fit(line, p, n_Luft_3)
# params_3, covariance_matrix_3 = np.polyfit(p, n_Luft_3, deg=1, cov=True)


errors_3 = np.sqrt(np.diag(covariance_matrix_3))

for name, value, error in zip('ab', params_3, errors_3):
    print(f'{name} = {value:.8f} ± {error:.8f}')

plt.plot(p, n_Luft_3, 'm+', label = 'Messwerte 3')
plt.plot(p_plot, line(p_plot, params_3[0], params_3[1]),'m-', label='Fit 3')
# plt.plot(p_plot, p_plot*params_3[0] + params_3[1],'m-', label='Fit 3')



params_4, covariance_matrix_4 = curve_fit(line, p, n_Luft_4)
# params_4, covariance_matrix_4 = np.polyfit(p, n_Luft_4, deg=1, cov=True)

errors_4 = np.sqrt(np.diag(covariance_matrix_4))

for name, value, error in zip('ab', params_4, errors_4):
    print(f'{name} = {value:.8f} ± {error:.8f}')

plt.plot(p, n_Luft_4, 'k+', label = 'Messwerte 4')
plt.plot(p_plot, line(p_plot, params_4[0], params_4[1]),'k-', label='Fit 4')
# plt.plot(p_plot,p_plot*params_4[0] + params_4[1],'k-', label='Fit 4')


plt.xlabel(r"$p \, / \,$mbar")
plt.ylabel('n')
plt.grid()
plt.tight_layout()
plt.legend()
plt.savefig('build/lorentz_lorenz.pdf')



a = unp.uarray([params_1[0], params_2[0], params_3[0], params_4[0]], [errors_1[0], errors_2[0], errors_3[0], errors_4[0]])
b = unp.uarray([params_1[1], params_2[1], params_3[1], params_4[1]], [errors_1[1], errors_2[1], errors_3[1], errors_4[1]])

n_atm = a*1013/(288.15*R) + b
print('\n')
for i in range(4):
    print(f'{n_atm[i]:.8f}')

n_ohne = (a[0] + a[1])/2*1013/(288.15*R) + (b[0] + b[1])/2
print("n_ohne", n_ohne)
print("abw", (331 - 292)/292)

n_mit = (a[2] + a[3])/2*1013/(288.15*R) + (b[2] + b[3])/2
print("n_mit", n_mit)
print("abw", (310 - 292)/292)

a_mittel = sum(a)/len(a)
b_mittel = sum(b)/len(b)
print("a ",a_mittel)
print("b ",b_mittel)

print("Aufgabe aus Anleitung: ", a_mittel*1013/(288.15*R) + b_mittel )
print("abw", (320 - 292)/292)
# print("Aufgabe aus Anleitung: ", a_mittel*1013*T/288.15 + b_mittel )
