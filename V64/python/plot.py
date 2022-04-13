import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

import os

if os.path.exists("../build") == False:
    os.mkdir("../build")


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

delta   = 10                #Platten sind jeweils um 10° geneigt
theta   = 10                #Platten wurden von 0 bis -10° verschoben
D       = 1                 #Dicke der Platten in mm
lamda   = 632.990 *10**(-6) #Wellenlänge Laser in mm

n_Glas = (2*delta*theta*D)/(2*delta*theta*D - M*lamda)
print("Brechungsindex Glas: ", n_Glas)
print("Mittelwert n_Glas: ", np.mean(n_Glas))


#Brechungsindex Luft
p, M_1, M_2, M_3, M_4 = np.genfromtxt('data/Brechungsindex_Luft_Messung.txt', unpack=True)
#M_1 und M_2 ohne Haube, M_3 und M4 mit Haube

#L = ufloat(100,1) #Länge Gaskammer in mm
L = 100

n_Luft_1 = (M_1*lamda)/L +1
print("n_Luft erster Durchgang ohne Haube: ", n_Luft_1)
print("Mittelwert erster Durchgang ohne Haube: ", np.mean(n_Luft_1))

n_Luft_2 = (M_2*lamda)/L +1
print("n_Luft zweiter Durchgang ohne Haube: ", n_Luft_2)
print("Mittelwert zweiter Durchgang ohne Haube: ", np.mean(n_Luft_2))

n_Luft_ohne = np.concatenate([n_Luft_1,n_Luft_2])
print("Mittelwert ohne Haube: ", np.mean(n_Luft_ohne))


n_Luft_3 = (M_3*lamda)/L +1
print("n_Luft dritter Durchgang mit Haube: ", n_Luft_3)
print("Mittelwert dritter Durchgang mit Haube: ", np.mean(n_Luft_3))
print("Abweichung ohne Haube: ", np.mean(n_Luft_ohne)/1.000292)

n_Luft_4 = (M_4*lamda)/L +1
print("n_Luft vierter Durchgang mit Haube: ", n_Luft_4)
print("Mittelwert vierter Durchgang mit Haube: ", np.mean(n_Luft_4))

n_Luft_mit = np.concatenate([n_Luft_3,n_Luft_4])
print("Mittelwert mit Haube: ", np.mean(n_Luft_mit))
print("Abweichung mit Haube: ", np.mean(n_Luft_mit)/1.000292)

#Lorentz-Lorenz Gesetz

T = 21.1 + 273.15 #Temperatur an dem Tag in Kelvin
R = 8.31446261815324 #universelle Gaskonstante
# p Druck in Gaskammer
# a müsste 3A/2 sein, wobei A diese molare refraktivität ist


def line(p, a, b):
    return a*p /(T*R) + b 

#erster Durchgang ohne Haube
popt, pcov = curve_fit(line, p, n_Luft_1)
errors = np.sqrt(np.diag(pcov))
print("a1 =", popt[0])
print("a1 fehler =",errors[0])
print("b1 =", popt[1])
print("b1 fehler =",errors[1])

a1 = ufloat(popt[0],errors[0])
b1 = ufloat(popt[1],errors[1])

plt.plot(p, n_Luft_1, 'r+', label = 'Messwerte 1')
plt.plot(p, line(p, popt[0], popt[1]),'r-', label='Fit 1')

#zweiter Durchgang ohne Haube
popt, pcov = curve_fit(line, p, n_Luft_2)
errors = np.sqrt(np.diag(pcov))
print("a2 =", popt[0])
print("a2 fehler =",errors[0])
print("b2 =", popt[1])
print("b2 fehler =",errors[1])

a2 = ufloat(popt[0],errors[0])
b2 = ufloat(popt[1],errors[1])

plt.plot(p, n_Luft_2, 'b+', label = 'Messwerte 2')
plt.plot(p, line(p, popt[0], popt[1]),'b-', label='Fit 2')

#dritter Durchgang mit Haube
popt, pcov = curve_fit(line, p, n_Luft_3)
errors = np.sqrt(np.diag(pcov))
print("a3 =", popt[0])
print("a3 fehler =",errors[0])
print("b3 =", popt[1])
print("b3 fehler =",errors[1])

a3 = ufloat(popt[0],errors[0])
b3 = ufloat(popt[1],errors[1])

plt.plot(p, n_Luft_3, 'm+', label = 'Messwerte 3')
plt.plot(p, line(p, popt[0], popt[1]),'m-', label='Fit 3')

#vierter Durchgang mit Haube
popt, pcov = curve_fit(line, p, n_Luft_4)
errors = np.sqrt(np.diag(pcov))
print("a4 =", popt[0])
print("a4 fehler =",errors[0])
print("b4 =", popt[1])
print("b4 fehler =",errors[1])

a4 = ufloat(popt[0],errors[0])
b4 = ufloat(popt[1],errors[1])


plt.plot(p, n_Luft_4, 'k+', label = 'Messwerte 4')
plt.plot(p, line(p, popt[0], popt[1]),'k-', label='Fit 4')

plt.xlabel(r"$p \, / \,$mbar")
plt.ylabel('n')
plt.grid()
plt.tight_layout()
plt.legend()
plt.savefig('build/lorentz_lorenz.pdf')

plt.clf()

a = np.array([a1, a2, a3, a4])
b = np.array([b1, b2, b3, b4])
a_mittel = np.mean(a)
b_mittel = np.mean(b)
print("a ",a_mittel)
print("b ",b_mittel)

print("Aufgabe aus Anleitung: ", a_mittel*1.013 /(288.15*R) + b_mittel )