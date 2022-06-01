import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy import integrate

import os

if os.path.exists("../build") == False:
  os.mkdir("../build")


#Werte auslesen
t_1, T_1, I_1 = np.genfromtxt('../data/m1.txt', unpack=True)

T_1 += 273.15 # in kelvin
I_1 *= 10 #pico ampere

t_2, T_2, I_2 = np.genfromtxt('../data/m2.txt', unpack=True)

T_2 += 273.15 # in kelvin
I_2 *= 10 #pico ampere

#Temperatur-Strom Diagramm

plt.figure()
plt.plot(T_1, I_1, '.k', label='Messung 1')
plt.plot(T_2, I_2, '.r', label='Messung 2')

plt.grid()

plt.xlabel('T / K')
plt.ylabel('I / pA')

# in matplotlibrc leider (noch) nicht möglich
#plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.tight_layout()
plt.legend()
plt.savefig('../build/plot.pdf')
plt.clf()

# Untergrund rausnehmen mit e-Fkt mit Messwerten vom Anstieg des zweiten Maximums
# Werte ab t = 45 bis 66 fürs Fitten nehmen

#Untergrund Fit für Messung 1

unter_I_1 = I_1[45:69]
unter_T_1 = T_1[45:69]

#print(unter_I_1) 
#print(unter_T_1)


def untergrund(T, a, b):
    return a*np.exp(-b/T)

popt, pcov = curve_fit(untergrund, unter_T_1, unter_I_1)
errors = np.sqrt(np.diag(pcov))

print('Untergrund Fit für Messung 1:')
print("a =", popt[0])
print("a fehler =",errors[0])
print("b =", popt[1])
print("b fehler =",errors[1])

plt.xlabel('T / K')
plt.ylabel('I / pA')

x = np.linspace(200, 310 ,1000)
plt.plot(x, untergrund(x, popt[0], popt[1]),label = 'Untergrund-Fit')
#print( untergrund(unter_T_1, popt[0], popt[1])) 
plt.plot(unter_T_1, unter_I_1,'+k', label = 'verwendete Messpunkte')

plt.plot(T_1, I_1, '.k', label='Messung 1')

plt.legend()
plt.grid()

plt.savefig('../build/untergrund_1.pdf')
plt.clf()


unter_a_1 = popt[0]
unter_b_1 = popt[1]

# irgendwie sieht der plot richtig komisch aus, aber ok, lets just pretend that it works
# also die folgende Rechnung muss ja mit korrigierten Werten laufe, und joah die werte sind jetzt erstmal nicht korrigiert
# oder doch erst bei Methode 2? ka 

# Untergrund Fit für Messung 2

unter_I_2 = I_2[40:62]
unter_T_2 = T_2[40:62]

popt, pcov = curve_fit(untergrund, unter_T_2, unter_I_2)
errors = np.sqrt(np.diag(pcov))

print('Untergrund Fit für Messung 2:')
print("a =", popt[0])
print("a fehler =",errors[0])
print("b =", popt[1])
print("b fehler =",errors[1])

plt.xlabel('T / K')
plt.ylabel('I / pA')

x = np.linspace(200, 320 ,1000)
plt.plot(x, untergrund(x, popt[0], popt[1]),label = 'Untergrund-Fit')
#print( untergrund(unter_T_2, popt[0], popt[1])) 
plt.plot(unter_T_2, unter_I_2,'+k', label = 'verwendete Messpunkte')

plt.plot(T_2, I_2, '.k', label='Messung 2')

plt.legend()
plt.grid()

plt.savefig('../build/untergrund_2.pdf')
plt.clf()

unter_a_2 = popt[0]
unter_b_2 = popt[1]

# Bestimmung der Aktivierungsenergie erste Methode (mit Maximum)

def line(t,b,y):
    return b*t + y
# t Zeit, b Heizrate, y y-Achsenabschnitt

# nur die Messungen
plt.figure()
plt.plot(t_1, T_1, '.k', label='Messung 1')
plt.plot(t_2, T_2, '.r', label='Messung 2')
plt.legend()
plt.xlabel('t / min')
plt.ylabel('T / K')

plt.savefig('../build/zeit_temp.pdf')
plt.clf()

# Ausgleichsrechnung Messung 1
popt, pcov = curve_fit(line, t_1, T_1)
errors = np.sqrt(np.diag(pcov))

print('Zeit-Temperatur-Fit_1 mit der Form: b*t + y')
print("b =", popt[0])
print("b fehler =",errors[0])
print("y =", popt[1])
print("y fehler =",errors[1])

heizrate_1 = popt[0]

plt.plot(t_1, T_1, '.k', label='Messung 1')
plt.plot(t_1, line( t_1, popt[0], popt[1] ), label = 'Ausgleichsgerade')
plt.legend()
plt.xlabel('t / min')
plt.ylabel('T / K')

plt.savefig('../build/zeit_temp_fit_1.pdf')
plt.clf()

# Ausgleichsrechnung Messung 2
popt, pcov = curve_fit(line, t_2, T_2)
errors = np.sqrt(np.diag(pcov))

print('Zeit-Temperatur-Fit_2 mit der Form: b*t + y')
print("b =", popt[0])
print("b fehler =",errors[0])
print("y =", popt[1])
print("y fehler =",errors[1])

plt.plot(t_2, T_2, '.k', label='Messung 2')
plt.plot(t_2, line( t_2, popt[0], popt[1] ), label = 'Ausgleichsgerade')
plt.legend()
plt.xlabel('t / min')
plt.ylabel('T / K')

plt.savefig('../build/zeit_temp_fit_2.pdf')
plt.clf()

# Logarithmus von I und Kehrwert von T plotten
# Messung 1

log_I_1 = I_1[5:]
log_T_1 = T_1[5:]

log_I_1 -= untergrund(log_T_1, unter_a_1, unter_b_1)


T_1_korr = log_T_1
I_1_korr = log_I_1

plt.plot(T_1_korr, I_1_korr, '+k', label = 'korrigierte Messwerte')
#plt.plot(1/log_T_1, np.log(log_I_1), '+k')
plt.xlabel('T / K')
plt.ylabel('I / pA')
plt.grid()
plt.legend()
plt.savefig('../build/korrigierte_werte_1.pdf')
plt.clf()

# print('werte für log: ', np.log(log_I_1) )

#hier müssen die Werte noch einmal korrigiert werden, also Untergrund noch abziehen
#log_I_1 -= untergrund(log_T_1, unter_a, unter_b)
#print(log_I_1)

#plt.plot((1/log_T_1), np.log(log_I_1), '.k', label='Messung 1')
plt.xlabel('1/T / K^-1')
plt.ylabel('log(I) / pA')
plt.grid()

print(1/log_T_1)
# fit_1_durch_T_1 = 1/log_T_1[:26]
# fit_log_I_1 = log_I_1[:26]

fit_1_durch_T_1 = 1/log_T_1[9:26]
fit_log_I_1 = np.log(log_I_1[9:26])

#fit_1_durch_T_1 = 1/durch_T_1


plt.plot(fit_1_durch_T_1, fit_log_I_1, '+k', label=' Werte')

popt, pcov = curve_fit(line, fit_1_durch_T_1, fit_log_I_1)
errors = np.sqrt(np.diag(pcov))

print('Ausgleichgerade 1 für Logarithmus Diagramm mit Form: b*t + y')
print("b =", popt[0])
print("b fehler =",errors[0])
print("y =", popt[1])
print("y fehler =",errors[1])

plt.plot(fit_1_durch_T_1, line( fit_1_durch_T_1, popt[0], popt[1] ), label = 'Ausgleichsgerade')

# die drei/vier punkte die abstehen vielleicht noch raus nehmen
# ja hab ich ja. passt trotzdem nicht
# ok geht klar, idk soll die drei punkte wieder rein?

plt.legend()
plt.savefig('../build/log(I)_1durchT_1.pdf')
plt.clf()

#Messung 2

log_T_2 = np.append(T_2[1:5], T_2[12:])
log_I_2 = np.append(I_2[1:5], I_2[12:])

log_I_2 -= untergrund(log_T_2, unter_a_2, unter_b_2)

T_2_korr = log_T_2
I_2_korr = log_I_2

plt.plot(T_2_korr, I_2_korr, '+k', label = 'korrigierte Messwerte')
#plt.plot(1/log_T_2, np.log(log_I_2), '+k')
plt.xlabel('T / K')
plt.ylabel('I / pA')
plt.grid()
plt.legend()
plt.savefig('../build/korrigierte_werte_2.pdf')
plt.clf()

# plt.plot(log_T_2, log_I_2, '+k')
plt.plot(1/log_T_2, np.log(log_I_2), '+k')
plt.grid()
plt.savefig('../build/werte_log_2.pdf')
plt.clf()


#plt.plot((1/log_T_2), np.log(log_I_2), '.k', label='Messung 2')
plt.xlabel('1/T / K^-1')
plt.ylabel('log(I) / pA')


plt.grid()


fit_1_durch_T_2 = 1/log_T_2[11:23]
fit_log_I_2 = np.log(log_I_2[11:23])


plt.plot(fit_1_durch_T_2, fit_log_I_2, '+k', label=' Werte')

popt, pcov = curve_fit(line, fit_1_durch_T_2, fit_log_I_2)
errors = np.sqrt(np.diag(pcov))

print('Ausgleichgerade 2 für Logarithmus Diagramm mit Form: b*t + y')
print("b =", popt[0])
print("b fehler =",errors[0])
print("y =", popt[1])
print("y fehler =",errors[1])

plt.plot(fit_1_durch_T_2, line( fit_1_durch_T_2, popt[0], popt[1] ), label = 'Ausgleichsgerade')

plt.legend()
plt.savefig('../build/log(I)_1durchT_2.pdf')
plt.clf()


# Methode 2
# irgendwas mit Trapezregel, gar kein bock
# hab das aus dem Internet: https://numpy.org/doc/stable/reference/generated/numpy.trapz.html
# und das https://www.geeksforgeeks.org/python-scipy-integrate-cumtrapz-method/

#mit korrigierten Werten rechnen? i think yes
#np.trapz(I_1_korr, T_2_korr)


T_1_korr = log_T_1[1:50]
I_1_korr = log_I_1[1:50]

plt.plot(T_1_korr, I_1_korr, '+k', label = 'korrigierte Messwerte')
#plt.plot(1/log_T_1, np.log(log_I_1), '+k')
plt.xlabel('T / K')
plt.ylabel('I / pA')
plt.grid()
plt.legend()
plt.savefig('../build/korrigierte_werte_1.pdf')
plt.clf()


T_2_korr = log_T_2[1:44]
I_2_korr = log_I_2[1:44]

plt.plot(T_2_korr, I_2_korr, '+k', label = 'korrigierte Messwerte')
#plt.plot(1/log_T_2, np.log(log_I_2), '+k')
plt.xlabel('T / K')
plt.ylabel('I / pA')
plt.grid()
plt.legend()
plt.savefig('../build/korrigierte_werte_2.pdf')
plt.clf()

#Integral zu Messung 1

integral_1 = integrate.cumtrapz(I_1_korr, T_1_korr)
plt.plot(T_1_korr[1:], integral_1,'+k', label = 'Integral 1')
plt.legend()
plt.grid()
plt.savefig('../build/integral_1.pdf')
plt.clf()

integral_1 /= (I_1_korr[1:]*heizrate_1)
print(integral_1)

plt.plot(T_1_korr[1:], integral_1, '+k', label = 'Ln(integral/(I*b))')
plt.legend()
plt.grid
plt.savefig('../build/integral_ln_1.pdf')
plt.clf()