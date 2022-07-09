import matplotlib.pyplot as plt
# from matplotlib.axes.Axes import axhspan
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
import uncertainties as unc

import os

if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("data")==False:
    os.mkdir("data")
#------------------------------------------------------------------------------Allgemeines---------------------------------------------------

def logausdruck(p, s_p, E, O):
    p = unp.uarray(p, s_p)
    z = p - E 
    n = O - E 
    return unp.log(z/n)

def line(x, m, n):
    return m*x + n

def fehler(exp, theo):
    dif = theo - exp
    return dif / theo

def middel(a,b,c):
    arr = unp.uarray(np.zeros(np.size(a)),np.zeros(np.size(a)))
    for i in range(0,np.size(a)):
        bar = np.mean([a[i], b[i], c[i]])
        arr[i] = unc.ufloat(bar ,( (a[i] - bar)**2 + (b[i] - bar)**2 + (c[i] - bar)**2)/ np.sqrt(3*2))
    return arr

def middeldeux(a,b):
    arr = unp.uarray(np.zeros(np.size(a)),np.zeros(np.size(a)))
    for i in range(0,np.size(a)):
        bar = np.mean([a[i], b[i]])
        arr[i] = unc.ufloat(bar ,( (a[i] - bar)**2 + (b[i] - bar)**2 )/ np.sqrt(1*2))
    return arr

def evaku(t, la,g0, g1, g2, name):
    ln_line_1 = unp.nominal_values(la)[:g1]
    ln_line_2 = unp.nominal_values(la)[g1:g2]
    ln_line_3 = unp.nominal_values(la)[g2:]

    t_line_1 = t[:g1]
    t_line_2 = t[g1:g2]
    t_line_3 = t[g2:]

    ln_line = [ln_line_1, ln_line_2, ln_line_3]
    t_line = [t_line_1, t_line_2, t_line_3]
    nums = range(3)
    m = np.zeros(3)
    sig_m = np.zeros(3)
    n = np.zeros(3)
    sig_n = np.zeros(3)

    for i, secs, tims in zip(nums, ln_line, t_line):

        popt, pcov = curve_fit(line, tims, secs)
        errors = np.sqrt(np.diag(pcov))

        print('Messung', name, f' mit Abschnitt {i+1}')
        print(f'Steigung: \n m = {popt[0]} \pm {errors[0]} \n')
        m[i] = popt[0]
        sig_m[i] = errors[0]
        print(f'Achsenabschnitt: \n n = {popt[1]} \pm {errors[1]} \n')
        n[i] = popt[1]
        sig_n[i] = errors[1]
    
    t_plot1 = np.linspace(-10 + g0, g1*10 + 50, 1000)
    t_plot2 = np.linspace(g1*10 - 50, g2*10 +50, 1000)
    t_plot3 = np.linspace(g2*10-50, 610, 1000)
    if len(t)<27:
        t_plot1 = np.linspace(-10 + g0, g1*5 + 10, 1000)
        t_plot2 = np.linspace(g1*5 - 10, g2*10 +10, 1000)
        t_plot3 = np.linspace(g2-10, 130, 1000)

    plt.figure()
    # plt.errorbar(t, unp.nominal_values(la), yerr=unp.std_devs(la), fmt='.k', label='Nicht benutzte Messdaten')
    plt.errorbar(t_line_1, ln_line_1, yerr=unp.std_devs(la)[:g1], fmt='.', color='red', label='Daten des 1. linearen Bereiches')
    plt.plot(t_plot1, line(t_plot1, m[0], n[0]), '-b', label='1. Ausgleichsrechnung')
    plt.errorbar(t_line_2, ln_line_2, yerr=unp.std_devs(la)[g1:g2], fmt='.', color='orange', label='Daten des 2. linearen Bereiches')
    plt.plot(t_plot2, line(t_plot2, m[1], n[1]), '-g', label='2. Ausgleichsrechnung')
    plt.errorbar(t_line_3, ln_line_3, yerr=unp.std_devs(la)[g2:], fmt='.', color='yellow', label='Daten des 3. linearen Bereiches')
    plt.plot(t_plot3, line(t_plot3, m[2], n[2]), '-m', label='3.Ausgleichsrechung')
    # plt.rc('axes', labelsize=size_label)
    plt.xlabel(r'$t \mathbin{/} \si{\second}$')
    plt.ylabel(r'$\ln(\frac{p-p_\text{E}}{p_0-p_\text{E}})$')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('build/evaku'+name+'.pdf')

    print("--------Saugvermögen--", name, "------- \n ")
    V = 34 #liter
    sig_V = 3.4
    if len(t)<20:
        V = 33 #liter
        sig_V =3.3

    S = - m*V
    sig_S = np.sqrt(V**2*sig_m**2 + m**2*sig_V**2)

    
    if len(t)<20:
        theo = 77 #liters per second
        abw = fehler(S, theo)
        for i in nums:
            print(f'Bereich{i+1} S= {S[i]} \pm {sig_S[i]}, abweichung {abw[i]}')
    else :
        abw_lower = fehler(S, 4.6*10/36)
        abw_higher = fehler(S, 5.5*10/36)
        for i in nums:
            print(f'Bereich{i+1} S= {S[i]} \pm {sig_S[i]}, abweichung {abw_lower[i]}, {abw_higher[i]}')
    print('\n')
    
    return unp.uarray(S, sig_S)

def leckrate(t, p, sig_p, p_G, go, drueck, name):
    t_plot = np.linspace(-10, 210, 1000)
    if name =='turbo':
        t_plot = np.linspace(-5, 125, 1000)

    if drueck=='50mbar' or drueck=='100mbar':
        popt, pcov = curve_fit(line, t[:go], p[:go])
    else:
        popt, pcov = curve_fit(line, t, p)        
    errors = np.sqrt(np.diag(pcov))

    print('Gleichgewichtsdruck',drueck, name)
    print(f'Steigung: \n m = {popt[0]} \pm {errors[0]} \n')
    m = ufloat(popt[0], errors[0])
    print(f'Achsenabschnitt: \n n = {popt[1]} \pm {errors[1]} \n')
    n = ufloat(popt[1], errors[1])

    if drueck=='50mbar' or drueck=='100mbar':
        plt.figure()
        plt.errorbar(t[:go], p[:go], yerr=sig_p[:go], fmt='.k', label='Für die Ausgleichsrechnung benutzten Messwerte')
        plt.errorbar(t[go:], p[go:], yerr=sig_p[go:], fmt='.b', label='Restlichen Messwerte')
        plt.plot(t_plot, line(t_plot, m.n, n.n ), '-r', label='Ausgleichsgerade')
        plt.rc('axes', labelsize=size_label)
        plt.xlabel(r'$t \mathbin{/} \si{\second}$')
        plt.ylabel(r'$p \mathbin{/} \si{\milli\bar}$')
        plt.legend()
        plt.tight_layout()
        plt.savefig('build/leck_'+name+'_'+drueck+'.pdf')
    
    elif drueck=='10mbar' or drueck=='05mbar':
        plt.figure()
        plt.errorbar(t, p, yerr=sig_p, fmt='.k', label='Messwerte')
        plt.plot(t_plot, line(t_plot, m.n, n.n ), '-r', label='Ausgleichsgerade')
        plt.rc('axes', labelsize=size_label)
        plt.xlabel(r'$t \mathbin{/} \si{\second}$')
        plt.ylabel(r'$p \mathbin{/} \si{\milli\bar}$')
        plt.legend()
        plt.tight_layout()
        plt.savefig('build/leck_'+name+'_'+drueck+'.pdf')

    elif name =='turbo':
        plt.figure()
        plt.errorbar(t, p, yerr=sig_p, fmt='.k', label='Messwerte')
        plt.plot(t_plot, line(t_plot, m.n, n.n ), '-r', label='Ausgleichsgerade')
        plt.rc('axes', labelsize=size_label)
        plt.xlabel(r'$t \mathbin{/} \si{\second}$')
        plt.ylabel(r'$p \mathbin{/} \SI{e-6}{\milli\bar}$')
        plt.legend()
        plt.tight_layout()
        plt.savefig('build/leck_'+name+'_'+drueck+'.pdf')

    ###Saugvermögen
    print("--------Saugvermögen--", name,drueck, "------- \n ")

    V = 34 #liter
    sig_V = 3.4
    if name=='turbo':
        V = 33 #liter
        sig_V =3.3

    S = m.n*V/p_G.n
    sig_S = np.sqrt((V/p_G.n)**2*m.s**2 + (m.n/p_G.n)**2*sig_V**2 + (m.n*V/(p_G.n**2))**2*p_G.s**2)

    if name=='turbo':
        theo = 77 #liters per second
        abw = fehler(S, theo)
        print(f'Gleichgewichtsdruck:{p_G.n} S= {S} \pm {sig_S}, abweichung {abw}')
    else :
        abw_lower = fehler(S, 4.6*10/36)
        abw_higher = fehler(S, 5.5*10/36)
        print(f'Gleichgwichtsdruck:{p_G.n} S= {S} \pm {sig_S}, abweichung {abw_lower}, {abw_higher}')
    print('\n')

    return unp.uarray(S, sig_S)

#-------------------------------------------------------------Evakuierungskurven-------------------------------------------------------------
print('___________________________________Evakuierungskurve___________________________________________________________________________ \n')

# firstplot
t_plot = np.linspace(0, 10, 100)
plt.figure()
plt.plot(t_plot, 3*t_plot+4, 'k.', label='tryy')
# plt.rc('axes', labelsize=size_label)
plt.xlabel(r'$t \mathbin{/} \si{\second}$')
plt.ylabel(r'$\ln(\frac{p-p_\text{E}}{p_0-p_\text{E}})$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('build/firstplot.pdf')

#------------------------------------------------------------ Drehschieberpumpe--------------------------------------------------------------
print('_____________________________________Drehschieberpumpe_________________________________________________________________________ \n')

t_eva, p1_eva, p2_eva, p3_eva = np.genfromtxt("data/Dreh_Evakuierung.txt", unpack=True)

numbers = np.arange(len(t_eva))
p_0 = ufloat(1000, 500)
p_E = ufloat(3.85e-3, 1.155e-3)

sig_p1_eva = np.zeros(len(p1_eva))
sig_p2_eva = np.zeros(len(p2_eva))
sig_p3_eva = np.zeros(len(p3_eva))

drucks = [p1_eva, p2_eva, p3_eva]
sigms = [sig_p1_eva, sig_p2_eva, sig_p3_eva]

for p, sig in zip(drucks, sigms):
    for n in range(len(t_eva)):
        if p[n] < 100:
            sig[n] = p[n]*0.3
        elif p[n] >= 100:
            sig[n] = p[n]*0.5
        elif p[n] == 100:
            print("well, hm")
        else:
            print('definitiv ein fehler')

#Messung 3 verschieben
p3_eva = p3_eva[2:]
sig_p3_eva = sig_p3_eva[2:]

#Mittelwert über die Drucke
p_mean = unp.uarray(np.zeros(len(t_eva)),np.zeros(len(t_eva)))
p_mean[:-2] = middel(p1_eva[:-2], p2_eva[:-2], p3_eva)
p_mean[-2:] = middeldeux(p1_eva[-2:], p2_eva[-2:])

sig_p_sys = np.zeros(len(p_mean))
for n in range(0, len(p_mean)):
    if unp.nominal_values(p_mean)[n] < 100:
            sig_p_sys[n] = unp.nominal_values(p_mean)[n]*0.3
    elif unp.nominal_values(p_mean)[n] >= 100:
            sig_p_sys[n] = unp.nominal_values(p_mean)[n]*0.5

#Logarithmusausdruck rechnen
ln = logausdruck(unp.nominal_values(p_mean), sig_p_sys,p_E, p_0)

#Saugvermögen bestimmen
S_eva = evaku(t_eva, ln, 0, 13, 31, 'dreh')

#Tabelle für die Auswertung
#for table
p3_eva_ft = np.zeros(len(t_eva))
sig_p3_eva_ft = np.zeros(len(t_eva))
for i in range(0, len(p3_eva)):
    p3_eva_ft[i] = p3_eva[i]
    sig_p3_eva_ft[i] = sig_p3_eva[i]
np.savetxt('data/tab_dreh_eva.txt', np.column_stack([t_eva, p1_eva, sig_p1_eva, p2_eva, sig_p2_eva, p3_eva_ft, sig_p3_eva_ft, unp.nominal_values(p_mean), sig_p_sys, unp.nominal_values(ln), unp.std_devs(ln)]), fmt=['%.1f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f','%.2f', '%.2f', '%.2f', '%.2f','%.2f'], header='t, m1, m2, m3, mean, ln')

#-------------------------------------------------------------Turbopumpe---------------------------------------------------------------------
print('_____________________________________Turbopumpe________________________________________________________________________________ \n')

print('EINHEIT IST NANO BAR ALSO 10^-9 bar 10^-6 mbar \n')

t_eva_t, p1_eva_t, p2_eva_t, p3_eva_t = np.genfromtxt("data/Turbo_Evakuierung.txt", unpack=True)
p1_eva_t *= 10
p2_eva_t *= 10
p3_eva_t *= 10

numbers = np.arange(len(t_eva_t))
p_0 = ufloat(5000, 1500)
p_E = ufloat(10.9, 3.27)

sig_p1_eva_t = np.zeros(len(p1_eva_t))
sig_p2_eva_t = np.zeros(len(p2_eva_t))
sig_p3_eva_t = np.zeros(len(p3_eva_t))

drucks = [p1_eva_t, p2_eva_t, p3_eva_t]
sigms = [sig_p1_eva_t, sig_p2_eva_t, sig_p3_eva_t]

for p, sig in zip(drucks, sigms):
    for n in range(len(t_eva_t)):
        sig[n] = p[n]*0.3

#Mittelwert berechen
p_mean_t = middel(p1_eva_t, p2_eva_t, p3_eva_t)

sig_p_sys_t = np.zeros(len(p_mean_t))
for n in range(0, len(p_mean_t)):
    sig_p_sys_t[n] = unp.nominal_values(p_mean_t)[n]*0.3

#Logarithmusfunktion anwenden
ln_t = logausdruck(unp.nominal_values(p_mean_t), sig_p_sys_t, p_E, p_0)

#Saugvermögen bestimmen
S_eva_t = evaku(t_eva_t, ln_t, 0, 5, 9, 'turbo')


np.savetxt('data/tab_turbo_eva.txt', np.column_stack([t_eva_t, p1_eva_t, sig_p1_eva_t, p2_eva_t, sig_p2_eva_t, p3_eva_t, sig_p3_eva_t, unp.nominal_values(p_mean_t), sig_p_sys_t, unp.nominal_values(ln_t), unp.std_devs(ln_t)]), fmt=['%.1f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f','%.2f', '%.2f', '%.2f', '%.2f','%.2f'], header='t, m1, m2, m3, mean, ln')

#-------------------------------------------------------------Leckratenmessung---------------------------------------------------------------
print('_____________________________________Leckratenmessung__________________________________________________________________________ \n')

size_label = 19

# firstplot
t_plot = np.linspace(0, 10, 100)
plt.figure()
plt.plot(t_plot, 3*t_plot+4, 'k.', label='tryy')
plt.rc('axes', labelsize=size_label)
plt.xlabel(r'$t \mathbin{/} \si{\second}$')
plt.ylabel(r'$\ln(\frac{p-p_\text{E}}{p_0-p_\text{E}})$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('build/firstplot.pdf')

#-------------------------------------------------------------Drehschieberpumpe--------------------------------------------------------------
print('_____________________________________Drehschieberpumpe_________________________________________________________________________ \n')

names = ['05', '10', '50', '100'] 
numbs = [0, 1, 2, 3]
p_g = [ ufloat(0.5, 0.15), ufloat(10, 3), ufloat(50, 15), ufloat(100, 50)]
gogos = [21, 21, 6, 4]

S = np.zeros(4)
sig_S = np.zeros(4) 

for i, x, pg, go in zip(numbs, names, p_g, gogos):

    t_i, p1_i, p2_i, p3_i = np.genfromtxt('data/Dreh_Leckrate_'+x+'.txt', unpack=True)

    if x=='05':
        p1_i *= 0.001
        p2_i *= 0.001
        p3_i *= 0.001

    sig_p1_i = np.zeros(len(p1_i))
    sig_p2_i = np.zeros(len(p2_i))
    sig_p3_i = np.zeros(len(p3_i))

    drucks = [p1_i, p2_i, p3_i]
    sigms = [sig_p1_i, sig_p2_i, sig_p3_i]

    for p, sig in zip(drucks, sigms):
        for n in range(len(t_i)):
            if p[n] < 100:
                sig[n] = p[n]*0.3
            elif p[n] >= 100:
                sig[n] = p[n]*0.5
            else:
                print('definitiv ein fehler')

    p_i = middel(p1_i, p2_i, p3_i)

    sig_p_i = np.zeros(len(p_i))
    for n in range(0, len(p_i)):
        if unp.nominal_values(p_i)[n] < 100:
                sig_p_i[n] = unp.nominal_values(p_i)[n]*0.3
        elif unp.nominal_values(p_i)[n] >= 100:
                sig_p_i[n] = unp.nominal_values(p_i)[n]*0.5

    S_i = leckrate(t_i, unp.nominal_values(p_i), sig_p_i, pg, go, ''+x+'mbar', 'dreh')
    S[i] = unp.nominal_values(S_i)
    sig_S[i] = unp.std_devs(S_i)

    np.savetxt('data/tab_dreh_leck_'+x+'mbar.txt', np.column_stack([t_i, p1_i, sig_p1_i, p2_i, sig_p2_i, p3_i, sig_p3_i, unp.nominal_values(p_i),sig_p_i]),fmt=['%.1f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f','%.3f', '%.3f', '%.3f'], header='t, m1, m2, m3, middel')


#-------------------------------------------------------------Turbopumpe---------------------------------------------------------------------
print('_____________________________________Turbopumpe________________________________________________________________________________ \n')
print('EINHEIT IST NANO BAR ALSO 10^-9 bar 10^-6 mbar \n')

names = ['1e-4', '2e-4', '5e-5', '7e-5'] 
numbs = [0, 1, 2, 3]
p_g = [ ufloat(100, 50), ufloat(200, 60), ufloat(50, 15), ufloat(70, 21)]
gogos = [20, 20, 20, 20]

S_t = np.zeros(4)
sig_S_t = np.zeros(4) 

for i, x, pg, go in zip(numbs, names, p_g, gogos):

    t_i, p1_i, p2_i, p3_i = np.genfromtxt('data/Turbo_Leckrate_'+x+'.txt', unpack=True)

    if x=='5e-5':
        p1_i *= 100
        p2_i *= 100
        p3_i *= 100
    else:
        p1_i *= 1000
        p2_i *= 1000
        p3_i *= 1000

    sig_p1_i = np.zeros(len(p1_i))
    sig_p2_i = np.zeros(len(p2_i))
    sig_p3_i = np.zeros(len(p3_i))

    drucks = [p1_i, p2_i, p3_i]
    sigms = [sig_p1_i, sig_p2_i, sig_p3_i]

    for p, sig in zip(drucks, sigms):
        for n in range(len(t_i)):
            sig[n] = p[n]*0.3

    p_i = middel(p1_i, p2_i, p3_i)

    sig_p_i = np.zeros(len(p_i))
    for n in range(0, len(p_i)):
        sig_p_i[n] = unp.nominal_values(p_i)[n]*0.3


    S_i = leckrate(t_i, unp.nominal_values(p_i), sig_p_i, pg, go, ''+x+'mbar', 'turbo')
    S_t[i] = unp.nominal_values(S_i)
    sig_S_t[i] = unp.std_devs(S_i)

    np.savetxt('data/tab_turbo_leck_'+x+'mbar.txt', np.column_stack([t_i, p1_i, sig_p1_i, p2_i, sig_p2_i, p3_i, sig_p3_i, unp.nominal_values(p_i),sig_p_i]),fmt=['%.1f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f','%.3f', '%.3f', '%.3f'], header='t, m1, m2, m3, middel')


#-------------------------------------------------------------Abbildung für die Diskussion---------------------------------------------------

plt.figure()
#Theoriebereich
fig,ax = plt.subplots()
ax.axhspan(1.277, 1.5277, xmin=0, xmax=1, alpha=0.5, color='black', label='Theoriewert' )
#Evakuierungskurve
plt.errorbar(407.5, unp.nominal_values(S_eva)[0], yerr=unp.std_devs(S_eva)[0], xerr=(800+15)/2, fmt='+', color='red', label='Evakuierung 1' )
plt.errorbar(8.25,  unp.nominal_values(S_eva)[1], yerr=unp.std_devs(S_eva)[1], xerr=(15+1.5)/2, fmt='+', color='orange', label='Evakuierung 2' )
plt.errorbar(0.975, unp.nominal_values(S_eva)[2], yerr=unp.std_devs(S_eva)[2], xerr=(0.45+1.5)/2, fmt='+', color='yellow', label='Evakuierung 3' )
#Leckratenmessung
plt.errorbar(100, S[3], yerr=sig_S[3], xerr=50, fmt='+',color='blue', label=r'$p_\text{G} = \SI{100(50)}{\milli\bar}$' )
plt.errorbar( 50, S[2], yerr=sig_S[2], xerr=15, fmt='+',color='green', label=r'$p_\text{G} = \SI{50(15)}{\milli\bar}$' )
plt.errorbar( 10, S[1], yerr=sig_S[1], xerr=3, fmt='+',color='magenta', label=r'$p_\text{G} = \SI{10(3)}{\milli\bar}$' )
plt.errorbar(0.5, S[0], yerr=sig_S[0], xerr=0.15, fmt='+',color='brown', label=r'$p_\text{G} = \SI{5.0(15)e-1}{\milli\bar}$' )
plt.rc('axes', labelsize=size_label)
plt.xlabel(r'$p \mathbin{/} \si{\milli\bar}$')
plt.ylabel(r'$S \mathbin{/} \left(\si{\litre\per\second}\right)$')
plt.xlim(0, 450)
plt.legend()
plt.tight_layout()
plt.savefig('build/compare_dreh.pdf')

np.savetxt('data/tab_dreh_saug.txt', np.column_stack([unp.nominal_values(S_eva)[0], unp.std_devs(S_eva)[0], unp.nominal_values(S_eva)[1], unp.std_devs(S_eva)[1], unp.nominal_values(S_eva)[2], unp.std_devs(S_eva)[2], S[3], sig_S[3], S[2], sig_S[2], S[1], sig_S[1], S[0], sig_S[0]]), fmt=['%.4f', '%.4f', '%.4f', '%.4f', '%.4f','%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f'])   


plt.figure()
#Theoriebereich
plt.hlines(77, xmin=0, xmax=5000, color='black', label='Theoriewert')
#Evakuierungskurve
plt.errorbar(980, unp.nominal_values(S_eva_t)[0], yerr=unp.std_devs(S_eva_t)[0], xerr=(5000+40)/2, fmt='+', color='red', label='Evakuierung 1' )
plt.errorbar(29,   unp.nominal_values(S_eva_t)[1], yerr=unp.std_devs(S_eva_t)[1], xerr=(40+18)/2, fmt='+', color='orange', label='Evakuierung 2' )
plt.errorbar(15.5, unp.nominal_values(S_eva_t)[2], yerr=unp.std_devs(S_eva_t)[2], xerr=(18+13)/2, fmt='+', color='yellow', label='Evakuierung 3' )
#Leckratenmessung
plt.errorbar(100, S[0], yerr=sig_S_t[0], xerr=30, fmt='+',color='blue', label=r'$p_\text{G} = \SI{1.0(3)e-4}{\milli\bar}$' )
plt.errorbar(200, S[1], yerr=sig_S_t[1], xerr=60, fmt='+',color='green', label=r'$p_\text{G} = \SI{2.0(6)e-4}{\milli\bar}$' )
plt.errorbar(50, S[2], yerr=sig_S_t[2], xerr=15, fmt='+',color='magenta', label=r'$p_\text{G} = \SI{5.0(25)e-5}{\milli\bar}$' )
plt.errorbar(70, S[3], yerr=sig_S_t[3], xerr=21, fmt='+',color='brown', label=r'$p_\text{G} = \SI{7.0(21)e-5}{\milli\bar}$' )
plt.rc('axes', labelsize=size_label)
plt.xlabel(r'$p \mathbin{/} \SI{e-6}{\milli\bar}$')
plt.ylabel(r'$S \mathbin{/} \left(\si{\litre\per\second}\right)$')
plt.xlim(0, 1000)
plt.legend()
plt.tight_layout()
plt.savefig('build/compare_turbo.pdf')

np.savetxt('data/tab_turbo_saug.txt', np.column_stack([unp.nominal_values(S_eva_t)[0], unp.std_devs(S_eva_t)[0], unp.nominal_values(S_eva_t)[1], unp.std_devs(S_eva_t)[1], unp.nominal_values(S_eva_t)[2], unp.std_devs(S_eva_t)[2], S_t[0], sig_S_t[0], S_t[1], sig_S_t[1], S_t[2], sig_S_t[2], S_t[3], sig_S_t[3]]), fmt=['%.4f', '%.4f', '%.4f', '%.4f', '%.4f','%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f'])   