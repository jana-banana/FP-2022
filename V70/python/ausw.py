import matplotlib.pyplot as plt
# from matplotlib.axes.Axes import axhspan
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
import uncertainties as unc

import os

if os.path.exists("../build") == False:
    os.mkdir("../build")
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
        arr[i] = unc.ufloat(np.mean([a[i], b[i], c[i]]), np.std([a[i], b[i], c[i]])/ np.sqrt(3))
    return arr

def andererfehler(a,b,c, d, e, f):
    err = unp.uarray(np.zeros(np.size(a)),np.zeros(np.size(a)))
    for i in range(0,np.size(a)):
        err[i] = unc.ufloat(np.mean([a[i], b[i], c[i]]), 1/3*np.sqrt(d[i]**2 + e[i]**2 + f[i]**2))
    return err

def evaku(t, la,g0, g1, g2, name):
    ln_line_1 = unp.nominal_values(la)[g0:g1]
    ln_line_2 = unp.nominal_values(la)[g1:g2]
    ln_line_3 = unp.nominal_values(la)[g2:]

    t_line_1 = t[g0:g1]
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
        print(f'Steigung: \n m = {popt[0]:.4f} \pm {errors[0]:.4f} \n')
        m[i] = popt[0]
        sig_m[i] = errors[0]
        print(f'Achsenabschnitt: \n n = {popt[1]:.4f} \pm {errors[1]:.4f} \n')
        n[i] = popt[1]
        sig_n[i] = errors[1]
    
    t_plot1 = np.linspace(-10 + g0, g1*10 + 50, 1000)
    t_plot2 = np.linspace(g1*10 - 50, g2*10 +50, 1000)
    t_plot3 = np.linspace(g2*10-50, 610, 1000)
    if len(t)<20:
        t_plot1 = np.linspace(-10 + g0, g1*10 + 10, 1000)
        t_plot2 = np.linspace(g1*10 - 10, g2*10 +10, 1000)
        t_plot3 = np.linspace(g2-10, 130, 1000)

    plt.figure()
    plt.errorbar(t, unp.nominal_values(la), yerr=unp.std_devs(la), fmt='.k')
    plt.errorbar(t_line_1, ln_line_1, yerr=unp.std_devs(la)[g0:g1], fmt='.', color='red')
    plt.plot(t_plot1, line(t_plot1, m[0], n[0]), '-b')
    plt.errorbar(t_line_2, ln_line_2, yerr=unp.std_devs(la)[g1:g2], fmt='.', color='orange')
    plt.plot(t_plot2, line(t_plot2, m[1], n[1]), '-g')
    plt.errorbar(t_line_3, ln_line_3, yerr=unp.std_devs(la)[g2:], fmt='.', color='yellow')
    plt.plot(t_plot3, line(t_plot3, m[2], n[2]), '-m')
    plt.xlabel('t / s')
    plt.ylabel('logarithmusausdruck')
    plt.savefig('../build/evaku'+name+'.pdf')

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
            print(f'Bereich{i+1} S= {S[i]:.4f} \pm {sig_S[i]:.4f}, abweichung {abw[i]}')
    else :
        abw_lower = fehler(S, 4.6*10/36)
        abw_higher = fehler(S, 5.5*10/36)
        for i in nums:
            print(f'Bereich{i+1} S= {S[i]:.4f} \pm {sig_S[i]:.4f}, abweichung {abw_lower[i]}, {abw_higher[i]}')
    print('\n')
    
    return unp.uarray(S, sig_S)

def leckrate(t, p, p_G, drueck, name):
    t_plot = np.linspace(-10, 210, 1000)
    if name =='turbo':
        t_plot = np.linspace(-5, 125, 1000)

    popt, pcov = curve_fit(line, t, unp.nominal_values(p))
    errors = np.sqrt(np.diag(pcov))

    print('Gleichgewichtsdruck',drueck, name)
    print(f'Steigung: \n m = {popt[0]} \pm {errors[0]} \n')
    m = ufloat(popt[0], errors[0])
    print(f'Achsenabschnitt: \n n = {popt[1]} \pm {errors[1]} \n')
    n = ufloat(popt[1], errors[1])

    plt.figure()
    plt.errorbar(t, unp.nominal_values(p), yerr=unp.std_devs(p), fmt='.k', label='Messwerte')
    plt.plot(t_plot, line(t_plot, m.n, n.n ), '-r', label='Ausgleichsgerade')
    plt.xlabel('t/s')
    plt.ylabel('p / mbar')
    plt.legend()
    plt.tight_layout()
    plt.savefig('../build/leck_'+name+'_'+drueck+'.pdf')

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


#------------------------------------------------------------ Drehschieberpumpe--------------------------------------------------------------
print('_____________________________________Drehschieberpumpe_________________________________________________________________________ \n')
print('_____________________Evakuierungskurve_______________________________ \n')
t_eva, p1_eva, p2_eva, p3_eva = np.genfromtxt("../data/Dreh_Evakuierung.txt", unpack=True)

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
        elif p[n] > 100:
            sig[n] = p[n]*0.5
        elif p[n] == 100:
            print("well, hm")
        else:
            print('definitiv ein fehler')



ln1 = logausdruck(p1_eva, sig_p1_eva ,p_E, p_0)
ln2 = logausdruck(p2_eva, sig_p2_eva ,p_E, p_0)
ln3 = logausdruck(p3_eva, sig_p3_eva ,p_E, p_0)


#####Auswertungsfunktion anwenden

S_1 = evaku(t_eva, ln1, 0, 13, 31, 'dreh_1')
S_2 = evaku(t_eva, ln2, 0, 13, 31, 'dreh_2')
S_3 = evaku(t_eva, ln3, 2, 13, 31, 'dreh_3')

##Mittelwerte ausrechnen
print('--------Mittelwerte der einzelnen Abschnitte--------- \n ')
smean = middel(unp.nominal_values(S_1), unp.nominal_values(S_2), unp.nominal_values(S_3))

abw_lower = fehler(unp.nominal_values(smean), 4.6*10/36)
abw_higher = fehler(unp.nominal_values(smean), 5.5*10/36)
for i in range(3):
    print(f'Bereich{i+1} S= {unp.nominal_values(smean)[i]} \pm {unp.std_devs(smean)[i]}, abweichung {abw_lower[i]}, {abw_higher[i]}')
print('\n')

smeann = andererfehler(unp.nominal_values(S_1), unp.nominal_values(S_2), unp.nominal_values(S_3), unp.std_devs(S_1), unp.std_devs(S_2), unp.std_devs(S_3))

abw_lower = fehler(unp.nominal_values(smeann), 4.6*10/36)
abw_higher = fehler(unp.nominal_values(smeann), 5.5*10/36)
for i in range(3):
    print(f'Bereich{i+1} S= {unp.nominal_values(smeann)[i]} \pm {unp.std_devs(smeann)[i]}, abweichung {abw_lower[i]}, {abw_higher[i]}')
print('\n')

np.savetxt('../data/tab_dreh_eva.txt', np.column_stack([t_eva, p1_eva, sig_p1_eva, unp.nominal_values(ln1), unp.std_devs(ln1),p2_eva, sig_p2_eva, unp.nominal_values(ln2), unp.std_devs(ln2), p3_eva, sig_p3_eva, unp.nominal_values(ln3), unp.std_devs(ln3)]),fmt=['%.1f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f','%.4f', '%.4f', '%.4f', '%.4f','%.4f', '%.4f', '%.4f'], header='t, m1, lnm1, m2, lnm2, m3, lnm3')


print('__________________Leckratenmessung_________________________')

t_100, p1_100, p2_100, p3_100 = np.genfromtxt('../data/Dreh_Leckrate_100.txt', unpack=True)

sig_p1_100 = np.zeros(len(p1_100))
sig_p2_100 = np.zeros(len(p2_100))
sig_p3_100 = np.zeros(len(p3_100))

drucks = [p1_100, p2_100, p3_100]
sigms = [sig_p1_100, sig_p2_100, sig_p3_100]

for p, sig in zip(drucks, sigms):
    for n in range(len(t_100)):
        if p[n] < 100:
            sig[n] = p[n]*0.3
        elif p[n] > 100:
            sig[n] = p[n]*0.5
        elif p[n] == 100:
            sig[n] = p[n]*0.3
        else:
            print('definitiv ein fehler')

p_100 = middel(p1_100, p2_100, p3_100)
p_G_100 = ufloat(100, 30)

t_plot_100 = np.linspace(-10, 110, 1000)

S_100 = leckrate(t_100, p_100, p_G_100, '100mbar', 'dreh')

np.savetxt('../data/tab_dreh_leck_100.txt', np.column_stack([t_100, p1_100, sig_p1_100, p2_100, sig_p2_100, p3_100, sig_p3_100, unp.nominal_values(p_100), unp.std_devs(p_100)]),fmt=['%.1f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f','%.2f', '%.4f', '%.4f'], header='t, m1, m2, m3, middel')

# popt, pcov = curve_fit(line, t_100[0:11], unp.nominal_values(p_100)[0:11])
# errors = np.sqrt(np.diag(pcov))
# 
# print('Gleichgewichtsdruck 100 Dreh') 
# print(f'Steigung: \n m = {popt[0]} \pm {errors[0]} \n')
# m = ufloat(popt[0], errors[0])
# print(f'Achsenabschnitt: \n n = {popt[1]} \pm {errors[1]} \n')
# n = ufloat(popt[1], errors[1])
# 
# plt.figure()
# plt.errorbar(t_100[0:11], unp.nominal_values(p_100)[0:11], yerr=unp.std_devs(p_100)[0:11], fmt='.k', label='benutzte Messwerte')
# plt.errorbar(t_100[11:], unp.nominal_values(p_100)[11:], yerr=unp.std_devs(p_100)[11:], fmt='.b', label='nicht benutzte Messwerte')
# plt.plot(t_plot_100, line(t_plot_100, m.n, n.n ), '-r', label='Ausgleichsgerade')
# plt.xlabel('t/s')
# plt.ylabel('p / mbar')
# plt.legend()
# plt.tight_layout()
# plt.savefig('../build/leck_dreh100.pdf')
# 
##Saugvermögen
# print("--------Saugvermögen--dreh 100------- \n ")
# V = 34 #liter
# sig_V = 3.4
# 
# S = m.n*V/p_G_100.n
# sig_S = np.sqrt((V/p_G_100.n)**2*m.s**2 + (m.n/p_G_100.n)**2*sig_V**2 + (m.n*V/(p_G_100.n**2))**2*p_G_100.s**2)
# 
# abw_lower = fehler(S, 4.6*10/36)
# abw_higher = fehler(S, 5.5*10/36)
# print(f'Gleichgwichtsdruck:{p_G_100.n} S= {S} \pm {sig_S}, abweichung {abw_lower}, {abw_higher}')
# print('\n')
# 
# S_100 = unp.uarray(S, sig_S)

t_50, p1_50, p2_50, p3_50 = np.genfromtxt('../data/Dreh_Leckrate_50.txt', unpack=True)

sig_p1_50 = np.zeros(len(p1_50))
sig_p2_50 = np.zeros(len(p2_50))
sig_p3_50 = np.zeros(len(p3_50))

drucks = [p1_50, p2_50, p3_50]
sigms = [sig_p1_50, sig_p2_50, sig_p3_50]

for p, sig in zip(drucks, sigms):
    for n in range(len(t_50)):
        if p[n] < 100:
            sig[n] = p[n]*0.3
        elif p[n] > 100:
            sig[n] = p[n]*0.5
        elif p[n] == 100:
            print("well, hm")
        else:
            print('definitiv ein fehler')

p_50 = middel(p1_50, p2_50, p3_50)
p_G_50 = ufloat(50, 15)

S_50 = leckrate(t_50, p_50, p_G_50, '50mbar', 'dreh')

np.savetxt('../data/tab_dreh_leck_50.txt', np.column_stack([t_50, p1_50, sig_p1_50, p2_50, sig_p2_50, p3_50, sig_p3_50, unp.nominal_values(p_50), unp.std_devs(p_50)]),fmt=['%.1f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f','%.2f', '%.4f', '%.4f'], header='t, m1, m2, m3, middel')

t_10, p1_10, p2_10, p3_10 = np.genfromtxt('../data/Dreh_Leckrate_10.txt', unpack=True)

sig_p1_10 = np.zeros(len(p1_10))
sig_p2_10 = np.zeros(len(p2_10))
sig_p3_10 = np.zeros(len(p3_10))

drucks = [p1_10, p2_10, p3_10]
sigms = [sig_p1_10, sig_p2_10, sig_p3_10]

for p, sig in zip(drucks, sigms):
    for n in range(len(t_10)):
        if p[n] < 100:
            sig[n] = p[n]*0.3
        elif p[n] > 100:
            sig[n] = p[n]*0.5
        elif p[n] == 100:
            print("well, hm")
        else:
            print('definitiv ein fehler')

p_10 = middel(p1_10, p2_10, p3_10)
p_G_10 = ufloat(10, 3)

S_10 = leckrate(t_10, p_10, p_G_10, '10mbar', 'dreh')

np.savetxt('../data/tab_dreh_leck_10.txt', np.column_stack([t_10, p1_10, sig_p1_10, p2_10, sig_p2_10, p3_10, sig_p3_10, unp.nominal_values(p_10), unp.std_devs(p_10)]),fmt=['%.1f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f','%.2f', '%.4f', '%.4f'], header='t, m1, m2, m3, middel')


t_05, p1_05, p2_05, p3_05 = np.genfromtxt('../data/Dreh_Leckrate_05.txt', unpack=True)

sig_p1_05 = np.zeros(len(p1_05))
sig_p2_05 = np.zeros(len(p2_05))
sig_p3_05 = np.zeros(len(p3_05))

drucks = [p1_05, p2_05, p3_05]
sigms = [sig_p1_05, sig_p2_05, sig_p3_05]

for p, sig in zip(drucks, sigms):
    for n in range(len(t_05)):
        if p[n] < 100000:
            sig[n] = p[n]*0.3
        elif p[n] > 100000:
            sig[n] = p[n]*0.5
        elif p[n] == 100000:
            print("well, hm")
        else:
            print('definitiv ein fehler')

p_05 = middel(p1_05, p2_05, p3_05)
p_05 *= 0.001
p_G_05 = ufloat(0.5, 0.15)

S_05 = leckrate(t_05, p_05, p_G_05, '05mbar', 'dreh')

np.savetxt('../data/tab_dreh_leck_05.txt', np.column_stack([t_05, p1_05, sig_p1_05, p2_05, sig_p2_05, p3_05, sig_p3_05, unp.nominal_values(p_05), unp.std_devs(p_05)]),fmt=['%.1f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f','%.4f', '%.4f', '%.4f'], header='t, m1, m2, m3, middel')

####Abbildung für die Diskussion

plt.figure()
#Theoriebereich
fig,ax = plt.subplots()
ax.axhspan(1.277, 1.5277, xmin=0, xmax=1, alpha=0.5, color='black', label='Theoriewert' )
#Evakuierungskurve
plt.errorbar(407.5, unp.nominal_values(smean)[0], yerr=unp.std_devs(smean)[0], xerr=(800+15)/2, fmt='+', color='red', label='Evakuierung 1' )
plt.errorbar(8.25, unp.nominal_values(smean)[1], yerr=unp.std_devs(smean)[1], xerr=(15+1.5)/2, fmt='+', color='orange', label='Evakuierung 2' )
plt.errorbar(0.975, unp.nominal_values(smean)[2], yerr=unp.std_devs(smean)[2], xerr=(0.45+1.5)/2, fmt='+', color='yellow', label='Evakuierung 3' )
#Leckratenmessung
plt.errorbar(550, unp.nominal_values(S_100), yerr=unp.std_devs(S_100), xerr=(1000+100)/2, fmt='+',color='blue', label='Leckrate 100 mbar' )
plt.errorbar(525, unp.nominal_values(S_50), yerr=unp.std_devs(S_50), xerr=(1000+50)/2, fmt='+',color='green', label='Leckrate 50 mbar' )
plt.errorbar(35, unp.nominal_values(S_10), yerr=unp.std_devs(S_10), xerr=(10+60)/2, fmt='+',color='magenta', label='Leckrate 10 mbar' )
plt.errorbar(0.75, unp.nominal_values(S_05), yerr=unp.std_devs(S_05), xerr=(1+0.5)/2, fmt='+',color='brown', label='Leckrate 0.5 mbar' )

plt.xlabel('p/ mbar')
plt.ylabel('S/ (L/s)')
plt.legend()
plt.tight_layout()
plt.savefig('../build/compare_dreh.pdf')

np.savetxt('../data/tab_dreh_saug.txt', np.column_stack([unp.nominal_values(smean)[0], unp.std_devs(smean)[0], unp.nominal_values(smean)[1], unp.std_devs(smean)[1], unp.nominal_values(smean)[2], unp.std_devs(smean)[2], unp.nominal_values(S_100), unp.std_devs(S_100), unp.nominal_values(S_50), unp.std_devs(S_50), unp.nominal_values(S_10), unp.std_devs(S_10), unp.nominal_values(S_05), unp.std_devs(S_05)]), fmt=['%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f','%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f'])   
#------------------------------------------------------------ Turbopumpe--------------------------------------------------------------
print('_____________________________________Turbopumpe_________________________________________________________________________ \n')
print('_____________________Evakuierungskurve_______________________________ \n')
print('EINHEIT IST NANO BAR ALSO 10^-9 bar 10^-6 mbar')

t_eva_t, p1_eva_t, p2_eva_t, p3_eva_t = np.genfromtxt("../data/Turbo_Evakuierung.txt", unpack=True)
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

ln1_t = logausdruck(p1_eva_t, sig_p1_eva_t, p_E, p_0)
ln2_t = logausdruck(p2_eva_t, sig_p2_eva_t, p_E, p_0)
ln3_t = logausdruck(p3_eva_t, sig_p3_eva_t, p_E, p_0)


#####Auswertungsfunktion anwenden

S_1_t = evaku(t_eva_t, ln1_t, 0, 3, 6, 'turbo_1')
S_2_t = evaku(t_eva_t, ln2_t, 0, 3, 6, 'turbo_2')
S_3_t = evaku(t_eva_t, ln3_t, 0, 3, 6, 'turbo_3')

##Mittelwerte ausrechnen
print('--------Mittelwerte der einzelnen Abschnitte--------- \n ')
smean = middel(unp.nominal_values(S_1_t), unp.nominal_values(S_2_t), unp.nominal_values(S_3_t))

theo = 77 #liters per second
abw = fehler(smean, theo)
for i in range(3):
    print(f'Bereich{i+1} S= {unp.nominal_values(smean)[i]} \pm {unp.std_devs(smean)[i]}, abweichung {abw[i]}')
print('\n')

np.savetxt('../data/tab_turbo_eva.txt', np.column_stack([t_eva_t, p1_eva_t, sig_p1_eva_t, unp.nominal_values(ln1_t), unp.std_devs(ln1_t),p2_eva_t, sig_p2_eva_t, unp.nominal_values(ln2_t), unp.std_devs(ln2_t), p3_eva_t, sig_p3_eva_t, unp.nominal_values(ln3_t), unp.std_devs(ln3_t)]), fmt=['%.1f', '%.2f', '%.2f', '%.4f', '%.4f', '%.2f','%.2f', '%.4f', '%.4f', '%.2f','%.2f', '%.4f', '%.4f'], header='t, m1. lnm1, m2, lnm2, m3, lnm3')


print('__________________Leckratenmessung_________________________')

t_1, p1_1, p2_1, p3_1 = np.genfromtxt('../data/Turbo_Leckrate_1e-4.txt', unpack=True)

sig_p1_1 = np.zeros(len(p1_1))
sig_p2_1 = np.zeros(len(p2_1))
sig_p3_1 = np.zeros(len(p3_1))

drucks = [p1_1, p2_1, p3_1]
sigms = [sig_p1_1, sig_p2_1, sig_p3_1]

for p, sig in zip(drucks, sigms):
    for n in range(len(t_1)):
        sig[n] = p[n]*0.3
       

p_1 = middel(p1_1, p2_1, p3_1)
p_1 *= 1000
p_G_1 = ufloat(100, 30)

S_1_t = leckrate(t_1, p_1, p_G_1, '100nbar', 'turbo')

np.savetxt('../data/tab_turbo_leck_1.txt', np.column_stack([t_1, p1_1*1000, sig_p1_1*1000, p2_1*1000, sig_p2_1*1000, p3_1*1000, sig_p3_1*1000, unp.nominal_values(p_1), unp.std_devs(p_1)]),fmt=['%.1f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f','%.2f', '%.4f', '%.4f'], header='t, m1, m2, m3, middel')

t_2, p1_2, p2_2, p3_2 = np.genfromtxt('../data/Turbo_Leckrate_2e-4.txt', unpack=True)

sig_p1_2 = np.zeros(len(p1_2))
sig_p2_2 = np.zeros(len(p2_2))
sig_p3_2 = np.zeros(len(p3_2))

drucks = [p1_2, p2_2, p3_2]
sigms = [sig_p1_2, sig_p2_2, sig_p3_2]

for p, sig in zip(drucks, sigms):
    for n in range(len(t_2)):
        sig[n] = p[n]*0.3
       

p_2 = middel(p1_2, p2_2, p3_2)
p_2 *= 1000
p_G_2 = ufloat(200, 60)

S_2_t = leckrate(t_2, p_2, p_G_2, '200nbar', 'turbo')
np.savetxt('../data/tab_turbo_leck_2.txt', np.column_stack([t_2, p1_2*1000, sig_p1_2*1000, p2_2*1000, sig_p2_2*1000, p3_2*1000, sig_p3_2*1000, unp.nominal_values(p_2), unp.std_devs(p_2)]),fmt=['%.1f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f','%.2f', '%.4f', '%.4f'], header='t, m1, m2, m3, middel')

t_5, p1_5, p2_5, p3_5 = np.genfromtxt('../data/Turbo_Leckrate_5e-5.txt', unpack=True)

sig_p1_5 = np.zeros(len(p1_5))
sig_p2_5 = np.zeros(len(p2_5))
sig_p3_5 = np.zeros(len(p3_5))

drucks = [p1_5, p2_5, p3_5]
sigms = [sig_p1_5, sig_p2_5, sig_p3_5]

for p, sig in zip(drucks, sigms):
    for n in range(len(t_5)):
        sig[n] = p[n]*0.3
       

p_5 = middel(p1_5, p2_5, p3_5)
p_5 *= 100
p_G_5 = ufloat(50, 15)

S_5_t = leckrate(t_5, p_5, p_G_5, '50nbar', 'turbo')
np.savetxt('../data/tab_turbo_leck_5.txt', np.column_stack([t_5, p1_5*100, sig_p1_5*100, p2_5*100, sig_p2_5*100, p3_5*100, sig_p3_5*100, unp.nominal_values(p_5), unp.std_devs(p_5)]),fmt=['%.1f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f','%.2f', '%.4f', '%.4f'], header='t, m1, m2, m3, middel')

t_7, p1_7, p2_7, p3_7 = np.genfromtxt('../data/Turbo_Leckrate_7e-5.txt', unpack=True)

sig_p1_7 = np.zeros(len(p1_7))
sig_p2_7 = np.zeros(len(p2_7))
sig_p3_7 = np.zeros(len(p3_7))

drucks = [p1_7, p2_7, p3_7]
sigms = [sig_p1_7, sig_p2_7, sig_p3_7]

for p, sig in zip(drucks, sigms):
    for n in range(len(t_7)):
        sig[n] = p[n]*0.3
       

p_7 = middel(p1_7, p2_7, p3_7)
p_7 *= 1000
p_G_7 = ufloat(70, 21)

S_7_t = leckrate(t_7, p_7, p_G_7, '70nbar', 'turbo')
np.savetxt('../data/tab_turbo_leck_7.txt', np.column_stack([t_7, p1_7*1000, sig_p1_7*1000, p2_7*1000, sig_p2_7*1000, p3_7*1000, sig_p3_7*1000, unp.nominal_values(p_7), unp.std_devs(p_7)]),fmt=['%.1f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f','%.2f', '%.4f', '%.4f'], header='t, m1, m2, m3, middel')

plt.figure()
#Theoriebereich
plt.hlines(77, xmin=0, xmax=5000, color='black', label='Theoriewert')
#Evakuierungskurve
plt.errorbar(2520, unp.nominal_values(smean)[0], yerr=unp.std_devs(smean)[0], xerr=(5000+40)/2, fmt='+', color='red', label='Evakuierung 1' )
plt.errorbar(29, unp.nominal_values(smean)[1], yerr=unp.std_devs(smean)[1], xerr=(40+18)/2, fmt='+', color='orange', label='Evakuierung 2' )
plt.errorbar(15.5, unp.nominal_values(smean)[2], yerr=unp.std_devs(smean)[2], xerr=(18+13)/2, fmt='+', color='yellow', label='Evakuierung 3' )
#Leckratenmessung
plt.errorbar(2550, unp.nominal_values(S_1_t), yerr=unp.std_devs(S_1_t), xerr=(5000+100)/2, fmt='+',color='blue', label='Leckrate 100 nbar' )
plt.errorbar(2600, unp.nominal_values(S_2_t), yerr=unp.std_devs(S_2_t), xerr=(5000+200)/2, fmt='+',color='green', label='Leckrate 200 nbar' )
plt.errorbar(97.5, unp.nominal_values(S_5_t), yerr=unp.std_devs(S_5_t), xerr=(50+145)/2, fmt='+',color='magenta', label='Leckrate 50 nbar' )
plt.errorbar(1135, unp.nominal_values(S_7_t), yerr=unp.std_devs(S_7_t), xerr=(70+2200)/2, fmt='+',color='brown', label='Leckrate 70 nbar' )

plt.xlabel('p/ nbar')
plt.ylabel('S/ (L/s)')
plt.legend()
plt.tight_layout()
plt.savefig('../build/compare_turbo.pdf')

np.savetxt('../data/tab_turbo_saug.txt', np.column_stack([unp.nominal_values(smean)[0], unp.std_devs(smean)[0], unp.nominal_values(smean)[1], unp.std_devs(smean)[1], unp.nominal_values(smean)[2], unp.std_devs(smean)[2], unp.nominal_values(S_1_t), unp.std_devs(S_1_t), unp.nominal_values(S_2_t), unp.std_devs(S_2_t), unp.nominal_values(S_5_t), unp.std_devs(S_5_t), unp.nominal_values(S_7_t), unp.std_devs(S_7_t)]), fmt=['%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f', '%.4f'])   