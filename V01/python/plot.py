import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

import os

if os.path.exists("build") == False:
    os.mkdir("build")

print('##################################################################################################################')
print('Justieren der Verzögerungsleitungen')
print('\n')

del_t, Counts = np.genfromtxt('data/Verzoegerungsleitung_Justage.txt', unpack=True)
del_t_links = del_t[3:17]
Counts_links = Counts[3:17]

del_t_plateau = np.delete(del_t[17:25],6)
Counts_plateau = np.delete(Counts[17:25],6)
mittel_plateau = np.mean(Counts_plateau)

del_t_rechts = np.delete(del_t[23:],1 )
Counts_rechts = np.delete(Counts[23:],1)

def line(t,b,y):
    return b*t + y


popt, pcov = curve_fit(line, del_t_links, Counts_links)
errors = np.sqrt(np.diag(pcov))
print('Verzögerung links mit der Form: b*t + y')
print("b =", popt[0])
print("b fehler =",errors[0])
b_links = ufloat(popt[0],errors[0])
print("y =", popt[1])
print("y fehler =",errors[1])
y_links = ufloat(popt[1],errors[1])
t_h_p_links = ((mittel_plateau/2) - y_links)/b_links
print('t_1/2_links : ', t_h_p_links)
print('\n')

nst_links = -popt[1]/popt[0]

x_links = np.linspace(nst_links,0,100)

#plt.plot(del_t, Counts,'k+',label='Messpunkte')
#plt.plot(del_t_links, Counts_links,'r+',label='Gerade links')
plt.errorbar(del_t_links, Counts_links, yerr=np.sqrt(Counts_links),fmt='r.',ecolor='r', label='Messwerte links')
plt.plot(x_links, line(x_links,popt[0],popt[1]),'r',label ='Fit Links' ) #linspace
#plt.xlabel(r'$\alpha \:/\: \si{\ohm}$')
#plt.ylabel(r'$y \:/\: \si{\micro\joule}$')


#plt.plot(del_t_rechts, Counts_rechts,'m+',label='Gerade rechts')
plt.errorbar(del_t_rechts, Counts_rechts, yerr=np.sqrt(Counts_rechts),fmt='m.',ecolor='m', label='Messwerte rechts')
popt, pcov = curve_fit(line, del_t_rechts, Counts_rechts)
errors = np.sqrt(np.diag(pcov))
print('Verzögerung rechts mit der Form: b*t + y')
print("b =", popt[0])
print("b fehler =",errors[0])
b_rechts = ufloat(popt[0],errors[0])
print("y =", popt[1])
print("y fehler =",errors[1])
y_rechts = ufloat(popt[1],errors[1])
t_h_p_rechts = ((mittel_plateau/2) - y_rechts)/b_rechts
print('t_1/2_rechts : ', t_h_p_rechts)
print('\n')

x_rechts = np.linspace(2.5,12.5,100)
plt.plot(x_rechts, line(x_rechts,popt[0],popt[1]),'m',label ='Fit rechts' ) #linspace


#plt.plot(del_t_plateau, Counts_plateau,'g+',label='Plateau')
plt.errorbar(del_t_plateau, Counts_plateau, yerr=np.sqrt(Counts_plateau),fmt='g.',ecolor='g', label='Messwerte Plateau')
print('mittel_plateau', mittel_plateau ,'+/-', np.sqrt(mittel_plateau))
plt.axhline(y=mittel_plateau, color='g', linestyle='-', label='Plateaumittelwert')
print('\n')
plt.axhline(y=mittel_plateau/2, color='g', linestyle='--', label='halber Plateaumittelwert')
print('t_1/2 : ', (t_h_p_rechts+t_h_p_links))

#plt.xlabel(r'$\Delta t/\si{\nano\second}$')
plt.xlabel('t / ns')
plt.ylabel('Ereignisse')
plt.grid()
plt.xlim(nst_links, 12.5)

# in matplotlibrc leider (noch) nicht möglich
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('build/Verzoegerungsleitung.pdf')
plt.clf()

print('##################################################################################################################')
print('Kalibration')
print('\n')

del_t_mca, Channel_alle, Counts_MCA = np.genfromtxt('data/Kalibration_MCA.txt', unpack=True)
del_t_MCA = np.delete(del_t_mca, 32)
Channel = np.delete(Channel_alle, 32)
plt.plot(Channel, del_t_MCA, 'k+', label='Messwerte')

popt, pcov = curve_fit(line, Channel, del_t_MCA)
errors = np.sqrt(np.diag(pcov))
print('Kalibration mit der Form: b*t + y')
print("b =", popt[0])
print("b fehler =",errors[0])
b_MCA = popt[0]
print("y =", popt[1])
print("y fehler =",errors[1])
y_MCA = popt[1]
print('\n')

x_kali=np.linspace(0,450,100)
plt.plot(x_kali, line(x_kali,popt[0],popt[1]),'r',label ='Fit' ) #linspace

plt.xlim(0,450)
plt.ylabel('t / ns')
plt.xlabel('Kanal')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('build/Kalibration_MCA.pdf')
plt.clf()

print('##################################################################################################################')
print('theoretischer Untergrund')
print('\n')


print('Channel-Intervall : [', popt[1] ,';' , (popt[0]*10)+popt[1] ,']')
T_S = 10*10**(-6)
t_Mess = 272201
N_Start = ufloat(4560816,2135.61)
I_Mess = N_Start/t_Mess
print('I_Mess : ', I_Mess)
mu = I_Mess*T_S
U_ges = mu * N_Start* unp.exp(-mu)
print('U_ges : ', U_ges)
U_norm = U_ges/(445-4)
print('U_norm : ', U_norm)

print('##################################################################################################################')
print('Bestimmung der Lebensdauer kosmischer Myonen')
print('\n')

Counts_Mess= np.genfromtxt('data/Myonen.Spe', unpack=True)
sieg_N = np.sqrt(Counts_Mess)

x_channel = []
x_zeit = []
length = len(Counts_Mess)
print(length)
for i in range(length):
    x_channel.append(i)
    x_zeit.append(line(i, b_MCA, y_MCA))
#print(x_channel)
#x_zeit = line(x_channel,b_MCA,y_MCA)
#print(x_zeit)

#plt.plot(x_channel, Counts_Mess,'g.',label='Messwerte')
plt.errorbar(x_channel, Counts_Mess, yerr=sieg_N,fmt='.',ecolor=None,markersize = 1.3, elinewidth = 0.42, label='Messwerte')

plt.ylabel('Anzahl an Ereignissen')
plt.xlabel('Kanal')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('build/Messung_Channel_Counts.pdf')
plt.clf()

plt.errorbar(x_zeit, Counts_Mess, yerr=sieg_N,fmt='.',ecolor=None,markersize = 1.3, elinewidth = 0.42, label='Messwerte')

x_zeit_raus = np.concatenate((x_zeit[:4],x_zeit[446:]), axis=None )
Counts_Mess_raus = np.concatenate((Counts_Mess[:4],Counts_Mess[446:]), axis=None )
plt.errorbar(x_zeit_raus,Counts_Mess_raus, yerr = np.sqrt(Counts_Mess_raus),fmt='k.',markersize = 1.3, elinewidth = 0.42,label='abgeschnittene Messwerte')

def e_Funktion(t,N,lamda,U_n):
    return N * np.exp( -lamda* np.array(t) ) + U_n

x_zeit_fit = x_zeit[4:446]
Counts_Mess_fit = Counts_Mess[4:446]
#plt.errorbar(x_zeit_fit, Counts_Mess_fit, yerr=sieg_N[4:446],fmt='.',ecolor='black', label='Messwerte')

popt, pcov = curve_fit(e_Funktion, x_zeit_fit, Counts_Mess_fit)
errors = np.sqrt(np.diag(pcov))
print('e-Funktion mit der Form: N*exp(-lambda*t)+U_n')
print("N =", popt[0])
print("N fehler =",errors[0])
print("lambda =", popt[1])
print("lambda fehler =",errors[1])
lambda_n = ufloat(popt[1],errors[1])
print("U_n =", popt[2])
print("U_n fehler =",errors[2])
U_num = ufloat(popt[2],errors[2])
print('\n')
print('Lebensdauer numerisch: ', 1/lambda_n)
tau_num = 1/lambda_n
print('\n')

x_zeit_lin = np.linspace(0,12,1000)
plt.plot(x_zeit_lin, e_Funktion(x_zeit_lin,popt[0],popt[1],popt[2]),'r',label ='Fit' ) #linspace

plt.ylabel('Anzahl an Ereignissen')
#plt.xlabel('t / $10^{-6}s$')
plt.xlabel('t / $\mu$s')

plt.xlim(0,12)
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('build/Messung_Zeit_Counts_num.pdf')
plt.clf()

Counts_Mess_fit = Counts_Mess_fit - unp.nominal_values(U_norm)

def eFkt_theo(t,N,lamda):
        return N * np.exp( -lamda* np.array(t) )

plt.errorbar(x_zeit, Counts_Mess, yerr=sieg_N,fmt='.',ecolor=None,markersize = 1.3, elinewidth = 0.42, label='bereinigte Messwerte')
plt.errorbar(x_zeit_raus,Counts_Mess_raus, yerr = np.sqrt(Counts_Mess_raus),fmt='k.',markersize = 1.3, elinewidth = 0.42,label='abgeschnittene Messwerte')

popt, pcov = curve_fit(eFkt_theo, x_zeit_fit, Counts_Mess_fit)
errors = np.sqrt(np.diag(pcov))
print('eFkt_theo mit der Form: N*exp(-lambda*t)')
print("N =", popt[0])
print("N fehler =",errors[0])
print("lambda =", popt[1])
print("lambda fehler =",errors[1])
lambda_t = ufloat(popt[1],errors[1])
print('\n')
print('Lebensdauer theoretisch: ', 1/lambda_t)
tau_theo = 1/lambda_t
print('\n')

plt.plot(x_zeit_lin, eFkt_theo(x_zeit_lin,popt[0],popt[1]),'r',label ='Fit' ) #linspace

plt.xlim(0,12)
plt.ylabel('Anzahl an Ereignissen')
plt.xlabel('t / $\mu$s')
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('build/Messung_Zeit_Counts_theo.pdf')
plt.clf()

print('Untergrund ABweichung: ', U_num/U_norm)
print('Lebensdauer Abweichung: ', tau_num/tau_theo)
print('Abweichung Lebensdauer theo Literatur: ', tau_theo/2.2)
print('Abweichung Lebensdauer num Literatur: ', tau_num/2.2)