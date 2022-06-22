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
p_7 *= 100
p_G_7 = ufloat(70, 21)

S_7 = leckrate(t_7, p_7, p_G_7, '70nbar', 'turbo')