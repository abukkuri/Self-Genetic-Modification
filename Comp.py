import numpy as np
import matplotlib.pyplot as plt
import random

num_sims = 100

dd_thres = 0.2 #0.2, 0.4, 0.6
mut = 0.05
dose = 0.4 #.6, .8, 1

g = 0.02
K = 200
r = 0.6
l1 = 1
b1 = 1
c21 = 0.7
c12 = 0.2
rho = 0.7

tend = 800

norm_macros = []
pacc_macros = []
strat_macros = []
nstrat_macros = []

norm_macroe = []
pacc_macroe = []
strat_macroe = []
nstrat_macroe = []

t_macro = []
extincts = []
extincte = []

def drug_death(current_strat,m1): #here, if death due to drug is >dd_thres, we won't transition into 2N+
        return m1/(l1+b1*current_strat)

def run_this():
        normale = [100]
        pacce = [0]
        strate = [0]

        normals = [100]
        paccs = [0]
        strats = [0]
        nstrats = [0]

        ee = False
        es = False

        t = [0]

        while t[-1] < tend:
                if t[-1] > 0 and t[-1] < 600:
                        m1 = dose
                else:
                        m1 = 0

                current_ne = normale[-1]
                current_pacce = pacce[-1]
                current_strate = strate[-1]

                current_ns = normals[-1]
                current_paccs = paccs[-1]
                current_strats = strats[-1]

                #broken into birth, death, switching

                rates = [(current_ns*r), current_ns*(m1/(l1+b1*current_strats)+r*((current_paccs+current_ns+current_ne+current_pacce)/K)),current_ns*(c21*m1/(l1+b1*current_strats)+g),0,0,current_paccs*c12,(current_ne*r), current_ne*(m1/(l1+b1*current_strate)+r*((current_paccs+current_ns+current_ne+current_pacce)/K)),current_ne*(c21*m1/(l1+b1*current_strate)+g),0,0,current_pacce*c12]
                rs = [(current_ns*r), current_ns*(m1/(l1+b1*current_strats)+r*((current_paccs+current_ns+current_ne+current_pacce)/K)),current_ns*(c21*m1/(l1+b1*current_strats)+g),0,0,current_paccs*c12]
                re = [(current_ne*r), current_ne*(m1/(l1+b1*current_strate)+r*((current_paccs+current_ns+current_ne+current_pacce)/K)),current_ne*(c21*m1/(l1+b1*current_strate)+g),0,0,current_pacce*c12]

                rate_sum = sum(rates)

                if sum(rs)<=0:
                        es = True
                if sum(re)<=0:
                        ee = True
                if rate_sum <= 0:
                        break

                tau = np.random.exponential(scale=1/rate_sum)

                t.append(t[-1] + tau)

                rand = random.uniform(0,1)

                #Normal cell division event
                if rand * rate_sum <= rates[0]:
                        paccs.append(paccs[-1])
                        normale.append(normale[-1])
                        pacce.append(pacce[-1])
                        strate.append(strate[-1])
                        if random.uniform(0,1)>mut:
                                strats.append(strats[-1])
                                nstrats.append(nstrats[-1])
                                normals.append(normals[-1] + 1)
                        else:
                                new_strat = strats[-1]+np.random.normal(0,0.01)
                                if drug_death(new_strat,m1)<drug_death(strats[-1],m1):
                                        strats.append(new_strat)
                                        nstrats.append(new_strat)
                                        normals.append(normals[-1] + 1)
                                elif drug_death(new_strat,m1)==drug_death(strats[-1],m1):
                                        strats.append(strats[-1])
                                        nstrats.append(nstrats[-1])
                                        normals.append(normals[-1] + 1)
                                else:
                                        strats.append(strats[-1])
                                        nstrats.append(nstrats[-1])
                                        normals.append(normals[-1])

                #Normal cell death event
                elif rand * rate_sum > rates[0] and rand * rate_sum <= sum(rates[:2]):
                        normals.append(normals[-1] - 1)
                        paccs.append(paccs[-1])
                        strats.append(strats[-1])
                        nstrats.append(nstrats[-1])
                        normale.append(normale[-1])
                        pacce.append(pacce[-1])
                        strate.append(strate[-1])

                #Normal switch event
                elif rand * rate_sum > sum(rates[:2]) and rand * rate_sum <= sum(rates[:3]):
                        normals.append(normals[-1]-1)
                        normale.append(normale[-1])
                        pacce.append(pacce[-1])
                        strate.append(strate[-1])
                        if random.uniform(0,1)<rho:
                                paccs.append(paccs[-1] + 1)
                        else:
                                paccs.append(paccs[-1])

                        strats.append(strats[-1])
                        nstrats.append(nstrats[-1])

                #PACC division event
                elif rand * rate_sum > sum(rates[:3]) and rand * rate_sum <= sum(rates[:4]):
                        normals.append(normals[-1])
                        paccs.append(paccs[-1] + 1)
                        strats.append(strats[-1])
                        nstrats.append(nstrats[-1])
                        normale.append(normale[-1])
                        pacce.append(pacce[-1])
                        strate.append(strate[-1])

                #PACC death event
                elif rand * rate_sum > sum(rates[:4]) and rand * rate_sum <= sum(rates[:5]):
                        normals.append(normals[-1])
                        paccs.append(paccs[-1]-1)
                        strats.append(strats[-1])
                        nstrats.append(nstrats[-1])
                        normale.append(normale[-1])
                        pacce.append(pacce[-1])
                        strate.append(strate[-1])

                #PACC switch event
                elif rand * rate_sum > sum(rates[:5]) and rand * rate_sum <= sum(rates[:6]):
                        normale.append(normale[-1])
                        pacce.append(pacce[-1])
                        strate.append(strate[-1])
                        if random.uniform(0,1)>mut: #using same mutation rate, but larger breadth
                                nstrats.append(nstrats[-1])
                        else:
                                new_strat = nstrats[-1]+np.random.normal(0,0.05)
                                if drug_death(new_strat,m1)<drug_death(nstrats[-1],m1):
                                        nstrats.append(new_strat)
                                else:
                                        nstrats.append(nstrats[-1])

                        if drug_death(nstrats[-1],m1)>dd_thres:
                                strats.append(strats[-1])
                                normals.append(normals[-1])
                                paccs.append(paccs[-1])
                        else:
                                strats.append(nstrats[-1])
                                normals.append(normals[-1]+2)
                                paccs.append(paccs[-1]-1)

                elif rand * rate_sum > sum(rates[:6]) and rand * rate_sum <= sum(rates[:7]):
                        normals.append(normals[-1])
                        paccs.append(paccs[-1])
                        strats.append(strats[-1])
                        pacce.append(pacce[-1])
                        if random.uniform(0, 1) > mut:
                                strate.append(strate[-1])
                                normale.append(normale[-1] + 1)
                        else:
                                new_s = strate[-1] + np.random.normal(0, 0.01)
                                if drug_death(new_s, m1) < drug_death(strate[-1], m1):
                                        strate.append(new_s)
                                        normale.append(normale[-1] + 1)
                                elif drug_death(new_s, m1) == drug_death(strate[-1], m1):
                                        strate.append(strate[-1])
                                        normale.append(normale[-1] + 1)
                                else:
                                        strate.append(strate[-1])
                                        normale.append(normale[-1])

                elif rand * rate_sum > sum(rates[:7]) and rand * rate_sum <= sum(rates[:8]):
                        normals.append(normals[-1])
                        paccs.append(paccs[-1])
                        strats.append(strats[-1])
                        normale.append(normale[-1] - 1)
                        pacce.append(pacce[-1])
                        strate.append(strate[-1])

                elif rand * rate_sum > sum(rates[:8]) and rand * rate_sum <= sum(rates[:9]):
                        normals.append(normals[-1])
                        paccs.append(paccs[-1])
                        strats.append(strats[-1])
                        normale.append(normale[-1] - 1)
                        if random.uniform(0, 1) < rho:
                                pacce.append(pacce[-1] + 1)
                        else:
                                pacce.append(pacce[-1])

                        strate.append(strate[-1])

                elif rand * rate_sum > sum(rates[:9]) and rand * rate_sum <= sum(rates[:10]):
                        normals.append(normals[-1])
                        paccs.append(paccs[-1])
                        strats.append(strats[-1])
                        normale.append(normale[-1])
                        pacce.append(pacce[-1] + 1)
                        strate.append(strate[-1])

                elif rand * rate_sum > sum(rates[:10]) and rand * rate_sum <= sum(rates[:11]):
                        normals.append(normals[-1])
                        paccs.append(paccs[-1])
                        strats.append(strats[-1])
                        normale.append(normale[-1])
                        pacce.append(pacce[-1] - 1)
                        strate.append(strate[-1])

                elif rand * rate_sum > sum(rates[:11]) and rand * rate_sum <= sum(rates[:12]):
                        normals.append(normals[-1])
                        paccs.append(paccs[-1])
                        strats.append(strats[-1])
                        pacce.append(pacce[-1] - 1)
                        if random.uniform(0, 1) > mut:  # using same mutation rate, but larger breadth
                                strate.append(strate[-1])
                                normale.append(normale[-1] + 2)
                        else:
                                new_s = strate[-1] + np.random.normal(0, 0.05)
                                if drug_death(new_s, m1) < drug_death(strate[-1], m1):
                                        strate.append(new_s)
                                        normale.append(normale[-1] + 2)
                                elif drug_death(new_s, m1) == drug_death(strate[-1], m1):
                                        strate.append(strate[-1])
                                        normale.append(normale[-1] + 2)
                                else:
                                        strate.append(strate[-1])
                                        normale.append(normale[-1])

        norm_macros.append(normals)
        pacc_macros.append(paccs)
        strat_macros.append(strats)
        nstrat_macros.append(nstrats)

        norm_macroe.append(normale)
        pacc_macroe.append(pacce)
        strat_macroe.append(strate)

        t_macro.append(t)

        if ee == True:
                extincte.append(1)

        if es == True:
                extincts.append(1)

for i in range(num_sims):
        run_this()

lwi = .3

extinct_eventse = len(extincte)
extinct_eventss = len(extincts)
#avg_extinct = np.mean(extinct)
#std_extinct = np.std(extinct)

print('Evolutionary Tracking Extinctions: ' + str(extinct_eventse))
print('Self-Genetic Modification Extinctions: ' + str(extinct_eventss))
#print(avg_extinct)
#print(std_extinct)


plt.figure()
plt.subplot(211)
plt.title('Self-Genetic Modification vs Evolutionary Tracking: Low Dose')
plt.plot(t_macro[0], norm_macros[0], label='2N+ SGM', c='r', lw=lwi)
plt.plot(t_macro[0], pacc_macros[0], label='PACC SGM', c='b', lw=lwi)
plt.plot(t_macro[0], norm_macroe[0], label='2N+ ET', c='y', lw=lwi)
plt.plot(t_macro[0], pacc_macroe[0], label='PACC ET', c='mediumslateblue', lw=lwi)

'''
for i in range(1,num_sims):
    plt.plot(t_macro[i],norm_macros[i], c='r',lw=lwi)
    plt.plot(t_macro[i],pacc_macros[i], c='b',lw=lwi)
    plt.plot(t_macro[i],norm_macroe[i], c='y',lw=lwi)
    plt.plot(t_macro[i],pacc_macroe[i], c='mediumslateblue',lw=lwi)
plt.legend()
plt.grid(True)
ax = plt.gca()
ax.axvspan(0, 600, facecolor=(0.5,0.5,0.5))
ax.axvspan(100, 200, facecolor=(0.5,0.5,0.5))
ax.axvspan(300, 400, facecolor=(0.5,0.5,0.5))
ax.axvspan(500, 600, facecolor=(0.5,0.5,0.5))
ax.axvspan(700, 800, facecolor=(0.5,0.5,0.5))
ax.axvspan(900, 1000, facecolor=(0.5,0.5,0.5))
plt.ylabel('Population Size')
plt.subplot(212)
plt.plot(t_macro[0],strat_macros[0], label='Drug Resist SGM',c='k',lw=lwi)
plt.plot(t_macro[0],strat_macroe[0], label='Drug Resist ET',c='m',lw=lwi)
plt.legend()
for i in range(1, num_sims):
    plt.plot(t_macro[i],strat_macros[i], c='k',lw=lwi)
    plt.plot(t_macro[i],strat_macroe[i], c='m',lw=lwi)
plt.grid(True)
ax = plt.gca()
ax.axvspan(0, 600, facecolor=(0.5,0.5,0.5))
ax.axvspan(100, 200, facecolor=(0.5,0.5,0.5))
ax.axvspan(300, 400, facecolor=(0.5,0.5,0.5))
ax.axvspan(500, 600, facecolor=(0.5,0.5,0.5))
ax.axvspan(700, 800, facecolor=(0.5,0.5,0.5))
ax.axvspan(900, 1000, facecolor=(0.5,0.5,0.5))
plt.ylabel('Drug Resistance')
plt.xlabel('Time')
plt.tight_layout()
plt.show()
'''