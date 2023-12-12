import numpy as np
import matplotlib.pyplot as plt
import random

num_sims = 10

dd_thres = 0.2 #0.2, 0.4, 0.6
mut = 0.05
dose = 1 #.6, .8, 1

g = 0.02
K=100
r = 0.6
l1 = 1
b1 = 1
c21 = 0.7
c12 = 0.2
rho = 0.7

tend = 1000

norm_macro = []
pacc_macro = []
strat_macro = []
nstrat_macro = []
t_macro = []
extinct = []

def drug_death(current_strat,m1):
        return m1/(l1+b1*current_strat)

def run_this():
        normal = [10]
        pacc = [0]
        strat = [0]
        nstrat = [0]
        t = [0]

        while t[-1] < tend:
                if t[-1] > 200 and t[-1] < 800:
                        m1 = dose
                else:
                        m1 = 0

                '''if t[-1] > 900:
                        m1 = dose
                elif t[-1] > 800:
                        m1 = 0
                elif t[-1] > 700:
                        m1 = dose
                elif t[-1] > 600:
                        m1 = 0
                elif t[-1] > 500:
                        m1 = dose
                elif t[-1] > 400:
                        m1 = 0
                elif t[-1] > 300:
                        m1 = dose
                elif t[-1] > 200:
                        m1 = 0
                elif t[-1] > 100:
                        m1 = dose
                else:
                        m1 = 0'''

                current_n = normal[-1]
                current_pacc = pacc[-1]
                current_strat = strat[-1]

                #broken into birth, death, switching

                rates = [(current_n*r), current_n*(m1/(l1+b1*current_strat)+r*((current_pacc+current_n)/K)),current_n*(c21*m1/(l1+b1*current_strat)+g),0,0,current_pacc*c12]

                rate_sum = sum(rates)

                if rate_sum == 0:
                        extinct.append(t[-1])
                        break

                tau = np.random.exponential(scale=1/rate_sum)

                t.append(t[-1] + tau)

                rand = random.uniform(0,1)

                #Normal cell division event
                if rand * rate_sum <= rates[0]:
                        pacc.append(pacc[-1])
                        if random.uniform(0,1)>mut:
                                strat.append(strat[-1])
                                nstrat.append(nstrat[-1])
                                normal.append(normal[-1] + 1)
                        else:
                                new_strat = strat[-1]+np.random.normal(0,0.01)
                                if drug_death(new_strat,m1)<drug_death(strat[-1],m1):
                                        strat.append(new_strat)
                                        nstrat.append(new_strat)
                                        normal.append(normal[-1] + 1)
                                elif drug_death(new_strat,m1)==drug_death(strat[-1],m1):
                                        strat.append(strat[-1])
                                        nstrat.append(nstrat[-1])
                                        normal.append(normal[-1] + 1)
                                else:
                                        strat.append(strat[-1])
                                        nstrat.append(nstrat[-1])
                                        normal.append(normal[-1])

                #Normal cell death event
                elif rand * rate_sum > rates[0] and rand * rate_sum <= sum(rates[:2]):
                        normal.append(normal[-1] - 1)
                        pacc.append(pacc[-1])
                        strat.append(strat[-1])
                        nstrat.append(nstrat[-1])

                #Normal switch event
                elif rand * rate_sum > sum(rates[:2]) and rand * rate_sum <= sum(rates[:3]):
                        normal.append(normal[-1]-1)
                        if random.uniform(0,1)<rho:
                                pacc.append(pacc[-1] + 1)
                        else:
                                pacc.append(pacc[-1])

                        strat.append(strat[-1])
                        nstrat.append(nstrat[-1])

                #PACC division event
                elif rand * rate_sum > sum(rates[:3]) and rand * rate_sum <= sum(rates[:4]):
                        normal.append(normal[-1])
                        pacc.append(pacc[-1] + 1)
                        strat.append(strat[-1])
                        nstrat.append(nstrat[-1])

                #PACC death event
                elif rand * rate_sum > sum(rates[:4]) and rand * rate_sum <= sum(rates[:5]):
                        normal.append(normal[-1])
                        pacc.append(pacc[-1]-1)
                        strat.append(strat[-1])
                        nstrat.append(nstrat[-1])

                #PACC switch event
                elif rand * rate_sum > sum(rates[:5]) and rand * rate_sum <= sum(rates[:6]):
                        if random.uniform(0,1)>mut: #using same mutation rate, but larger breadth
                                nstrat.append(nstrat[-1])
                        else:
                                new_strat = nstrat[-1]+np.random.normal(0,0.05)
                                if drug_death(new_strat,m1)<drug_death(nstrat[-1],m1):
                                        nstrat.append(new_strat)
                                else:
                                        nstrat.append(nstrat[-1])

                        if drug_death(nstrat[-1],m1)>dd_thres:
                                strat.append(strat[-1])
                                normal.append(normal[-1])
                                pacc.append(pacc[-1])
                        else:
                                strat.append(nstrat[-1])
                                normal.append(normal[-1]+2)
                                pacc.append(pacc[-1]-1)

        norm_macro.append(normal)
        pacc_macro.append(pacc)
        strat_macro.append(strat)
        nstrat_macro.append(nstrat)
        t_macro.append(t)

for i in range(num_sims):
        run_this()

lwi = 0.3

extinct_events = len(extinct)
avg_extinct = np.mean(extinct)
std_extinct = np.std(extinct)

print(extinct_events)
print(avg_extinct)
print(std_extinct)


plt.figure()
plt.subplot(211)
plt.title('Minor Self-Genetic Modification')
plt.plot(t_macro[0],norm_macro[0], label='2N+', c='r',lw=lwi)
plt.plot(t_macro[0],pacc_macro[0], label='PACC', c='b',lw=lwi)
for i in range(1,num_sims):
    plt.plot(t_macro[i],norm_macro[i], c='r',lw=lwi)
    plt.plot(t_macro[i],pacc_macro[i], c='b',lw=lwi)
plt.legend()
plt.grid(True)
ax = plt.gca()
ax.axvspan(200, 800, facecolor=(0.5,0.5,0.5))
'''ax.axvspan(100, 200, facecolor=(0.5,0.5,0.5))
ax.axvspan(300, 400, facecolor=(0.5,0.5,0.5))
ax.axvspan(500, 600, facecolor=(0.5,0.5,0.5))
ax.axvspan(700, 800, facecolor=(0.5,0.5,0.5))
ax.axvspan(900, 1000, facecolor=(0.5,0.5,0.5))'''
plt.ylabel('Population Size')
plt.subplot(212)
plt.plot(t_macro[0],strat_macro[0], c='k',lw=lwi)
for i in range(1, num_sims):
    plt.plot(t_macro[i],strat_macro[i], c='k',lw=lwi)
plt.grid(True)
ax = plt.gca()
ax.axvspan(200, 800, facecolor=(0.5,0.5,0.5))
'''ax.axvspan(100, 200, facecolor=(0.5,0.5,0.5))
ax.axvspan(300, 400, facecolor=(0.5,0.5,0.5))
ax.axvspan(500, 600, facecolor=(0.5,0.5,0.5))
ax.axvspan(700, 800, facecolor=(0.5,0.5,0.5))
ax.axvspan(900, 1000, facecolor=(0.5,0.5,0.5))'''
plt.ylabel('Drug Resistance')
plt.xlabel('Time')
plt.tight_layout()
plt.show()
