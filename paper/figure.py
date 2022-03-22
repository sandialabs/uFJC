import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from ufjc import uFJC


num_points = 250
plt.rc('font', size=10)
mpl.rcParams['figure.dpi'] = 300
plt.figure()
plt.subplot(2, 1, 1)

for varepsilon in [10, 25, 50, 100]:
    model = uFJC(potential='morse', N_b=8, varepsilon=varepsilon)
    eta = np.linspace(0, model.eta_max, num_points)
    plt.plot(model.gamma(eta), eta, label=r'$\varepsilon=$'+str(varepsilon))
    plt.plot(model.gamma(model.eta_max), model.eta_max, 'kx')
plt.xlim([0, 1.7])
plt.ylim([0, 26])
plt.xlabel(r'$\gamma(\eta)$')
plt.ylabel(r'$\eta$')
plt.legend()

plt.subplot(2, 1, 2)
for N_b in [3, 5, 10, 25]:
    model = uFJC(potential='morse', N_b=N_b)
    gamma = np.linspace(0, 1.6, num_points)
    plt.plot(gamma, model.nondim_g_eq(gamma), label=r'$N_b=$'+str(N_b))
plt.xlim([0, 1.1])
plt.ylim([0, 5.2])
plt.xlabel(r'$\gamma$')
plt.ylabel(r'$\mathscr{g}_\mathrm{eq}(\gamma)$')
plt.legend()

plt.tight_layout()
plt.savefig('paper/figure.png')
