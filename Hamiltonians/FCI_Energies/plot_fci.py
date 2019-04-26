import numpy as np
import matplotlib.pyplot as plt


fontsize = 12
plt.rcParams["font.size"] = fontsize
plt.rcParams["axes.labelsize"] = fontsize
plt.rcParams["axes.titlesize"] = fontsize
# 'serif' | 'sans-serif' | 'cursive' | 'fantasy' | 'monospace'
fontfamily = 'sans-serif'
plt.rcParams["font.family"] = fontfamily

fci_sto3g = np.genfromtxt('H2_sto-3g_FCI_energy.txt')
fci_631g = np.genfromtxt('H2_6-31g_FCI_energy.txt')

fig, ax = plt.subplots(figsize=[8,8])

ax.plot(fci_sto3g[:,0],fci_sto3g[:,1],c='blue',label='sto-3g')
ax.plot(fci_631g[:,0],fci_631g[:,1],c='green',label='6-31g')

ax.set_xlabel('H-H separation (angstrom)')
ax.set_ylabel('FCI energy (hartrees)')
ax.set_title('H2 FCI PES in two bases')
ax.legend()

plt.savefig('fci_energy_sto-3g_6-31g.png')
#plt.show()
plt.close()
