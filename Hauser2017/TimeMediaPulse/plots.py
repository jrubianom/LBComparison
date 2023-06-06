import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Errors.dat", unpack=True, skiprows=1)

refinements = np.unique(data[0])

error_ref = np.zeros_like(refinements)
error_trans = np.zeros_like(refinements)

time = np.zeros_like(refinements)
std_dev = np.zeros_like(refinements)

for i in range(len(refinements)):
    error_ref[i] = np.mean(data[4], where=data[0]==refinements[i])
    error_trans[i] = np.mean(data[7], where=data[0]==refinements[i])

    time[i] = np.mean(data[1], where=data[0]==refinements[i])
    std_dev[i] = np.std(data[1], where=data[0]==refinements[i])

std_err = std_dev/np.sqrt(len(data[0]/len(refinements)))

# Figsize based on PRL column width and the Golden Ratio
fig1, ax1 = plt.subplots(ncols=1, nrows=1, figsize=(3.375, 2.086), dpi=500)

ax1.set_xlabel('Relative error')
ax1.set_ylabel('CPU time [s]')
ax1.set_title('Reflected pulse')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.errorbar(error_ref, time, yerr=std_dev, fmt='k.')

plt.savefig('CPU_time_vs_rel_error_ref.jpg')

fig2, ax2 = plt.subplots(ncols=1, nrows=1, figsize=(3.375, 2.086), dpi=500)

ax2.set_xlabel('Relative error')
ax2.set_ylabel('CPU time [s]')
ax2.set_title('Transmited pulse')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.errorbar(error_trans, time, yerr=std_dev, fmt='k.')

plt.savefig('CPU_time_vs_rel_error_trans.jpg')

plt.show()