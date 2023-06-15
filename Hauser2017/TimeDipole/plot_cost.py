import numpy as np
import matplotlib.pyplot as plt

i_arr, times = np.loadtxt('fileTime.dat', unpack=True)

refinements, counts = np.unique(i_arr, return_counts=True)

y = np.zeros_like(refinements)
std_dev = np.zeros_like(y)
xB = np.zeros_like(y)
xE = np.zeros_like(y)

for i in refinements:
    iter = int(i)
    L = 100+10*(iter-1)

    try:
        dataB = np.loadtxt('fileB_'+str(iter)+'.dat', unpack=True)
        dataE = np.loadtxt('fileE_'+str(iter)+'.dat', unpack=True)
    except FileNotFoundError:
        print('broken')
        break

    idxs = dataB[0]>=L/2+L/10
    # Simulation data
    expB = dataB[1][idxs]
    expE = dataE[1][idxs]
    # Theory
    teoB = dataB[2][idxs]
    teoE = dataE[2][idxs]
    # Compute cost function
    costB = np.sum((expB-teoB)**2)/teoB.size
    costE = np.sum((expE-teoE)**2)/teoE.size
    xB[iter-1] = costB
    xE[iter-1] = costE
    # Average time and compute error
    y[iter-1] = np.mean(times, where=i_arr==i)
    std_dev[iter-1] = np.std(times, where=i_arr==i)

std_err = std_dev/np.sqrt(counts[0])

# Figsize based on PRL column width and the Golden Ratio
fig1, (axB, axE) = plt.subplots(ncols=2, nrows=1, sharey='col')

axB.set_xlabel('Cost')
axB.set_ylabel('Time [s]')
axB.set_title('Dipole')
axB.set_xscale('log')
axB.set_yscale('log')
axB.errorbar(xB, y, yerr=std_err, fmt='k.')

axE.set_xlabel('Cost')
axE.set_ylabel('Time [s]')
axE.set_title('Dipole')
axE.set_xscale('log')
axE.set_yscale('log')
axE.errorbar(xE, y, yerr=std_err, fmt='k.')

plt.savefig('CPU_time_vs_cost_DIP.jpg')