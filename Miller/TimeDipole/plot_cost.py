import numpy as np
from scipy.stats import linregress
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

fitB = linregress(np.log(xB), np.log(y))
fittedB = lambda x: x**(fitB.slope)*np.exp(fitB.intercept)

fitE = linregress(np.log(xE), np.log(y))
fittedE = lambda x: x**(fitE.slope)*np.exp(fitE.intercept)

print('LB HV\nDipole')
print('Magnetic field [power, coeff, r] = ', fitB.slope, fitB.intercept, fitB.rvalue)

print('LB HV\nDipole')
print('Electric field [power, coeff, r] = ', fitE.slope, fitE.intercept, fitE.rvalue)

fig1, (axB, axE) = plt.subplots(ncols=2, nrows=1, sharey='col')

axB.set_xlabel('Cost')
axB.set_ylabel('Time [s]')
axB.set_title('Magnetic field error (HV)')
axB.set_xscale('log')
axB.set_yscale('log')
axB.errorbar(xB, y, yerr=std_err, fmt='k.')
string = '{0:.3e}'.format(np.exp(fitB.intercept)) + '$\\cdot\\epsilon' + '^{' + '{0:.3f}'.format(fitB.slope) + '}$'
axB.plot(xB, fittedB(xB), 'k-', label=string)
axB.legend()

axE.set_xlabel('Cost')
axE.set_ylabel('Time [s]')
axE.set_title('Electric field error (HV)')
axE.set_xscale('log')
axE.set_yscale('log')
axE.errorbar(xE, y, yerr=std_err, fmt='k.')
string = '{0:.3e}'.format(np.exp(fitE.intercept)) + '$\\cdot\\epsilon' + '^{' + '{0:.3f}'.format(fitE.slope) + '}$'
axE.plot(xE, fittedE(xE), 'k-', label=string)
axE.legend()

plt.savefig('CPU_time_vs_cost_DIP_HV.jpg')

plt.show()