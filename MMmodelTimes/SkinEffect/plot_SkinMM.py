import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

i_arr, times = np.loadtxt('times.dat', unpack=True, skiprows=1)

refinements, counts = np.unique(i_arr, return_counts=True)

y = np.zeros_like(refinements)
std_dev = np.zeros_like(y)
x = np.zeros_like(y)

for i in refinements:
    iter = int(i)
    try:
        string = 'data_'+str(iter)+'.dat'
        data = np.loadtxt(string, unpack=True)
    except FileNotFoundError:
        print('broken')
        break
    idxs = data[0]>=0.25
    # Simulation data
    exp = data[1][idxs]
    # Theory
    teo = data[2][idxs]
    # Compute cost function
    cost = np.sum((exp-teo)**2/teo**2)/teo.size
    x[iter-1] = cost
    # Average time and compute error
    y[iter-1] = np.mean(times, where=i_arr==i)
    std_dev[iter-1] = np.std(times, where=i_arr==i)

std_err = std_dev/np.sqrt(counts[0])

fit = linregress(np.log(x), np.log(y))
fitted = lambda x: x**(fit.slope)*np.exp(fit.intercept)

print('LB MM\nSkin Effect')
print('[power, coeff, r] = ', fit.slope, fit.intercept, fit.rvalue)

fig1, ax1 = plt.subplots(ncols=1, nrows=1)
ax1.set_xlabel('Cost')
ax1.set_ylabel('Time [s]')
ax1.set_title('Skin Effect (MM)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.errorbar(x, y, yerr=std_err, fmt='k.')
string = '{0:.3e}'.format(np.exp(fit.intercept)) + '$\\cdot\\epsilon' + '^{' + '{0:.3f}'.format(fit.slope) + '}$'
ax1.plot(x, fitted(x), 'k-', label=string)
ax1.legend()

plt.savefig('CPU_time_vs_cost_SE_MM.jpg')

#plt.show()