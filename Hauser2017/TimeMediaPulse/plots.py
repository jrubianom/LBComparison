import numpy as np
from scipy.stats import linregress
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

fit_ref = linregress(np.log(error_ref), np.log(time))
fitted_ref = lambda e: e**fit_ref.slope*np.exp(fit_ref.intercept)

fit_trans = linregress(np.log(error_trans), np.log(time))
fitted_trans = lambda e: e**fit_trans.slope*np.exp(fit_trans.intercept)

print('LB MM\nGaussian pulse crossing an interface')
print('Time vs Reflected pulse error', '[power, coeff] = ', fit_ref.slope, fit_ref.intercept, fit_ref.rvalue)
print('Time vs Transmitted pulse error', '[power, coeff] = ', fit_trans.slope, fit_trans.intercept, fit_trans.rvalue)

fig1, ax1 = plt.subplots()

ax1.set_xlabel('Relative error')
ax1.set_ylabel('CPU time [s]')
ax1.set_title('Reflected pulse (HV)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.errorbar(error_ref, time, yerr=std_dev, fmt='k.')
string_ref = '{0:.3f}'.format(np.exp(fit_ref.intercept)) + '$\\epsilon_{ref}' + '^{' + '{0:.3f}'.format(fit_ref.slope) + '}$'
ax1.plot(error_ref, fitted_ref(error_ref), 'k-', label=string_ref)
ax1.legend()

plt.savefig('CPU_time_vs_rel_error_ref_HV.jpg')

fig2, ax2 = plt.subplots()

ax2.set_xlabel('Relative error')
ax2.set_ylabel('CPU time [s]')
ax2.set_title('Transmited pulse (MM)')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.errorbar(error_trans, time, yerr=std_dev, fmt='k.')
string_trans = '{0:.3f}'.format(np.exp(fit_trans.intercept)) + '$\\epsilon_{trans}' + '^{' + '{0:.3f}'.format(fit_trans.slope) + '}$'
ax2.plot(error_trans, fitted_trans(error_trans), 'k-', label=string_trans)
ax2.legend()

plt.savefig('CPU_time_vs_rel_error_trans_HV.jpg')

plt.show()