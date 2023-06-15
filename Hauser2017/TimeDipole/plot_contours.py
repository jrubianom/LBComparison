import numpy as np
import matplotlib.pyplot as plt

i_arr, t = np.loadtxt('fileTime.dat', unpack=True)

refinements = np.unique(i_arr)

for i in refinements:
    stringB = 'fileB_' + str(int(i)) + '.dat'
    stringE = 'fileE_' + str(int(i)) + '.dat'
    stringContour = 'contour_' + str(int(i)) + '.dat'

    dataB = np.loadtxt(stringB, unpack=True)
    dataE = np.loadtxt(stringB, unpack=True)
    dataContour = np.loadtxt(stringContour, unpack=True)

    fig1, (axB, axE) = plt.subplots(ncols=1, nrows=2, sharex='col')

    axB.set_xlabel('r')
    axB.set_ylabel('B')
    axB.set_title('$B(r, \\theta=\\pi/2)$')

    axE.set_xlabel('r')
    axE.set_ylabel('E')
    axE.set_title('$E(r, \\theta=\\pi/2)$')

    axB.plot(dataB[0], dataB[2], 'k-', label='Theoretical')
    axB.plot(dataB[0][0:-1:5], dataB[1][0:-1:5], 'kx', label='Simulated')

    axE.plot(dataE[0], dataE[2], 'k-', label='Theoretical')
    axE.plot(dataE[0][0:-1:5], dataE[1][0:-1:5], 'kx', label='Simulated')

    plt.savefig('Fields'+str(int(i))+'.jpg')

    fig2, axContour = plt.subplots()

    axContour.set_xlabel('x')
    axContour.set_ylabel('z')
    
    N = int(100+10*(i-1))
    axContour.contourf(dataContour[0].reshape((N,N)), dataContour[1].reshape((N,N)), np.abs(dataContour[2].reshape((N,N))), levels=25)

    plt.savefig('Contour.jpg')

    plt.show()