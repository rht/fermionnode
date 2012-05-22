#from __future__ import division
from pyqm import *
from scipy.misc import factorial
from scipy.optimize import fsolve, broyden1, broyden2, newton_krylov, newton
from scipy.constants import hbar
from itertools import permutations
import time

# QHE
# http://www.itp.phys.ethz.ch/education/lectures_fs11/solid/QHE.pdf
# http://online.physics.uiuc.edu/courses/phys598PTD/fall09/L16.pdf

#####
# QHO
#####
def hermite(n, x): # return H n(x)
    if n < 0:
        return 0.0
    elif n == 0:
        return 1.0
    else:
        return 2*x*hermite(n-1, x) - 2*(n-1)*hermite(n-2, x)


def QHO_psi_n(n, x, m, omega):
    """1D QHO eigenfunctions"""
    #http://en.wikipedia.org/wiki/Quantum_harmonic_oscillator
    #nu = m * omega / hbar
    nu = m * omega
    # normalization coefficient
    C =  (nu/pi)**(1/4) * sqrt(1/(2**n*factorial(n)))
    return C * exp(-nu* x**2 /2) * hermite(n, sqrt(nu)*x)


#eigenstate of electron in a 2D periodic box
def phi2D(n, length):
    def result(*x):
        #PBC
        return (sqrt(2 / length) *
                    cos(n[0]*2*float64(pi)*x[0]/length) *
                    cos(n[1]*2*float64(pi)*x[1]/length))
    return result



def E(n):
    #todo: make this 2D
    return n * n / (8 * length ** 2)


def antisymmetrize(Psis):
    def psi(xs):
        slatermatrix = matrix(array([[Psi(*x) for Psi in Psis]
                                     for x in xs]))
        slaterlogdet = slogdet(slatermatrix)
        #print(str(slaterlogdet) + " "+ str(slaterlogdet[0] * exp(slaterlogdet[1]) ))
        return slaterlogdet[0] * exp(slaterlogdet[1]) / sqrt(factorial(len(xs))), slatermatrix
    return psi


def wfgridcreator(wavefunction, X, Y):
    '''return grid of a "complicated" function'''
    meshsize = len(X)
    wfgrid = zeros(shape=shape(X))
    for i in range(meshsize):
        for j in range(meshsize):
            wfgrid[i, j] = wavefunction(X[i,j],Y[i,j])
    return wfgrid



### node search
def bruteforcezerofinder(functiongrid,length,meshsize,X,Y):
    zerolocations = []
    R = linspace(0,length,num=meshsize)
    mean_wfgrid = mean(abs(functiongrid))
    print(mean_wfgrid)
    for i in range(meshsize):
        for j in range(meshsize):
            #if abs(wavefunction([[i,j]] + otherelectrons)) < mean_wfgrid/100:
            if abs(functiongrid[i,j]) < mean_wfgrid/20:
                zerolocations.append([X[i,j],Y[i,j]])
    return array(zerolocations).T # array of x, array of y

###NOW, for finding the nodes
#root finding all star can be found here http://docs.scipy.org/doc/scipy-0.10.1/reference/optimize.html
##now, find nodes using fsolve or broyden1 or broyden2(broyden's "bad" method)
#wavefunctionwith4variablesfixed = lambda x: [float(wavefunction([x[0],x[1],5,2,4,5])),0]
#nodes = fsolve(wavefunctionwith4variablesfixed,array([.1,.2]))[0]
#print nodes
#plot(*nodes)


