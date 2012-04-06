from __future__ import division
from pylab import *
from scipy.misc import factorial
from scipy.optimize import fsolve, broyden1, broyden2
import mpl_toolkits.mplot3d.axes3d as plot3
'''use this http://code.enthought.com/chaco/'''
#random cool paper http://www.sciencemag.org/content/318/5852/949.full
#http://physicsandphysicists.blogspot.com/2007/11/simplest-2-slit-and-decoherence.html

meshsize = 50 
length = 2

#initialize coordinates
X = linspace(0,length,num=meshsize)
#K = linspace(-pi,pi,num=meshsize)

#eigenstate of electron in a box
def phi(n):
    @vectorize
    def result(x):
        return sqrt(2./length) * sin(n*pi*x/length)
    return result
def E(n):
    return n * n / (8* length**2)

def antisymmetrize(Psis):
    def psi(x):
        N = len(Psis)
        slatermatrix = matrix(array([[Psis[i](x[j]) for i in range(N)] for j in range(N)]  ))
        return det(slatermatrix) / sqrt(factorial(N))
    return psi

fig=figure()
ax = plot3.Axes3D(fig)
X, Y = meshgrid(X,X)
wavefunction = antisymmetrize([phi(i) for i in range(1,5)])
wfgrid = zeros(shape=shape(X))
for i in range(meshsize):
    for j in range(meshsize):
        wfgrid[i,j] = wavefunction([X[i,j],Y[i,j],5.,2.,4.,5.])

ax.plot_surface(X,Y,wfgrid)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('$\psi$')

#now, find nodes using fsolve or broyden1 or broyden2(broyden's "bad" method)
wavefunctionwith4variablesfixed = lambda x: wavefunction([x[0],x[1],5,2,4,5])
nodes = broyden1(wavefunctionwith4variablesfixed,[.1,.2])[0]
plot(*nodes)


#zerolocations = []
#for i in range(meshsize):
    #for j in range(meshsize):
        #if abs(wfgrid[i,j]) < 1e-25:
            #print wfgrid[i,j]
            #zerolocations.append([i,j])
            ##zerolocations.append((X[i,j],Y[i,j]))

#figure()
#plot(*array(zerolocations).T)

show()

