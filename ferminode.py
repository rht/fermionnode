from __future__ import division
from pylab import *
from scipy.misc import factorial
from scipy.optimize import fsolve, broyden1, broyden2, newton_krylov, newton
from itertools import permutations
import mpl_toolkits.mplot3d.axes3d as plot3
'''use this http://code.enthought.com/chaco/'''
#random cool paper http://www.sciencemag.org/content/318/5852/949.full
#http://physicsandphysicists.blogspot.com/2007/11/simplest-2-slit-and-decoherence.html

meshsize = 100 
length = 2
higheststate = 3
numelectrons = int(factorial(higheststate))  # number of permutations

#initialize coordinates
X = linspace(0,length,num=meshsize)
Y = X
#K = linspace(-pi,pi,num=meshsize)

#eigenstate of electron in a 2D box
def phi(n):
    #@vectorize
    def result(x):
        #PBC
        return sqrt(2./length) * cos(n[0]*2*pi*x[0]/length) * cos(n[1]*2*pi*x[1]/length)
    return result
def E(n):
    #todo: make this 2D
    return n * n / (8* length**2)

def randpos(): #random position chooser
    return length * random(2)

def antisymmetrize(Psis):
    def psi(x):
        N = len(Psis)
        slatermatrix = matrix(array([[Psis[i](x[j]) for i in range(N)] for j in range(N)]  ))
        return det(slatermatrix) / sqrt(factorial(N))
    return psi

fig=figure()
ax = plot3.Axes3D(fig)
X, Y = meshgrid(X,X)


#wavefunction = antisymmetrize([phi([1,1]),phi([1,2]),phi([2,1]),phi([2,2]),phi([1,3]),phi([3,1])])
wavefunction = antisymmetrize([phi(i) for i in permutations(range(1,higheststate+1))])
wfgrid = zeros(shape=shape(X))
otherelectrons = [randpos() for i in range(numelectrons-1)] # number of permutation -1
for i in range(meshsize):
    for j in range(meshsize):
        wfgrid[i,j] = wavefunction([[X[i,j],Y[i,j]]] + otherelectrons)

ax.plot_surface(X,Y,wfgrid)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('$\psi$')
savefig("plots/%delectrons-%dmeshsize-%dlength.png"%(numelectrons,meshsize,length))



###NOW, for finding the nodes
#root finding all star can be found here http://docs.scipy.org/doc/scipy-0.10.1/reference/optimize.html
##now, find nodes using fsolve or broyden1 or broyden2(broyden's "bad" method)
#wavefunctionwith4variablesfixed = lambda x: [float(wavefunction([x[0],x[1],5,2,4,5])),0]
#nodes = fsolve(wavefunctionwith4variablesfixed,array([.1,.2]))[0]
#print nodes
#plot(*nodes)

zerolocations = []
R = linspace(0,length,num=meshsize)
mean_wfgrid = mean(abs(wfgrid))
for i in R:
    for j in R:
        #print abs(wavefunction([[i,j]] + otherelectrons))
        if abs(wavefunction([[i,j]] + otherelectrons)) < mean_wfgrid/100:
            zerolocations.append([i,j])

figure()
zerox,zeroy = array(zerolocations).T
plot(zerox,zeroy,'bo')
otherelectronsplot = array(otherelectrons).T
plot(otherelectronsplot[0],otherelectronsplot[1],'go')
xlim(0,2)
ylim(0,2)
savefig("plots/node-%delectrons-%dmeshsize-%dlength.png"%(numelectrons,meshsize,length))
savetxt("plots/wfgrid-%delectrons-%dmeshsize-%dlength.txt"%(numelectrons,meshsize,length), wfgrid)

show()

