from __future__ import division
from pyqm import *
from scipy.misc import factorial
from scipy.optimize import fsolve, broyden1, broyden2, newton_krylov, newton
from itertools import permutations
import time
'''use this http://code.enthought.com/chaco/'''
#random cool paper http://www.sciencemag.org/content/318/5852/949.full
#http://physicsandphysicists.blogspot.com/2007/11/simplest-2-slit-and-decoherence.html

starttime = time.time()
meshsize = 100 
length = 2
higheststate = 4
numelectrons = int(factorial(higheststate))  # number of permutations

#initialize coordinates
X = linspace(0,length,num=meshsize)
Y = X
#K = linspace(-pi,pi,num=meshsize)

#eigenstate of electron in a 2D box
def phi(n):
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

###############
#Plotting the wavefunction cross section
###############
X, Y = meshgrid(X,X)


#wavefunction = antisymmetrize([phi([1,1]),phi([1,2]),phi([2,1]),phi([2,2]),phi([1,3]),phi([3,1])])
wavefunction = antisymmetrize([phi(i) for i in permutations(range(1,higheststate+1))])
wfgrid = zeros(shape=shape(X))
otherelectrons = [randpos() for i in range(numelectrons-1)] # number of permutation -1
#It will be hard to compare if the position of otherelectrons always change, so
#24 electrons
#otherelectrons = [array([ 1.4156458 ,  0.61850838]), array([ 1.10460647,  0.67602137]), array([ 1.69724783,  1.05551967]), array([ 1.1111193,  0.4182782]), array([ 0.81117111,  1.86703291]), array([ 1.12617196,  0.53396474]), array([ 1.42082103,  1.95538821]), array([ 1.7028005 ,  1.82304221]), array([ 0.13782677,  1.74442594]), array([ 1.91003568,  1.833085  ]), array([ 1.97984843,  1.97071588]), array([ 1.70103109,  0.8389556 ]), array([ 0.6029357 ,  1.24555397]), array([ 0.01928903,  0.51623989]), array([ 1.95794243,  1.85688722]), array([ 1.99433081,  0.81612825]), array([ 1.69033985,  0.64101047]), array([ 1.03335098,  0.16361   ]), array([ 1.92201494,  1.18299267]), array([ 0.26790975,  0.71625702]), array([ 0.56472414,  0.2837483 ]), array([ 1.42118901,  1.16702965]), array([ 0.4000738 ,  0.57757868])]

for i in range(meshsize):
    for j in range(meshsize):
        wfgrid[i,j] = wavefunction([[X[i,j],Y[i,j]]] + otherelectrons)

plot2Dwavefunction(X,Y,wfgrid)
title("wavefunction cross section for %d electrons" %(numelectrons))
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
print mean_wfgrid
for i in range(meshsize):
    for j in range(meshsize):
        #if abs(wavefunction([[i,j]] + otherelectrons)) < mean_wfgrid/100:
        if abs(wfgrid[i,j]) < mean_wfgrid/20:
            zerolocations.append([X[i,j],Y[i,j]])


######
## MAIN ENTREE
######
figure()
zerox,zeroy = array(zerolocations).T
plot(zerox,zeroy,'bo')
otherelectronsplot = array(otherelectrons).T
plot(otherelectronsplot[0],otherelectronsplot[1],'ro')
xlim(0,2)
ylim(0,2)
title("fermion nodes for %d electrons -- total time %fs" %(numelectrons,time.time()-starttime))
savefig("plots/%delectrons-%dmeshsize-%dlength-nodes.png"%(numelectrons,meshsize,length))
savetxt("plots/wfgrid-%delectrons-%dmeshsize-%dlength.txt"%(numelectrons,meshsize,length), wfgrid)

show()





####bisection method
def bisect(f, a, b, e):
	""" Determines zero between a and b using Bisection. """
	n = 0
	fa = f(a)
	if fa == 0.0: return (a, n)
	fb = f(b)
	if fb == 0.0: return (b, n)
		
	while (abs(a-b) > e):
		c = 0.5*(a+b)
		fc = f(c)
		
		if fc == 0.0: return (c, n)
		n = n + 1
		if fb*fc < 0.0:
			a = c
			fa = fc
		
		else:
			b = c
			fb = fc
	
		
	if fa < fb:
		return (a, n)
	else:
		return (b, n)
