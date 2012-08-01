#from __future__ import division
import sys
from pyqm import *
from scipy.misc import factorial
#from scipy.optimize import fsolve, broyden1, broyden2, newton_krylov, newton
from scipy.constants import hbar
from itertools import permutations
import time

# QHE
# http://www.itp.phys.ethz.ch/education/lectures_fs11/solid/QHE.pdf
# http://online.physics.uiuc.edu/courses/phys598PTD/fall09/L16.pdf



def phi1D(n,float length):
    '''eigenstate of electron in a 1D periodic box'''
    def result(x):
        #PBC
        return (sqrt(2. / length) *
                    cos(n*2*float64(pi)*x/length))
    return result


def E(n):
    #todo: make this 2D
    return n * n / (8 * length ** 2)


def antisymmetrize(Psis):
    def phi2D(int nx, int ny,float length, float x, float y):
        '''eigenstate of electron in a 2D periodic box.
           n is a tuple for quantum number (nx, ny)
        '''
        return (sqrt(2. / length) *
                        cos(nx*2*pi*x/length) *
                        cos(n[1]*2*pi*y/length))


    '''from the arrays of wavefunction functions, return their slater
    determinant and slater matrix'''
    def psi(xs):
        slatermatrix = matrix(array([[Psi(*x) for Psi in Psis]
                                     for x in xs]))
        slaterlogdet = slogdet(slatermatrix)
        return slaterlogdet[0] * exp(slaterlogdet[1]) / sqrt(factorial(len(xs))), slatermatrix
    return psi


#from scipy.weave import blitz
def wfgridcreator(wavefunction, X, Y):
    '''return grid of a "complicated" function'''
    meshsize = len(X)
    # I have tried:
    # 1. blitz: nope, my code is too abstract
    # 2. don't even wfgrid[0:-1,0:-1] = wavefunction(X[0:-1], Y[0:-1])
    # 3. one line list comprehension
        #wfgrid =  array([[wavefunction(X[i,j],Y[i,j]) for i in indices]
        #for j in indices]
        #4.83s
    # 4. one line
    #wfgrid = wavefunction(X, Y)
         #COOL, this solves the broadcasting headache, but 0.12s slower
         # http://stackoverflow.com/questions/11017347/list-comprehension-for-two-variable-loop-in-python-and-numpy
       
     #old method, too paranoid because of broadcasting
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


def compute_line_param(r1, r2):
    # find the parameter for y = mx + c from two points
    deltaxy = r2 - r1
    m = deltaxy[1] / deltaxy[0]
    c = r2[1] - m * r1[0]
    return m, c


def pathplotter(X, l, wavefunction, otherelectrons):
    # general solution can be found here
    # http://stackoverflow.com/questions/7878398/how-to-extract-an-arbitrary-line-of-values-from-a-numpy-array

    e1 = otherelectrons[0]
    e2 = otherelectrons[-1]
    m, c = compute_line_param(e1, e2)
    Y = m * X + c

    #l = (X[1] - X[0]) * sqrt(1 + m * m) * arange(len(X))
    #projected location of the electrons along the path

    def project_point2line(r, m, c):
        # project point r into line y = mx + c
        y1 = m * r[0] + c
        unitvec = array([1, m]) # here is the fancy way to do it in 1 line
        unitvec /= norm(unitvec) # unitvec = (lambda x:x/norm(x))(array([1, m]))

        deltaxy = (r[1] - y1) * unitvec
        #return array([r[0], y1]) + deltaxy
        return r[0] + deltaxy[0]
    index1 = square(project_point2line(e1, m, c) - X).argmin()
    index2 = square(project_point2line(e2, m, c) - X).argmin()
    wf = array([wavefunction(X[i], Y[i]) for i in range(len(X))])
    plot(l, wf)
    plot(l[index1], [0], 'ro')
    plot(l[index2], [0], 'bo')
    #axhline(0)
    grid()
    xlabel("path")
    ylabel("$\psi$")
    title("wavefunction line section")
    return wf

def zeroplotter(zerox, zeroy, otherelectrons):
    otherelectrons_plot = otherelectrons.T
    plot(zerox, zeroy, 'bo')
    plot(otherelectrons_plot[0],otherelectrons_plot[1],'ro')
    xlim(0,length)
    ylim(0,length)
    #title("fermion nodes for %d electrons -- total time %fs" %(N_electrons,time.time()-starttime))
    #savefig("plots/%delectrons-%dmeshsize-%dlength-nodes-%s.png"%(N_electrons,meshsize,length,timestamp))
    #savetxt("plots/wfgrid-%delectrons-%dmeshsize-%dlength-%s.txt"%(N_electrons,meshsize,length,timestamp), wfgrid)


class Electrons:
    def __init__(self, N_electrons, wf_1particle, length):
        self.starttime = time.time()
        self.N_electrons = N_electrons
        # generate position for background electrons, of dimension (n, 2)
        #self.pos = length * random((N_electrons - 1, 2))
        self.pos = array([[ 0.95045864,  1.39100804],
                [ 0.12410101,  1.73416674],
                [ 1.90899799,  0.76060809],
                [ 0.08314922,  0.32317986],
                [ 0.78967509,  0.12436469]])

        self.wf_1particle = wf_1particle

        # precompute other electrons' determinant once for all
        # what does this physically mean?
        self.precompute_det()


    def precompute_det(self):
        self.det_signed = []
        self.slatermatrix_signed = []
        for i in range(self.N_electrons):
            eachsign = pow(-1, 2 + i)
            # effective wf if you fix other electrons
            det, mat = antisymmetrize(delete(self.wf_1particle, i, 0))(self.pos)
            self.det_signed.append(eachsign * det)
            self.slatermatrix_signed.append(eachsign * mat)
        #return otherelectrons_det_signed, otherelectrons_slatermatrix_signed


    def updatepos(self, index, deltaxy):
        self.pos[index] += deltaxy
        #brute force precompute det after each update
        self.precompute_det()
        #on progress, method to update the det efficiently
        # update slatermatrix http://arxiv.org/pdf/0906.4354.pdf
        #self.slatermatrix_signed_inv = [m.I for m in otherelectrons_slatermatrix_signed]
        #datas = [(otherelectrons, otherelectrons_det_signed, otherelectrons_slatermatrix_signed_inv)]
        #deltar = zeros(shape=shape(otherelectrons))
        #deltar[1] = array([0,.1])
        #R = 1 + array([deltar.T * m for m in otherelectrons_slatermatrix_signed_inv])
        #newmatrixinv = (1 - array([m * deltar for m in otherelectrons_slatermatrix_signed_inv]) / R) * otherelectrons_slatermatrix_signed_inv
        #print(newmatrixinv)

    def changepos(self, index, pos):
        self.pos[index] = pos
        self.precompute_det()


    def eff_wavefunction(self):
        return (lambda x, y: sum([self.wf_1particle[i](x,y) *
            self.det_signed[i] for i in range(self.N_electrons)]))

    def surfaceplot(self, X, Y, wavefunction, save=False):
        '''Compute the 2D surface contour plot of a wavefunction'''
        wfgrid = wfgridcreator(wavefunction, X, Y)
        otherelectrons_plot = self.pos.T
        #pcolor(X,Y,wfgrid)
        wfcontour = contour(X, Y, wfgrid, levels=[0])
        #wfcontour = contour(X, Y, wfgrid)
        plot(otherelectrons_plot[0], otherelectrons_plot[1], 'ro')

        xlabel("x coordinate")
        ylabel("y coordinate")
        newtime = time.time()
        title("wavefunction cross section for %d electrons time %.02fs" %(self.N_electrons,
            newtime-self.starttime))
        self.starttime = newtime
        #colorbar()
        timestamp = time.strftime("%b%d%Y-%H-%M")
        if save:
            savefig("plots/report1/%delectrons-%dmeshsize-%dlength-%s.png"%(self.N_electrons,
                      meshsize, length,timestamp))
        return wfcontour


