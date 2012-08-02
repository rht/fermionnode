import sys
from pylab import *
#from pyqm import createvideo, createvideofromdirectory, tempdir
#from scipy.optimize import fsolve, broyden1, broyden2, newton_krylov, newton
from scipy.constants import hbar
from scipy.misc import factorial
from itertools import permutations, product
import time

# QHE
# http://www.itp.phys.ethz.ch/education/lectures_fs11/solid/QHE.pdf
# http://online.physics.uiuc.edu/courses/phys598PTD/fall09/L16.pdf


def phi2D(n, length):
    '''
       eigenstate of electron in a 2D periodic box.
       n is a tuple for quantum number (nx, ny)
       real part, hence e^ikx and e^-ikx should look the same
    '''
    def result(*x):
        #PBC
        #return (sqrt(2. / length) *
                    #cos(n[0]*2*float64(pi)*x[0]/length) *
                    #cos(n[1]*2*float64(pi)*x[1]/length))
        return (sqrt(2. / length) *
                    exp(1j*2*pi/length* dot(n, x)))
    return result

def phi2Dlessfancy(n, length, x):
    return (sqrt(2. / length) *
            exp(1j*2*pi/length*dot(n, x)))


def phi1D(n, length):
    '''eigenstate of electron in a 1D periodic box'''
    def result(x):
        #PBC
        return (sqrt(2. / length) *
                    exp(1j*n*2*float64(pi)*x/length))
    return result


def E(n):
    #todo: make this 2D
    return n * n / (8 * length ** 2)


def antisymmetrize(Psis):
    '''from the arrays of wavefunction functions, return their slater
    determinant and slater matrix'''
    def psi(xs):
        slatermatrix = matrix([[Psi(*x) for Psi in Psis]
                                     for x in xs])

        slaterlogdet = slogdet(slatermatrix)
        return slaterlogdet[0] * exp(slaterlogdet[1]) / sqrt(factorial(len(xs))), slatermatrix
    return psi


#from scipy.weave import blitz
#old
def wfgridcreator(wavefunction, X, Y):
    '''return grid of a "complicated" function, might be useful in the future'''
    #meshsize = len(X)
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



# more helper functions
def createvideo(figures, prefix=None):
    #http://dawes.wordpress.com/2007/12/04/animating-png-files/
    #http://stackoverflow.com/questions/4092927/generating-movie-from-python-without-saving-individual-frames-to-files
    #http://www.scipy.org/Cookbook/Matplotlib/Animations
    import tempfile
    directory = tempfile.gettempdir()
    pref = 0
    if prefix: pref = time.strftime("%b%d%Y")
    else: pref = ''

    os.spawnvp(os.P_WAIT, 'trash', ('trash', directory + '/*'))
    #http://forum.videohelp.com/threads/306745-Slow-motion-with-ffmpeg
    #http://ffmpeg.org/trac/ffmpeg/wiki/How%20to%20speed%20up%20/%20slow%20down%20a%20video
    command = ('ffmpeg','-i', directory + '/%03d.png', 'out%s.mp4' % pref, '-vcodec',
            'mpg4', '-vf', '"setpts=40.0*PTS"', '-y', '-r', '1')
    #command = ('convert', directory + '/%03d.png', 'out.gif')
    #command = ('mencoder', 'mf:/'+directory+'/%03.png', '-speed', '0.4', '-mf',
            #'w=800:h=600:fps=25:type=png', '-ovc', 'lavc', '-lavcopts',
            #'vcodec=mpeg4:mbd2:trell', '-oac', 'copy', '-o', 'output.avi' )
    # -y is for auto-overwrite
    #convert -delay 50 Th*.JPG anim.mpg

    for i in range(len(figures)):
        filename = directory + '/%s%03d.png'%(pref, i)
        figures[i].savefig(filename)
        #print('Wrote file '+ filename)
        clf()
    os.spawnvp(os.P_WAIT, 'ffmpeg', command)
    #os.spawnvp(os.P_WAIT, 'mencoder', command)
    #os.spawnvp(os.P_WAIT, 'convert', command)

def tempdir():
    import tempfile
    return tempfile.gettempdir()


def createvideofromdirectory(directory):
    command = ('ffmpeg','-i', directory + '/%03d.png', 'out.mp4', '-vcodec',
            'mpg4', '-vf', '"setpts=40.0*PTS"', '-y', '-r', '1')
    os.spawnvp(os.P_WAIT, 'ffmpeg', command)




class Electrons:
    def __init__(self, higheststate, length, meshsize=None):
        self.active_electron = 0
        self.length = length
        self.higheststate = higheststate
        self.meshsize = meshsize
        self.starttime = time.time()
#http://www.theo3.physik.uni-stuttgart.de/lehre/ss08/sst/Script-AHCN-Chap-6.pdf
        self.quantum_numbers = array([i for i in
            product(range(-self.higheststate, self.higheststate+1),
                repeat=2) if norm(i) <= self.higheststate])

        self.N_electrons = len(self.quantum_numbers)
        print "number of electrons %d" %self.N_electrons

        # wavefunction computation
        # array of wavefunctions used in total ground state wf
        self.wf_1particle = array([phi2D(i, length) for i in self.quantum_numbers])



        # generate position for background electrons, of dimension (n, 2)
        self.pos = length * random((self.N_electrons - 1, 2))
        print self.pos
        #self.pos = array([[ 1.17253759,  0.76248715],
                #[ 0.55659908,  0.8677818 ],
                #[ 1.70763184,  0.35243496],
                #[ 1.84601349,  0.84027975],
                #[ 0.37452466,  1.80036508],
                #[ 0.67844901,  0.73389094],
                #[ 0.89556134,  1.83766219],
                #[ 1.84195774,  1.58432459],
                #[ 0.5413078,   0.12694015],
                #[ 1.39434316,  1.21330462],
                #[ 1.49001942,  1.09599095],
                #[ 0.95245113,  0.26030348]])
        

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


    #performancetestjuly31
    #combine precompute_det and antisymmetrize together
    def precompute_det(self):
        self.det_signed = []
        self.slatermatrix_signed = []
        #crafted from the old wf_1particle
        almost_matrix = array([[phi2Dlessfancy(j, self.length, x) for j in
            self.quantum_numbers] for x in self.pos])

        for i in range(self.N_electrons):
            eachsign = pow(-1, 2 + i)
            # effective wf if you fix other electrons

            mat = matrix(delete(almost_matrix, i, 1))

            #Antisymmetrization
            slaterlogdet = slogdet(mat)
            det = (slaterlogdet[0] * exp(slaterlogdet[1]) /
                    sqrt(factorial(len(self.pos))))

            self.det_signed.append(eachsign * det)
            self.slatermatrix_signed.append(eachsign * mat)


    def updatepos(self, index, deltaxy):
        self.pos[index] += deltaxy
        #brute force precompute det after each update
        self.precompute_det()
        #on progress, method to update the det efficiently
        # -> not really required, because this is not the bottleneck of the
        # computation
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


    #performancetweakaug1
    def eff_wavefunction(self):
#        return vectorize(lambda x,y: dot(sqrt(2. / self.length) *
                    #exp(1j*2*pi/self.length*dot(self.quantum_numbers,(x,y))),
            #array(self.det_signed)))
        return vectorize(lambda x,y: dot(phi2D(self.quantum_numbers,
            self.length)(x,y),
            array(self.det_signed)))


    def plot_active_electron(self,color='go'):
        plot(self.pos[self.active_electron][0],
                self.pos[self.active_electron][1], color)


    def surfaceplot(self, X, Y, wavefunction, save=False):
        '''Compute the 2D surface contour plot of a wavefunction'''
        #wfgrid = wfgridcreator(wavefunction, X, Y)
        #performancetweakaug1
        wfgrid = wavefunction(X, Y)
        otherelectrons_plot = delete(self.pos, self.active_electron, 0).T
        #wfcontour = pcolor(X,Y,wfgrid)
        #wfcontour = contourf(X, Y, wfgrid, colors=('r','b'))
        wfcontour = contour(X, Y, wfgrid, levels=[0])
        plot(otherelectrons_plot[0], otherelectrons_plot[1], 'ro')
        self.plot_active_electron()

        xlabel("x coordinate")
        ylabel("y coordinate")
        newtime = time.time()
        title("wavefunction cross section for %d electrons time %.02fs" %(self.N_electrons,
            newtime-self.starttime))
        self.starttime = newtime
        #colorbar()
        timestamp = time.strftime("%b%d%Y-%H-%M-%S")
        if save:
            #savefig("plots/report2/%delectrons-%dmeshsize-%dlength-%s.png"%(self.N_electrons,
            savefig("plots/5electrons/%delectrons-%dmeshsize-%dlength-%s.png"%(self.N_electrons,
                      self.meshsize, self.length,timestamp))
        return wfcontour


