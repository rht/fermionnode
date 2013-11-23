#references
#Time evolution
##Numerical approaches to time evolution of complex quantum systems
##introduces Crank-Nicolshon scheme
##http://www.sciencedirect.com/science/article/pii/S0375960109004927

from pylab import *
#from fipy import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.misc import factorial
import scipy.sparse as sparse
import os
import time
import subprocess
import sh
from matplotlib import animation

#utilities
##1. timing.
class tic:
    #This can be used either in a context, or as a function
    #http://stackoverflow.com/questions/5849800/tic-toc-functions-analog-in-python
    def __init__(self):
        self.tic = time.time()
    def __enter__(self):
        self.tic = time.time()
    def __exit__(self):
        print('Elapsed: %s s' % (time.time() - self.tic))
    def __iter__(self):
        return self
    #def next(self):
        #self.toc = time.time()
        #print("time spent", self.toc - self.tic)
        #self.tic = self.toc
    def __call__(self):
        self.toc = time.time()
        print("Elapsed: %f s"% (self.toc - self.tic))
        self.tic = self.toc


##2. get source code
def getsourcecode(fn):
    import inspect
    print ''.join(inspect.getsourcelines(fn)[0])
##3. send plot to imgur
def imgur():
    savefig("/tmp/pic.png")
    print sh.tee(sh.imgur2("/tmp/pic.png"),"pbcopy")
    print sh.echo("copied to clipboard")






#physics
## QM
def groundstate(H):
    """
    Classical mechanics has been developed continuously from the time of Newton
    and applied to an ever-widening range of dynamical systems, including the
    electromagnetic field in interaction with matter. The underlying ideas and
    the laws governing their application form a simple and elegant scheme, which
    one would be inclined to think could not be seriously modified without
    having all its attractive features spoilt. Nevertheless it has been found
    possible to set up a new scheme, called quantum mechanics, which is more
    suitable for the description of phenomena on the atomic scale and which is
    in some respects more elegant and satisfying than the classical scheme. This
    possibility is due to the changes which the new scheme involves being of a
    very profound character and not clashing with the features of classical
    theory that make it so attractive, as a result of which all these features
    can be incorporated in the new scheme.
    """
    """find the gs of a Hamiltonian matrix"""
    Hinv = H.I
    E, wf = sparse.linalg.eigsh(H, k = 1)
    return 1/E, wf
commutator = lambda x,y: x * y - y * x

##Pauli spin matrices SU(2)
sx = mat("0 1; 1 0")
sy = mat("0 -1j; 1j 0")
sz = mat("1 0; 0 -1")
Pauli = [sx, sy, sz]
smin = mat("0 1; 0 0")
smax = mat("0 0; 1 0")


##Angular momentum matrices SO(3)
#http://quantummechanics.ucsd.edu/ph130a/130_notes/node247.html
sqrt2 = sqrt(2)
Lx = mat("0 1 0; 1 0 1; 0 1 0") / sqrt2
Ly = mat("0 1 0; -1 0 1; 0 -1 0") / sqrt2 / 1j
Lz = mat("1 0 0; 0 0 0; 0 0 -1") / sqrt2


#Gell-Mann matrices SU(3)
Gellmann = [mat("0 1 0; 1 0 0; 0 0 0"),
           mat("0 -1j 0; -1j 0 0; 0 0 0"),
           mat("1 0 0; 0 -1 0; 0 0 0"),
           mat("0 0 1; 0 0 0; 1 0 0"),
           mat("0 0 -1j; 0 0 0; 1j 0 0"),
           mat("0 0 0; 0 0 1; 0 1 0"),
           mat("0 0 0; 0 0 -1j; 0 1j 0"),
           mat("1 0 0; 0 1 0; 0 0 -2") / sqrt(3)]


##Gamma matrices
# weyl basis
weyl = matrix([[0, 1],[ -1, 0]])
pauli = lambda x: [sx, sy, sz][x - 1]
index = [0, 1, 2, 3]
gamma0 = kron(sx, eye(2))
gamma5 = -kron(sz, eye(2))
gamma = lambda x: gamma0 if x == 0 else kron(weyl, pauli(x))


#QFT
def fourdot(p, k): #four vector dot product
    return p[0] * k[0] - dot(p[1:], k[1:])
def fourvector(threevector, m): #generate on-shell fourvector only
    return array([sqrt(m ** 2 + norm(threevector) ** 2)] + list(threevector))
def fermion_u(p, spin, m):
    # Peskin 46
    #http://en.wikipedia.org/wiki/Dirac_spinor
    E, p1, p2, p3 = p
    if spin == 1:
        return sqrt(E + m) * array([1, 0, p3 / (E + m),
            (p1 + 1j * p2) / (E + m)])
    else:
        return sqrt(E + m) * array([0, 1, p1 - 1j * p2 / (E + m),
            -p3 / (E + m)])
def fermion_ubar(p, spin, m):
    return fermion_u(p, spin, m).T
def antifermion_v(p, spin, m):
    # Peskin 46
    #http://en.wikipedia.org/wiki/Dirac_spinor
    E, p1, p2, p3 = p
    "?? spin 1?"
    if spin == 1:
        return sqrt(E + m) * array([p1 - 1j * p2 / (E + m), -p3 / (E + m), 0, 1])
    else:
        return sqrt(E + m) * array([p3 / (E + m), (p1 + 1j * p2) / (E + m), 1,
            0])
def slash(k):
    return matrix(add.reduce([gamma(i) * k[i] for i in
        range(4)]))
def antifermion_vbar(p, spin, m):
    return antifermion_v(p, spin, m).T





#memes
def eigsort(H):
    #sort eigen vectors and values in the order of increasing eigenvalues
    vals, vecs = eig(H)
    return sort(vals), vecs.T[vals.argsort()]


def eighsort(H):
    vals, vecs = eigh(H)
    return sort(vals), vecs.T[vals.argsort()]


def speighsort(H):
    # sparse version of eighsort
    vals, vecs = sparse.linalg.eigsh(H, k=6, which='SM')
    return sort(vals), vecs.T[vals.argsort()]


def quickplot_surface(X, Y, Z):
    #plot_surface without worrying of specifying a 3d figure
    ax = figure().add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z)
    return ax

def renormalize(x): return x/norm(x)


def mesh1D(dx,nx):
    # create 1D mesh
    L = nx * dx
    #return linspace(-L,L,nx)
    return arange(-L/2.,L/2.,dx)


def mesh2D(dx,nx):
    L = nx * dx
    X = arange(-(L-1)/2, (L+1)/2,dx)
    return meshgrid(X,X) 


def Kmesh2D(dx,nx):
    Kx = linspace(-pi, pi, num = nx)
    return Kx, Kx


def gaussian(x, sigma):
    return  exp(-x*x/(2.*sigma**2))/sqrt(2.*pi*sigma**2)


def xtop(wavefn):
    '''convert position space wavefn into momentum space wavefn'''
    #http://docs.scipy.org/doc/numpy/reference/routines.fft.html#background-information
    return fftshift(abs(fftn(wavefn)))


def ptox(wavefn):
    return fftshift(abs(ifftn(wavefn)))


def derivative(f,dx):
    return convolve([1,-1],f,mode='same') / dx


def centralderivative(dx, nx, order=2):
    """
    Cubum autem in duos cubos, aut quadratoquadratum in duos quadratoquadratos,
    et generaliter nullam in infinitum ultra quadratum potestatem in duos
    eiusdem nominis fas est dividere cuius rei demonstrationem mirabilem sane
    detexi. Hanc marginis exiguitas non caperet.
    """
    def diaones(offset): return diag(ones(nx-abs(offset)), offset)
    if order is 2:
        mat = .5 * matrix(diaones(1) - diaones(-1)) 
        mat[0,-1] = -.5
        mat[-1,0] = .5
        mat /= dx
    else: #4th order
        mat = 1./12 * matrix(-diaones(2) + 8*diaones(1) - 8*diaones(-1) +
                diaones(-2))
        mat /= dx
    return mat


def spectralderivative(dx, nx):
    from scipy.linalg import toeplitz
    column = concatenate(([0], .5 * (-1) ** arange(1,nx) / tan(arange(1,nx) * dx / 2)))
    return toeplitz(column, -column)


def K1D(dx,nx, order=2, BC='pbc',enablesparse=0):
    """1D laplacematrix"""
    #ref http://www.mech.kth.se/~ardeshir/courses/literature/fd.pdf
    #dirichlet boundary condition
    #http://stackoverflow.com/questions/4843034/application-of-boundary-conditions-in-finite-difference-solution-for-the-heat-eq
    #2nd order equation 19
    #http://www.mathworks.com/matlabcentral/fileexchange/27279-laplacian-in-1d-2d-or-3d/content/laplacian.m
    # new since nov 21 2012, sparse
    # http://docs.scipy.org/doc/scipy/reference/sparse.html
    #http://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors_of_the_second_derivative
    if enablesparse:
        if BC=='pbc':
            K = sparse.diags([ones(nx-1), ones(nx-1), -2*ones(nx), 1, 1],
                             [1,-1,0,nx-1,-nx+1])
        elif BC=='dd':
            K = sparse.diags([ones(nx-1), ones(nx-1), -2*ones(nx)],
                             [1,-1,0])
        elif BC=='dn':
            K = sparse.diags([ones(nx-1), ones(nx-1),
                concatenate(-2*ones(nx-1), [-1])], [1,-1,0])
        elif BC=='nd':
            K = sparse.diags([ones(nx-1), ones(nx-1),
                concatenate([-1], -2*ones(nx-1))], [1,-1,0])
        elif BC=='nn':
            K = sparse.diags([ones(nx-1), ones(nx-1),
                concatenate([-1], -2*ones(nx-2), [-1])], [1,-1,0])
        return K / (dx*dx)

    else:
        def diaones(offset): return diag(ones(nx-abs(offset)), offset)
        if order == 2:
            K = matrix(diaones(1) + diaones(-1) + diaones(0)*-2)
            if BC=='pbc':
                K[0,-1] = K[-1,0] = 1
            elif (BC=='dd') or (BC=='fixed'):  # dirichlet-dirichlet
                # dirichlet-dirichlet also means fixed
                pass
            elif BC=='dn':
                K[-1,-1] = -1
            elif BC=='nd':
                K[0,0] = -1
            elif BC=='nn':
                K[0,0] = K[-1,-1] = -1
            return K / (dx*dx)
        #4th order
        else:
            return matrix( -diaones(2) + 16 * diaones(1) - 30*diaones(0) +
                    16*diaones(-1) - diaones(-2)
                    ) / (12*dx*dx)

#http://astrophysics.fic.uni.lodz.pl/100yrs/pdf/12/069.pdf
def laplace(x,dx):
    # python implementation of laplace matrix
    ##from scipy.ndimage.filters import laplace
    #the line above doesn't work because our wavefunction is complex
    return convolve(x,[1,-2,1],mode='same') / (dx*dx)


# 2 dimensional laplacematrix
#http://math.mit.edu/classes/18.086/2006/project1-Dominguez.pdf
#laplacian matrix
#http://stackoverflow.com/questions/1764859/how-to-compute-laplacian-of-a-field
#http://en.wikipedia.org/wiki/Laplacian_matrix
def K2D(dx, nx, BC='pbc'):
    #which means, the operator needs to be flattened
    K = K1D(dx, nx, BC=BC, enablesparse=1)
    I = sparse.identity(nx)
    return sparse.kron(K, I) + sparse.kron(I, K)


def direct_sum(A,B):
    res = zeros(add(A.shape, B.shape),dtype=complex)
    res[:A.shape[0],:A.shape[1]] = A
    res[A.shape[0]:,A.shape[1]:] = B
    return res


def kronecker(*args):
    _sparse_ = True
    kronoperator = sparse.kron if _sparse_ else kron

    if len(args) is 2:
        return kronoperator(*args)
    else:
        return kronoperator(args[0], kronecker(*args[1:]))


def K3D(dx, nx, BC='pbc'):
    K = K1D(dx, nx, BC=BC, enablesparse=1)
    I = sparse.identity(nx)
    return kronecker(K, I, I) + kronecker(I, K, I) + kronecker(I, I, K)


def KnD(dx, nx, n=3, BC='pbc'):
    K = K1D(dx, nx, BC=BC, sparse=1)
    I = sparse.identity(nx)
    #http://stackoverflow.com/questions/1550130/cloning-row-or-column-vectors
    #KI = [I,]*(n-1)
    #KI.append(K)
    KI = repeat(I[newaxis,:,:],(n-1),0)
    KI = concatenate((KI, [K,]), axis=0)
    return add.reduce([kronecker(*roll(KI,i, axis=0)) for i in range(n)])


def plot2Dspectrum(Hamiltonian, K):
    spectrum = array([ sort(real(eigvals(Hamiltonian(i)))) for i in K])
    figure()
    for i in range(shape(spectrum)[1]):
        plot(K, spectrum[:,i])


def animate1D(X,Z):
    fig = figure()
    #line,foo = plot(X,Z[0],X,potential(X))
    line, = plot(X,Z[0])
    for i in Z[1:]:
        waitforbuttonpress()
        line.set_ydata(i)
        draw()


def plotsurfaceanimate(X,Y,Z):
    fig = figure()
    ax = Axes3D(fig)
    ax.set_zlim3d(0,.12)
    for i in Z:
        ax.plot_surface(X,Y,i)
        draw()


def createvideo(xs=None,ys=None, prefix=None, outputdir='.'):
    pref = time.strftime("%H-%M-%b%d%Y")
    writer = animation.writers['ffmpeg'](fps=15)
    fig = figure()
    if xs is None:
        l, = plot(ys[0])
        with writer.saving(fig, outputdir+'/out'+prefix+pref+'.mp4', len(ys)-1):
            for i in range(1, len(ys)):
                l.set_ydata(ys[i])
                writer.grab_frame()
    else:
        l, = plot(xs[0], ys[0])
        with writer.saving(fig, outputdir+'/out'+prefix+pref+'.mp4', len(ys)-1):
            for i in range(1, len(ys)):
                l.set_data(xs[i], ys[i])
                writer.grab_frame()


def tempdir():
    import tempfile
    return tempfile.gettempdir()


def createvideofromdirectory(directory,prefix='',outputdir='.'):
    #http://dawes.wordpress.com/2007/12/04/animating-png-files/
    #http://stackoverflow.com/questions/4092927/generating-movie-from-python-without-saving-individual-frames-to-files
    #http://www.scipy.org/Cookbook/Matplotlib/Animations
    #note that os.spawnvp is deprecated
    command = ('ffmpeg','-i', directory + '/%03d.png', outputdir+
            '/%sout.mp4'%prefix, '-vcodec',
            'mpg4', '-vf', '"setpts=40.0*PTS"', '-y', '-r', '1')
    subprocess.call(command)


##########################
###New since April 2012###
##########################

def plot2Dwavefunction(X,Y, wf):
    ax = quickplot_surface(X, Y, wf)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('$\psi$')


#####
# QHO
#####
#obsolete, use scipy's "physicsist" hermite polynomial instead
#def hermite(n, x): # return H n(x)
    #if n < 0:
        #return 0.0
    #elif n == 0:
        #return 1.0
    #else:
        #return 2*x*hermite(n-1, x) - 2*(n-1)*hermite(n-2, x)
from scipy.special import hermite

def QHO_psi_n(n, x, m, omega):
    """1D QHO eigenfunctions"""
    #http://en.wikipedia.org/wiki/Quantum_harmonic_oscillator
    #nu = m * omega / hbar
    nu = m * omega
    # normalization coefficient
    C =  (nu/pi)**(1./4) * sqrt(1/(2**n*factorial(n)))
    return C * exp(-nu* x**2 /2) * hermite(n)(sqrt(nu)*x)

def QHO_psi2D(n, x, m, omega):
    """2D QHO eigenfunctions"""
    nu = m * omega
    C =  ((nu/pi)**(1./4)) ** 2 * sqrt(1/(2**array(n)*factorial(n)))
    return (C * exp(-nu* norm(x)**2 /2) * hermite(n[0])(sqrt(nu)*x[0]) *
            hermite(n[1])(sqrt(nu)*x[1]))



#############
# Solid State
#############
# todo
# done
# calculate reciprocal lattice vector
# density of state
# implement brillouin zone path

#common data
#lattice vectors
sc_lattices = array(mat("1 0 0; 0 1 0; 0 0 1"))
fcc_lattices = array(mat("0 .5 .5; .5 0 .5; .5 .5 0"))
bcc_lattices = array(mat("-.5 .5 .5; .5 -.5 .5; .5 .5 -.5"))
#basis vectors
fcc_basiss = array([[0,0,0]])
bcc_basiss = array(mat("0 0 0; .5 .5 .5"))



def get_reciprocals(vecs):
    # compute reciprocal lattice vectors
    # http://en.wikipedia.org/wiki/Reciprocal_lattice
    b = []
    nvec = len(vecs)
    for i in range(nvec):
        j = (i+1) % nvec
        k = (i+2) % nvec
        b.append(2*pi*cross(vecs[j], vecs[k]) / dot(vecs[i],
                cross(vecs[j], vecs[k])))
    return b


def get_density_of_state(m_H, E):
    # bonus
    # grosso section 1.4 equation 48b
    # could have been done more efficiently
    def green_function(E):
        return inv(E - m_H)
    return -1 / pi * imag(trace(green_function(E + 0.01j)))


def get_unitcells(lattices, basiss):
    unitcells = list(basiss)
    for v in lattices:
        for b in basiss:
            unitcells.append(b+v)
    return unitcells


def get_FK(h,k,l, lattices, unitcells):
    reciprocals = get_reciprocals(lattices)
    #unitcells = get_unitcells(lattices, basiss)
    K = reciprocals * array([h,k,l])
    f = 1
    S = f* sum([exp(-1j*sum(dot(K, v))) for v in unitcells])
    if imag(S) < 1e-10:
        if real(S) < 1e-10:
            S = 0
        else: S = real(S)
    elif real(S) < 1e-10:
        S = 1j*imag(S)
    return S


def get_kpath(crystal_type, point_names, Npoints=50):
    # see http://en.wikipedia.org/wiki/Brillouin_zone
    # brillouin zone critical points
    # http://en.wikipedia.org/wiki/File:Brillouin_Zone_(1st,_FCC).svg
    # orthorhombic lattice critical points
    # see https://www.msu.edu/~dodat/files/Brillouin_zone.pdf
    # see also http://users.mrl.uiuc.edu/floren/Thesis/Appendices.pdf
    # see also http://www.cryst.ehu.es/
    # see also http://www.aflowlib.org/
    # http://users.mrl.uiuc.edu/floren/Thesis/Appendices.pdf
    critical_points = {
    "orc" : {"Gamma" : [0, 0, 0],
             "X" : [1, 0, 0],
             "Y" : [0, 1, 0],
             "Z" : [0, 0, 1],
             "R" : [1, 1, 1],
             "S" : [1, 1, 0],
             "T" : [0, 1, 1],
             "U" : [1, 0, 1]},
    "sc" : {"Gamma" : [0, 0, 0],
            "X" : [.5, 0, 0],
            "M" : [.5, 0, .5],
            "R" : [.5, .5, .5]},
    "fcc" : {"Gamma" : [0, 0, 0], #http://lamp.tu-graz.ac.at/~hadley/ss1/bzones/fcc.php
             "X" : [1., 0, 0],
             "K" : [.75, 0, .75],
             "L" : [.5, .5, .5],
             "W" : [1, 0, .5],
             "U/K" : [1, 1, 0],
             "U" : [1, -.25, .25]},
    "bcc": {"Gamma" : [0, 0, 0],
            "H" : [0, 0, 1],
            "N" : [0, .5, .5],
            "P" : [.5, .5, .5]}
    }

    def create_path(initial_k, final_k):
        initial_k, final_k = 2 * pi * array(initial_k), 2 * pi * array(final_k)
        Delta_k = final_k - initial_k
        #Npoints = int(norm(Delta_k) / dk)
        delta_k = Delta_k / Npoints
        kpoints = [initial_k + i * delta_k for i in range(Npoints)]
        return kpoints

    all_kpoints = []
    all_critpoints = [0]
    points = [critical_points[crystal_type][i] for i in point_names]
    for i in range(len(points) - 1):
        all_kpoints += create_path(points[i], points[i+1])
        all_critpoints.append(len(all_kpoints))
    def criticalpointsplot():
        for i in range(len(point_names)):
            axvline(x=all_critpoints[i], label=point_names[i])
    return array(all_kpoints)

#Gell-Mann matrices SU(4)
#see http://prb.aps.org/pdf/PRB/v85/i19/e195132
#http://ocw.mit.edu/courses/mathematics/18-712-introduction-to-representation-theory-fall-2010/lecture-notes/
Gellmann4 = [direct_sum(G,array([0])) for G in Gellmann]
Gellmann4.append([mat("0 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 0"),
                  mat("0 0 0 -1j; 0 0 0 0; 0 0 0 0; 1j 0 0 0"),
                  mat("0 0 0 0; 0 0 0 1; 0 0 0 0; 0 1 0 0"),
                  mat("0 0 0 0; 0 0 0 -1j; 0 0 0 0; 0 1j 0 0"),
                  mat("0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 1 0"),
                  mat("0 0 0 0; 0 0 0 0; 0 0 0 -1j; 0 0 1j 0"),
                  mat("1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -3") / sqrt(6)])


