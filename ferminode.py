from libferminode import *
import matplotlib.animation as animation


'''use this http://code.enthought.com/chaco/'''
#random cool paper http://www.sciencemag.org/content/318/5852/949.full
#http://physicsandphysicists.blogspot.com/2007/11/simplest-2-slit-and-decoherence.html


def surfaceplotter(X, Y, wavefunction, otherelectrons, save=True):
    wfgrid = wfgridcreator(wavefunction, X, Y)
    otherelectrons_plot = otherelectrons.T
    #plot2Dwavefunction(X,Y,wfgrid)
    #pcolor(X,Y,wfgrid)
    wfcontour = contour(X, Y, wfgrid)
    #zerocontour = wfcontour.collections[3]
    plot(otherelectrons_plot[0], otherelectrons_plot[1],'ro')
    xlabel("x coordinate")
    ylabel("y coordinate")
    title("wavefunction cross section for %d electrons time %fs" %(numelectrons,
        time.time()-starttime))
    colorbar()
    timestamp = time.strftime("%b%d%Y-%H-%M")
    if save:
        savefig("plots/report1/%delectrons-%dmeshsize-%dlength-%s.png"%(numelectrons,
                  meshsize, length,timestamp))
    return wfcontour


def compute_line_param(r1, r2):
    # find the parameter of y = mx + c from two points
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
    def project(r, m, c):
        y1 = m * r[0] + c
        unitvec = array([1, m])
        unitvec /= norm(unitvec)
        deltaxy = (r[1] - y1) * unitvec
        #return array([r[0], y1]) + deltaxy
        return r[0] + deltaxy[0]
    index1 = square(project(e1, m, c) - X).argmin()
    index2 = square(project(e2, m, c) - X).argmin()
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
    #title("fermion nodes for %d electrons -- total time %fs" %(numelectrons,time.time()-starttime))
    #savefig("plots/%delectrons-%dmeshsize-%dlength-nodes-%s.png"%(numelectrons,meshsize,length,timestamp))
    #savetxt("plots/wfgrid-%delectrons-%dmeshsize-%dlength-%s.txt"%(numelectrons,meshsize,length,timestamp), wfgrid)


################
# Initialization
################

starttime = time.time()
meshsize = 100
length = 2  # size of the box
higheststate = 3
numelectrons = int(factorial(higheststate))  # number of permutations
dimension = 2

# initialize coordinates
X = linspace(0, length, num=meshsize)
Y = X

class Electrons:
    def __init__(self, numelectrons, wf_1particle):
        self.numelectrons = numelectrons
        # generate position for background electrons, of dimension (n, 2)
        self.pos = length * random((numelectrons - 1, 2))
        print(self.pos)
        self.wf_1particle = wf_1particle

        # precompute other electrons' determinant once for all
        # what does this physically mean?
        self.precompute_det()

    def precompute_det(self):
        self.det_signed = []
        self.slatermatrix_signed = []
        for i in range(self.numelectrons):
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


    def eff_wavefunction(self):
        return lambda x, y: sum([self.wf_1particle[i](x,y) *
            self.det_signed[i] for i in range(self.numelectrons)])




##########################
# wavefunction computation
##########################
# array of wavefunctions used in total ground state wf
wf_1particle = array([phi2D(i, length) for i in permutations(range(1,higheststate+1))])
#wf_1particle = append(wf_1particle, phi2D((5,1),length))

otherelectrons = Electrons(numelectrons, wf_1particle)



##################
# PLOTTING TIME!!!
# Plotting the wavefunction cross section
##################

# report 1
##########
#XX, YY = meshgrid(X, X)
#surfaceplotter(XX, YY, otherelectrons.eff_wavefunction(), otherelectrons.pos, save=False)


# report 2
##########
m, c = compute_line_param(otherelectrons.pos[0], otherelectrons.pos[1])
l = (X[1] - X[0]) * sqrt(1 + m * m) * arange(len(X))

wf1 = pathplotter(X, l, otherelectrons.eff_wavefunction(), otherelectrons.pos)
otherelectrons.updatepos(0, array([.1, .1 * m]))

# plot again the difference
wf2 = pathplotter(X, l, otherelectrons.eff_wavefunction(), otherelectrons.pos)

# again, plot again the difference
otherelectrons.updatepos(0, array([.1, .1 * m]))
wf3 = pathplotter(X, l, otherelectrons.eff_wavefunction(), otherelectrons.pos)

figure()
plot(l, wf2 - wf1)
show()



#datas = [(X, Y, wfgrid, otherelectrons_plot)]
#for i in range(10):
    #datas.append((X, Y, wfgrid, array(otherelectrons_plot)))

#datas = []
#for i in range(20):
    #zerox,zeroy = bruteforcezerofinder(wfgrid,length,meshsize,X,Y)
    #otherelectrons_plot += array([(0,.05)]).T
    #wfgrid = wfgridcreator(wavefunction, X, Y, list(otherelectrons_plot.T), meshsize)
    #datas.append((zerox, zeroy, array(otherelectrons_plot)))
#createvideo(datas, zeroplotter, directory="animatenode/")



##############################
# moving otherelectrons around
##############################

