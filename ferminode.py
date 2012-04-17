from libferminode import *
import matplotlib.animation as animation


'''use this http://code.enthought.com/chaco/'''
#random cool paper http://www.sciencemag.org/content/318/5852/949.full
#http://physicsandphysicists.blogspot.com/2007/11/simplest-2-slit-and-decoherence.html


def randpos(dimension):  # random position chooser
    return length * random(dimension)


def gridplotter(X, Y, wfgrid, otherelectrons_plot, save=True):
    #plot2Dwavefunction(X,Y,wfgrid)
    #pcolor(X,Y,wfgrid)
    wfcontour = contour(X, Y, wfgrid)
    #zerocontour = wfcontour.collections[3]
    plot(otherelectrons_plot[0], otherelectrons_plot[1],'ro')
    title("wavefunction cross section for %d electrons time %fs" %(numelectrons,time.time()-starttime))
    colorbar()
    timestamp = time.strftime("%d%b%Y-%H-%M")
    if save:
        savefig("plots/%delectrons-%dmeshsize-%dlength-%s.png"%(numelectrons, meshsize, length,timestamp))
    return wfcontour


def zeroplotter(zerox, zeroy, otherelectrons_plot):
    plot(zerox, zeroy, 'bo')
    plot(otherelectrons_plot[0],otherelectrons_plot[1],'ro')
    xlim(0,length)
    ylim(0,length)
    #title("fermion nodes for %d electrons -- total time %fs" %(numelectrons,time.time()-starttime))
    #savefig("plots/%delectrons-%dmeshsize-%dlength-nodes-%s.png"%(numelectrons,meshsize,length,timestamp))
    #savetxt("plots/wfgrid-%delectrons-%dmeshsize-%dlength-%s.txt"%(numelectrons,meshsize,length,timestamp), wfgrid)

starttime = time.time()
meshsize = 100
length = 2  # size of the box
higheststate = 3
numelectrons = int(factorial(higheststate))  # + 1 # number of permutations
dimension = 2

#initialize coordinates
X = linspace(0, length, num=meshsize)
Y = X
#K = linspace(-pi,pi,num=meshsize)
otherelectrons = array([randpos(dimension) for i in range(numelectrons - 1)])
#print(otherelectrons)
#otherelectrons = [array([.2,.3])]
otherelectrons_plot = otherelectrons


###############
# Plotting the wavefunction cross section
###############
X, Y = meshgrid(X, X)

wf_1particle = array([phi2D(i, length) for i in permutations(range(1,higheststate+1))])
#wf_1particle = append(wf_1particle, phi2D((5,1),length))

# precompute other electrons' determinant once for all
# what does this physically mean?
otherelectrons_det_signed = []
otherelectrons_slatermatrix_signed = []
for i in range(numelectrons):
    eachsign = pow(-1, 2 + i)
    det, mat = antisymmetrize(delete(wf_1particle, i, 0))(otherelectrons)
    otherelectrons_det_signed.append(eachsign * det)
    otherelectrons_slatermatrix_signed.append(eachsign * mat)
wavefunction = lambda x, y: sum([wf_1particle[i](x,y) * otherelectrons_det_signed[i] for i in range(numelectrons)])

print(wavefunction(.2,.3))


#############
# MAIN ENTREE
#############

#wfgrid = wfgridcreator(wavefunction, X, Y)
#gridplotter(X, Y, wfgrid, otherelectrons_plot, save=False)
#show()

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
# update slatermatrix http://arxiv.org/pdf/0906.4354.pdf
otherelectrons_slatermatrix_signed_inv = [m.I for m in otherelectrons_slatermatrix_signed]
datas = [(otherelectrons, otherelectrons_det_signed, otherelectrons_slatermatrix_signed_inv)]
deltar = zeros(shape=shape(otherelectrons))
deltar[1] = array([0,.1])
R = 1 + array([deltar.T * m for m in otherelectrons_slatermatrix_signed_inv])
#newmatrixinv = (1 - array([m * deltar for m in otherelectrons_slatermatrix_signed_inv]) / R) * otherelectrons_slatermatrix_signed_inv
#print(newmatrixinv)

