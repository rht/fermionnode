from libferminode import *
import matplotlib.animation as animation


'''use this http://code.enthought.com/chaco/'''
#random cool paper http://www.sciencemag.org/content/318/5852/949.full
#http://physicsandphysicists.blogspot.com/2007/11/simplest-2-slit-and-decoherence.html


def randpos(dimension):  # random position chooser
    return length * random(dimension)


def gridplotter(X, Y, wfgrid, otherelectronsplot, save=True):
    #plot2Dwavefunction(X,Y,wfgrid)
    #pcolor(X,Y,wfgrid)
    wfcontour = contour(X, Y, wfgrid)
    #zerocontour = wfcontour.collections[3]
    plot(otherelectronsplot[0], otherelectronsplot[1],'ro')
    title("wavefunction cross section for %d electrons time %fs" %(numelectrons,time.time()-starttime))
    colorbar()
    timestamp = time.strftime("%d%b%Y-%H-%M")
    if save:
        savefig("plots/%delectrons-%dmeshsize-%dlength-%s.png"%(numelectrons, meshsize, length,timestamp))
    return wfcontour


def zeroplotter(zerox, zeroy, otherelectronsplot):
    plot(zerox, zeroy, 'bo')
    plot(otherelectronsplot[0],otherelectronsplot[1],'ro')
    xlim(0,length)
    ylim(0,length)
    #title("fermion nodes for %d electrons -- total time %fs" %(numelectrons,time.time()-starttime))
    #savefig("plots/%delectrons-%dmeshsize-%dlength-nodes-%s.png"%(numelectrons,meshsize,length,timestamp))
    #savetxt("plots/wfgrid-%delectrons-%dmeshsize-%dlength-%s.txt"%(numelectrons,meshsize,length,timestamp), wfgrid)

starttime = time.time()
meshsize = 100
length = 2  # size of the box
higheststate = 4
numelectrons = int(factorial(higheststate)) + 1 # number of permutations
dimension = 2

#initialize coordinates
X = linspace(0, length, num=meshsize)
Y = X
#K = linspace(-pi,pi,num=meshsize)
otherelectrons = [randpos(dimension) for i in range(numelectrons - 1)]
#print(otherelectrons)
#otherelectrons = [array([.2,.3])]
otherelectronsplot = array(otherelectrons).T


###############
# Plotting the wavefunction cross section
###############
X, Y = meshgrid(X, X)

wf_1particle = array([phi2D(i, length) for i in permutations(range(1,higheststate+1))])
wf_1particle = append(wf_1particle, phi2D((5,1),length))

# precompute other electrons' determinant once for all
# what does this physically mean?
signed_otherelectrons_determinant = array([pow(-1, 2 + i) * antisymmetrize(delete(wf_1particle, i, 0))(otherelectrons) for i in range(numelectrons)])
print(signed_otherelectrons_determinant)
wavefunction = lambda x, y: sum([wf_1particle[i](x,y) * signed_otherelectrons_determinant[i] for i in range(numelectrons)])

print(wavefunction(.2,.3))
#wfgrid = wfgridcreator(wavefunction, X, Y)
#gridplotter(X, Y, wfgrid, otherelectronsplot, save=False)
#show()

#datas = [(X, Y, wfgrid, otherelectronsplot)]
#for i in range(10):
    #datas.append((X, Y, wfgrid, array(otherelectronsplot)))


#############
# MAIN ENTREE
#############


#datas = []
#for i in range(20):
    #zerox,zeroy = bruteforcezerofinder(wfgrid,length,meshsize,X,Y)
    #otherelectronsplot += array([(0,.05)]).T
    #wfgrid = wfgridcreator(wavefunction, X, Y, list(otherelectronsplot.T), meshsize)
    #datas.append((zerox, zeroy, array(otherelectronsplot)))
#createvideo(datas, zeroplotter, directory="animatenode/")

