from libferminode import *
import matplotlib.animation as animation
ion()


'''use this http://code.enthought.com/chaco/'''
#random cool paper http://www.sciencemag.org/content/318/5852/949.full
#http://physicsandphysicists.blogspot.com/2007/11/simplest-2-slit-and-decoherence.html


################
# Initialization
################

meshsize = 100
length = 2  # size of the box
higheststate = 5
dimension = 2

# initialize coordinates
X = linspace(0, length, num=meshsize)
Y = X

otherelectrons = Electrons(higheststate, length, meshsize=meshsize)
XX, YY = meshgrid(X, X)



#########################################
# Plotting the wavefunction cross section
#########################################
report0, report1, report2 = 0, 1, 0

# report 0
##########
if report0:
    # how does the wavefunction look like for 2 electrons in 1 dimension?
    oneD2electrons = lambda x, y: (phi1D(1, length)(x) * phi1D(0, length)(y)-
    phi1D(0, length)(x) * phi1D(1, length)(y))
    contour(XX, YY, wfgridcreator(oneD2electrons, XX, YY))
    xlabel("x1")
    ylabel("x2")
    savefig("plots/report0/ceperleyfig2.png")
    exit()


# report 1, in 1 line
##########
if report1:
    otherelectrons.surfaceplot(XX, YY, otherelectrons.eff_wavefunction(),
            save=False)
    raw_input()
    exit()

# report 2
##########

#1. Interactive plotting
# learned from here http://kitchingroup.cheme.cmu.edu/software/python/matplotlib/interactive-annotations-in-matplotlib
fig = figure()
#import copy
#figures = []
def onclick(event):
    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %(event.button,
    event.x, event.y, event.xdata, event.ydata)
    otherelectrons.changepos(0, array([event.xdata, event.ydata]))
    cla()
    otherelectrons.surfaceplot(XX, YY, otherelectrons.eff_wavefunction())
    draw()
    #figures.append(copy.deepcopy(fig))

otherelectrons.surfaceplot(XX, YY, otherelectrons.eff_wavefunction())
fig.canvas.mpl_connect('button_press_event', onclick)
raw_input()
#createvideo(figures)
exit()



#2. moving the electrons around

figures = []
for i in range(10):
    figures.append(figure())
    otherelectrons.updatepos(0, array([.1, .00001 * i]))
    otherelectrons.surfaceplot(XX, YY, otherelectrons.eff_wavefunction(),
            save=False)

createvideo(figures)
exit()




m, c = compute_line_param(otherelectrons.pos[0], otherelectrons.pos[1])
l = (X[1] - X[0]) * sqrt(1 + m * m) * arange(len(X))

figure()
wf1 = pathplotter(X, l, otherelectrons.eff_wavefunction(), otherelectrons.pos)
otherelectrons.updatepos(0, array([.1, .1 * m]))

# plot again the difference
wf2 = pathplotter(X, l, otherelectrons.eff_wavefunction(), otherelectrons.pos)

# again, plot again the difference
#otherelectrons.updatepos(0, array([.1, .1 * m]))
#wf3 = pathplotter(X, l, otherelectrons.eff_wavefunction(), otherelectrons.pos)

figure()
plot(l, wf2 - wf1)
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

