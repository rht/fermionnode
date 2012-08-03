import matplotlib
matplotlib.use('TkAgg')
#matplotlib.use('module://mplh5canvas.backend_h5canvas')
from libferminode import *
import matplotlib.animation as animation


################
# Initialization
################

meshsize = 100
length = 2  # size of the box
higheststate = 2
dimension = 2
X = linspace(0, length, num=meshsize)  # initialize coordinates
XX, YY = meshgrid(X, X)
otherelectrons = Electrons(higheststate, length, meshsize=meshsize)

"""
    If you are too lazy to look at libferminode.py,
    the useful methods inside the Electrons
    class are:
    1. updatepos(index, deltaxy)
    2. changepos(index, newxy) -> for moving the electrons around, that's it
    3. eff_wavefunction() -> the reduced wavefunction of the N electron system,
       plot this function if you want to see nodal surface of the electrons
    4. alternatively, you can simply plot with surfaceplot(X, Y,
    eff_wavefunctions)
"""


#########################################
# Plotting the wavefunction cross section
#########################################
report = 2


# report 0
# reproduce 2 electrons in 1D
#############################
if report is 0:
    # what does the wavefunction look like for 2 electrons in 1 dimension?
    # This appears in ceperley's paper somehow
    ion()
    oneD2electrons = lambda x, y: (phi1D(1, length)(x) * phi1D(0, length)(y)-
    phi1D(0, length)(x) * phi1D(1, length)(y))
    contour(XX, YY, oneD2electrons(XX, YY))
    xlabel("x1")
    ylabel("x2")
    savefig("plots/report0/ceperleyfig2.png")
    exit()


# report 1, in 1 line
# electrons swapping
#####################
if report is 1:
    #ion()
    otherelectrons.surfaceplot(XX, YY, otherelectrons.eff_wavefunction(),
            save=0)
    # now swap the electrons
    otherelectrons.pos[0], otherelectrons.pos[1] = otherelectrons.pos[1], otherelectrons.pos[0]
    #raw_input()
    otherelectrons.surfaceplot(XX, YY, otherelectrons.eff_wavefunction(),
            save=0)
    draw()
    #raw_input()
    # moral of the story: exchanging electrons doesn't change the feature of the
    # zeros
    show()
    exit()


# report 2
# interactive plotting
######################

# learned from here and many other google pages, http://kitchingroup.cheme.cmu.edu/software/python/matplotlib/interactive-annotations-in-matplotlib
if report is 2:
    fig = figure()
    directory = tempdir()
    i = 0
    def onclick(event):
        if event.inaxes is None: return
        if event.button != 1: return
        print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %(event.button,
        event.x, event.y, event.xdata, event.ydata)
        otherelectrons.changepos(otherelectrons.active_electron, array([event.xdata, event.ydata]))
        cla()
        otherelectrons.surfaceplot(XX, YY, otherelectrons.eff_wavefunction())
        draw()
        global i
        filename = directory + '/%03d.png' % i
        i+=1
        #savefig(filename)
        #figures.append(filename)

    def onkey(event):
        if event.key in ('q','Q'): exit()
        if event.key in 's': return
        if event.key in 'e': return
            #num = input("give me the number of electrons")
        elif int(event.key) < otherelectrons.N_electrons - 1:
            print "active electron %s" % event.key
            otherelectrons.plot_active_electron('ro')
            otherelectrons.active_electron = int(event.key)
            otherelectrons.plot_active_electron()
            draw()
    otherelectrons.surfaceplot(XX, YY, otherelectrons.eff_wavefunction())
    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('key_press_event', onkey)
    fig.canvas.mpl_connect('motion_notify_event', onclick)
    show()
    #raw_input()
    #createvideofromdirectory(directory)
    exit()


# report 3
# automatically moving the electrons around
###########################################
if report is 3:
    figures = []
    ioff()
    for i in range(10):
        figures.append(figure())
        otherelectrons.updatepos(0, array([.01 * i, 0]))
        otherelectrons.surfaceplot(XX, YY, otherelectrons.eff_wavefunction(),
                save=False)

    createvideo(figures,prefix=True)
    exit()


# report 4
# plot the wavefunction along a parametrized line
#################################################
if report is 4:
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
    show()
