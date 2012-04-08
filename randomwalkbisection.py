from pylab import *
from random import choice

length = 2.
def fn(x): return sqrt(sum((x-length/2)**2))-.5
meshsize = 50
X = linspace(0,length,num=meshsize+1)
grid = arange(meshsize)

class bisect:
    def __init__(self,fn=None):
        self.fn=lambda x: fn(array([X[x[0]],X[x[1]]]))
        self.cursor = array([choice(grid),choice(grid)])
        self.points = [array(self.cursor)]
        self.nodepoints = []
        self.nodewf = []
        self.atnode = False

    def _inbound(self):
        #return 0 <= self.cursor.all() < meshsize
        return (0 <= self.cursor[0] < meshsize) and (0 <= self.cursor[0] < meshsize)
    def _checksignflip(self):
        return self.fn(x1) * self.fn(x2) > 0

    def walk(self):
        cursorbefore = array(self.cursor)
        self.cursor += array([choice([-1,1]), choice([-1,1])])
        if not self._inbound():
            self.cursor = cursorbefore
        self.points.append(array(self.cursor))

    def search1stsignflip(self):
        fninit = self.fn(self.cursor)
        fntemp = float(fninit)
        while fntemp*fninit > 0:
            fntemp = self.fn(self.cursor)
            xtemp = array(self.cursor)
            print X[xtemp[0]],X[xtemp[1]],fntemp
            self.walk()
        self.nodepoints.append(xtemp)
        self.nodepoints.append(self.cursor)
        self.nodewf.append(fntemp)
        self.nodewf.append(self.fn(self.cursor))
        print "found"


    def searchnextsignflip(self):
        fninit = self.nodewf[-1]
        fntemp = float(fninit)
        xtemp = array(self.cursor)
        while fntemp*fninit > 0:
            self.cursor = array(xtemp)
            fntemp = self.fn(self.cursor)
            self.walk()
        self.nodepoints.append(self.cursor)
        self.nodewf.append(fntemp)
        self.atnode = True
        print "found next"

    def searchsignflips(self):
        if not self.atnode:
            self.search1stsignflip()
            self.atnode = True
        for i in range(10):
            self.searchnextsignflip()
            

a = bisect(fn)
a.searchsignflips()
#a.search1stsignflip()



plot(*array([X[i] for i in a.points]).T)
#plot(a.cursor[0],a.cursor[1],'ro',markerfacecolor='green')
gca().add_patch(Circle((length/2.,length/2.),radius=.5,fill=False))
xlim(0,length)
ylim(0,length)
show()
