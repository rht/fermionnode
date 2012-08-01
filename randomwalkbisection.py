from pylab import *
from random import choice
from mpl_toolkits.mplot3d import Axes3D
from time import time

t0 = time()
length = 2.
def fn(*x): return sqrt((x[0]-length/2)**2 + (x[1]-length/2)**2) - .5
meshsize = 200
X = linspace(0,length,num=meshsize+1)
grid = arange(meshsize)

class bisect:
    def __init__(self,fn=None):
        self.fn = lambda x: fn(X[x[0]],X[x[1]])
        self.cursor = array([choice(grid),choice(grid)])
        self.points = [array(self.cursor)]
        self.nodepoints = []
        self.atnode = False
        self.xdirection = 0
        self.ydirection = 0
        self.diagonal = 0

    def _inbound(self):
        return (0 <= self.cursor[0] < meshsize) and (0 <= self.cursor[1] < meshsize)
    def _checksignflip(self,x1,x2):
        return self.fn(x1) * self.fn(x2) < 0

    def walk(self):
        cursorbefore = array(self.cursor)
        self.cursor += array([choice([-1,1]), choice([-1,1])])
        if not self._inbound():
            self.cursor = cursorbefore
            self.walk()
        else:
            self.points.append(array(self.cursor))

    def search1stsignflip(self):
        cursorbefore = array(self.cursor)
        while not self._checksignflip(cursorbefore,self.cursor):
            cursorbefore = array(self.cursor)
            self.walk()
        self.nodepoints.append(cursorbefore)
        self.nodepoints.append(self.cursor)
        print("found")
        self.xdirection = sign(self.cursor[0] - cursorbefore[0])
        self.ydirection = sign(self.cursor[1] - cursorbefore[1])
        self.diagonal = 1 # 1 for diagonal shift
        self.atnode = True


    def searchnextsignflip(self):
        cursorbefore = self.nodepoints[-2]
        cursornow = self.nodepoints[-1]
        print(cursornow)
        #xmovement
        if self.diagonal:
            #diagonal
            if sign(cursornow[0]-cursorbefore[0]) == self.xdirection:
                cursorused = cursorbefore
                othercursor = cursornow
            else:
                cursorused = cursornow
                othercursor = cursorbefore
            self.cursor = cursorused + [2*self.xdirection,0]
            if self._checksignflip(othercursor,self.cursor):
                self.ydirection *= -1
            else:
                self.diagonal = 0
                self.xdirection *= -1
        else:
            #horizontal
            self.cursor += [self.xdirection,-self.ydirection]
            self.diagonal = 1
            if self._checksignflip(cursornow,self.cursor):
                self.xdirection *= -1
        self.nodepoints.append(self.cursor)

        print("found next")

    def searchsignflips(self):
        if not self.atnode:
            self.search1stsignflip()
        for i in range(20):
            self.searchnextsignflip()


a = bisect(fn)
a.searchsignflips()
#a.search1stsignflip()
A,B=meshgrid(X,X)

figure()
zc = contour(A,B,fn(A,B),levels=[0])
#colorbar()
print "foo", a.nodepoints[-1][0]
print "bar", X[a.nodepoints[-1]]
plot(*array([X[i] for i in a.nodepoints]).T)
plot(*array([X[i] for i in a.points]).T)
#plot(a.cursor[0],a.cursor[1],'ro',markerfacecolor='green')
xmin, ymin = min(array(a.nodepoints).T[0]), min(array(a.nodepoints).T[1])
xmax, ymax = max(array(a.nodepoints).T[0]), max(array(a.nodepoints).T[1])
print xmin, ymin, xmax, ymax
#xlim(0,length)
#xlim(xmin-1,xmax+1)
#ylim(0,length)
#ylim(ymin-1,ymax+1)
print(time()-t0)
show()
