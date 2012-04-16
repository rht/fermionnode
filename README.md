##Early results
![nodes](https://github.com/rht/fermionnode/raw/master/plots/5electrons-100meshsize-2length-nodes.png)

* as expected, all the 'position' of the other electrons lies along the node
* as of April 7, 24 electrons require 5mins of simulation
* comment: "the nodal surface for a real wavefunction is always closed loops. (Notice that the apparent ends at the edges of the box always have partners on the other side which is related by the periodic boundary conditions.) This is so because the nodes separate the regions of positive from regions of negative wavefunction and therefore can't end."
* as of April 8, 24 electrons require 2mins of simulation
  full performance report on April 8, 24e, 143s:
      * psi 5.882s ncalls 10000
      * result 114.826 ncalls 5760000
  6e require 13.37s
* performance report on April 16:
      * 6e require 2.73

##Todo
* implement diffusion monte carlo
* summarize the review paper
* really explore http://code.google.com/p/qutip/
* urgent: how the zeros of slater determinant is efficiently calculated. **Consult http://arxiv.org/abs/0804.2161 for optimization**. Alternatively, check out http://code.google.com/p/pi-qmc/. But then, scipy's det also uses LU decomposition
* wtf is 'bosonization' in ceperley's path integral? Pauli surface?
* How different is the usual antisymmetrization vs superposition? (and also entanglement). check out Neumaier's ans in physics.SE
* ground state and "overdamping" is a bad analogy, but I need to know how meaningless it really is.
* read http://prl.aps.org/abstract/PRL/v72/i15/p2442_1
* read http://altair.physics.ncsu.edu/articles.htm
* read Path integrals in the theory of condensed helium http://rmp.aps.org/abstract/RMP/v67/i2/p279_1
* read Quantum Monte Carlo simulations of solics http://rmp.aps.org/abstract/RMP/v73/i1/p33_1
* read A fast and efficient algorithm for Slater determinant updates in quantum Monte Carlo simulations http://jcp.aip.org/resource/1/jcpsa6/v130/i20/p204105_s1
* See Ceperley publications http://people.physics.illinois.edu/Ceperley/papers/index.htm
* read http://prl.aps.org/pdf/PRL/v72/i15/p2442_1
* read http://prb.aps.org/pdf/PRB/v44/i17/p9410_1
* see http://code.google.com/p/pi-qmc/source/browse/trunk/src/fixednode/FreeParticleNodes.cc -> not useful
* read http://physics.stackexchange.com/questions/11605/when-is-the-minus-sign-problem-in-quantum-simulations-an-obstacle



##Done
**Apr 6**

*  Plot the slater-det gs wavefunction when N-1 electrons are fixed
* Identify the nodes

**Apr 7**

* switch to use symmetric PBC instead of box boundary
* time the python code using time.time()

**Apr 16**

* implement reduced slater determinant

##Abandoned
* use scipy.weave.blitz() to improve the slater determinant calculation. Meh, doesn't help that much
* http://stackoverflow.com/questions/9337188/animating-a-contour-plot-in-matplotlib-using-funcanimation  not flexible. That's it
```python
fig = figure()
ims = []
for i in range(2):
    otherelectronsplot += array([(0,.05)]).T
    zerox,zeroy = bruteforcezerofinder(wfgrid,length,meshsize,X,Y)
    wfgrid = wfgridcreator(wavefunction, X, Y, list(otherelectronsplot.T), meshsize)
    ims.append(plot(zerox,zeroy,'bo'))
ani = animation.ArtistAnimation(fig, ims, interval=50, blit = True, repeat_delay = 1000)
ani.save('blah.mp4')
```


##Papers and references
* MAIN review paper Fermion nodes http://www.springerlink.com/content/m344615283325722/fulltext.pdf
* fermion monte carlo without fixed nodes http://jcp.aip.org/resource/1/jcpsa6/v131/i5/p054106_s1
* Structure of fermion nodes http://arxiv.org/pdf/cond-mat/0601485v1.pdf
* H. J. M. van Bemmel et al, "Fixed-node quantum Monte Carlo method for lattice fermions", [Phys. Rev. Lett. 72, 2442-2445 (1994)](http://prl.aps.org/abstract/PRL/v72/i15/p2442_1)

* Numerical sign problem http://en.wikipedia.org/wiki/Numerical_sign_problem
  Basically, integration of highly oscillatory functions give you headache (this sentence doesn't say anything).
  current strategies: meron-cluster algorithm, stochastic quantization(doesn't work for fermions), fixed node method
* more articles on fixed node method http://www.lorentz.leidenuniv.nl/~saarloos/Correlateds/fixednode.html
* [Antisymmetrizer](http://en.wikipedia.org/wiki/Antisymmetrizer)

###Electronic Wavefunction Calculations
* Electronic Wave Functions. I. A General Method of Calculation for the Stationary States of Any Molecular System http://rspa.royalsocietypublishing.org/content/200/1063/542.short


##Google keywords 
* fermion node
* fixed node quantum monte carlo
* free fermion node
* slater determinant monte carlo

##Node searching strategy
1. expensive: the most rudimentary one is by using cutoff method abs(wavefunction) < cutoff
2. http://www.mathworks.com/matlabcentral/answers/15082-multiple-solutions-using-fsolve basically, square the wavefunction, and find the minimum using fmin
3. bisection method. http://stackoverflow.com/questions/4326667/solving-equation-using-bisection-method -> but this example is only for 1D. more: http://stackoverflow.com/questions/3513660/multivariate-bisection-method
4. random walk bisection method. will be implemented soon
