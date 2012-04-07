##Early results
![contour](https://github.com/rht/fermionnode/raw/master/plots/5electrons-100meshsize-2length.png)

![nodes](https://github.com/rht/fermionnode/raw/master/plots/nodes-5electrons-100meshsize-2length.png)

* meshsize = 500, numelectrons = 10 is too slow
* as expected, all the 'position' of the other electrons lies along the node
* 24 electrons require 5mins of simulation

##Todo
* implement diffusion monte carlo
* google software packages, e.g. piqmc
* summarize the review paper
* really explore http://code.google.com/p/qutip/
* use weave.blitz() to improve the slater determinant calculation
* urgent: how the zeros of slater determinant is efficiently calculated

##Done

**Apr 6**

*  Plot the slater-det gs wavefunction when N-1 electrons are fixed
* Identify the nodes

**Apr 7**

* switch to use symmetric PBC instead of box boundary
* time the python code using time.time()


##Papers and references
* MAIN review paper Fermion nodes http://www.springerlink.com/content/m344615283325722/fulltext.pdf
* fermion monte carlo without fixed nodes http://jcp.aip.org/resource/1/jcpsa6/v131/i5/p054106_s1
* Structure of fermion nodes http://arxiv.org/pdf/cond-mat/0601485v1.pdf
* H. J. M. van Bemmel et al, "Fixed-node quantum Monte Carlo method for lattice fermions", ![Phys. Rev. Lett. 72, 2442-2445 (1994)](http://prl.aps.org/abstract/PRL/v72/i15/p2442_1)

* Numerical sign problem http://en.wikipedia.org/wiki/Numerical_sign_problem
  Basically, integration of highly oscillatory functions give you headache (this sentence doesn't say anything).
  current strategies: meron-cluster algorithm, stochastic quantization(doesn't work for fermions), fixed node method
* more articles on fixed node method http://www.lorentz.leidenuniv.nl/~saarloos/Correlateds/fixednode.html


##Read later
* http://prl.aps.org/abstract/PRL/v72/i15/p2442_1
* http://altair.physics.ncsu.edu/articles.htm

##Google keywords 
https://www.google.com/search?ix=aca&sourceid=chrome&ie=UTF-8&q=fermion+node
https://www.google.com/search?aq=f&ix=aca&sourceid=chrome&ie=UTF-8&q=%22fixed-node+quantum+Monte+Carlo%22
https://www.google.com/search?ix=aca&sourceid=chrome&ie=UTF-8&q=free+fermion+node
