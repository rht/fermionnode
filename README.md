##Early results
![contour](fermionnode/blob/master/5electrons-100meshsize-2length.png)

![nodes](fermionnode/blob/master/nodes-5electrons-100meshsize-2length.png)

* meshsize = 500, numelectrons = 10 is too slow
* as expected, all the 'position' of the other electrons lies along the node

##Todo
* implement diffusion monte carlo
* google software packages, e.g. piqmc
* summarize the review paper
* really explore http://code.google.com/p/qutip/

##Done
* Apr 6 Plot the slater-det gs wavefunction when N-1 electrons are fixed
* Apr 6 Identify the nodes


##Papers and references
* review paper Fermion nodes http://www.springerlink.com/content/m344615283325722/fulltext.pdf
* fermion monte carlo without fixed nodes
http://jcp.aip.org/resource/1/jcpsa6/v131/i5/p054106_s1
* Structure of fermion nodes http://arxiv.org/pdf/cond-mat/0601485v1.pdf
* Numerical sign problem http://en.wikipedia.org/wiki/Numerical_sign_problem
  Basically, integration of highly oscillatory functions give you headache (this sentence doesn't say anything).


##Read later
* http://prl.aps.org/abstract/PRL/v72/i15/p2442_1
* http://altair.physics.ncsu.edu/articles.htm

##Google keywords 
https://www.google.com/search?ix=aca&sourceid=chrome&ie=UTF-8&q=fermion+node
https://www.google.com/search?aq=f&ix=aca&sourceid=chrome&ie=UTF-8&q=%22fixed-node+quantum+Monte+Carlo%22
https://www.google.com/search?ix=aca&sourceid=chrome&ie=UTF-8&q=free+fermion+node
