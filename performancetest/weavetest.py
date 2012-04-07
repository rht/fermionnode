from timeit import Timer

exp = "a = 0\nfor i in range(10000):a+=i*i"
exp2 = "from scipy import weave,arange\na=arange(10000)\nweave.blitz('''sum(a*a)''')"

t1 = Timer(exp)
print t1.timeit(500)


#with blitz
t2 = Timer(exp2)
print t2.timeit(500)
