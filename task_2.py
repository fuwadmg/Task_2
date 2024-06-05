from urllib.request import urlopen
import scipy.special as sps
import matplotlib.pyplot as pp
import numpy as np
import re

class RCS:
    
 def __init__(self, D, f):
     r = D/2
     l = 299792458/f
     k = 2*np.pi/l
     res = 0
     for n in range(1,50):
         res+= (-1)**n*(n+0.5)*(self.b(n, k*r)-self.a(n,k*r))
     self.res = np.abs(res)**2*l**2/np.pi
    
 def h(self, n, x):
     return sps.spherical_jn(n, x)+1j*sps.spherical_yn(n, x)
   
 def a(self, n, x):
     return sps.spherical_jn(n, x)/self.h(n, x)

 def b(self, n, x):
     return (x*sps.spherical_jn(n-1, x)-n*sps.spherical_jn(n, x))/(x*self.h(n-1,x)-n*self.h(n,x))

class Result:
    
    def __init__(self, D, fmin, fmax):
        self.fmin = fmin
        self.fmax = fmax
        self.f = np.linspace(fmin,fmax,1001)
        self.res = RCS(D, self.f).res
        fig, ax = pp.subplots()
        ax.set_xlim(fmin,fmax)
        ax.set_xlabel("f, Гц")
        ax.set_ylabel("σ, м^2")
        ax.plot(self.f, self.res)
        pp.show()
        
    def print(self):
            d = open('res2.xml', 'x')
            text = ['freq', 'lamda', 'rcs']
            val = [self.f, 299792458/self.f, self.res]
            d.write('<?xml version="1.1" encoding="UTF-8" ?>\n<data>\n')
            for a in range(len(self.f)):
                d.write('\t<row>\n')
                for b in range(3):
                    d.write('\t\t<{0}>{1}</{0}>\n'.format(text[b], val[b][a]))
                d.write('\t</row>\n')
            d.write('</data>')
                        
s = str(urlopen('https://jenyay.net/uploads/Student/Modelling/task_rcs.csv').read())
n = '15'
reg = '[0-9]+.?[0-9]*e?-?[0-9]?'
r = re.search('{0}, {1}, {1}, {1}'.format(n, reg), s)
_, D, fmin, fmax = [float(s) for s in r.group().split(', ')]
res = Result(D, fmin, fmax)
res.print()
        
    
     
