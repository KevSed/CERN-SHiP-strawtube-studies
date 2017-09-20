# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
#from decorators import *
#import sys
#import re

def gauss(x,mean,sigma):
    return np.exp(-1/2*((x-mean)/sigma)**2)

#string1 = 'old'
#string2 = 'new'

string1 = 'nofield'
string2 = '10mT'

#params = [0.4005,43.03,0.91141,44.4343,1.8375,44.247,-5.1355,49.077,7.3157,43.974,0.5124,43.4044,6.044,43.419,-0.3628,42.9651]
params = [0.91141,44.4343,-1.4496,44.1311,-5.1355,49.077,0.2807,45.723,0.5124,43.4044,-1.16431,45.032,-0.3628,42.9651,2.2152,41.763]

x = np.linspace(-90, 90, 1000)

amu = '$\mu^{+}$'
mu = '$\mu^{-}$'

plt.plot(x, gauss(x,params[0],params[1]), 'r-', label=r'y-projection {} {}'.format(amu,string1), color='b')
plt.plot(x, gauss(x,params[2],params[3]), 'r-', label=r'y-projection {} {}'.format(amu,string2))
plt.xlim(x.min(), x.max())
#plt.ylim(f(x, r_1, r_2).min()-0.4, f(x, r_1, r_2).max()
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("hists/nofield/{}/gauss_comparison_y+_{}.pdf".format('comp',string2))
plt.close()

plt.plot(x, gauss(x,params[4],params[5]), 'r-', label=r'y-projection {} {}'.format(mu,string1), color='b')
plt.plot(x, gauss(x,params[6],params[7]), 'r-', label=r'y-projection {} {}'.format(mu,string2))
plt.xlim(x.min(), x.max())
#plt.ylim(f(x, r_1, r_2).min()-0.4, f(x, r_1, r_2).max()
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("hists/nofield/{}/gauss_comparison_y-_{}.pdf".format('comp',string2))
plt.close()

plt.plot(x, gauss(x,params[8],params[9]), 'r-', label=r'x-projection {} {}'.format(amu,string1), color='b')
plt.plot(x, gauss(x,params[10],params[11]), 'r-', label=r'x-projection {} {}'.format(amu,string2))
plt.xlim(x.min(), x.max())
#plt.ylim(f(x, r_1, r_2).min()-0.4, f(x, r_1, r_2).max()
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("hists/nofield/{}/gauss_comparison_x+_{}.pdf".format('comp',string2))
plt.close()

plt.plot(x, gauss(x,params[12],params[13]), 'r-', label=r'x-projection {} {}'.format(mu,string1), color='b')
plt.plot(x, gauss(x,params[14],params[15]), 'r-', label=r'x-projection {} {}'.format(mu,string2))
plt.xlim(x.min(), x.max())
#plt.ylim(f(x, r_1, r_2).min()-0.4, f(x, r_1, r_2).max()
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig("hists/nofield/{}/gauss_comparison_x-_{}.pdf".format('comp',string2))
plt.close()
