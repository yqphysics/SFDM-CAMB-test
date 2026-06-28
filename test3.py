from matplotlib import pyplot as plt
import numpy as np


x = np.linspace(0.1,10,1000)

sig = 0.1

x0 = 5

P1 = np.exp(- (np.log(x/x0))**2/(2*sig**2))/np.sqrt( 2* np.pi*sig**2)


P2 = np.exp(- (x-x0)**2/(2*sig**2))/np.sqrt( 2* np.pi*sig**2)




plt.plot(x,P1)
plt.plot(x,P2)
plt.show()