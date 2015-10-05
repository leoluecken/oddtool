
import matplotlib.pyplot as pl
import numpy as np

#name= 'LK2delays';
name='TwoMackeyGlass';

tx = np.genfromtxt("../%s_output/x.txt"%(name), unpack=True)
        
T = tx[0,:]
X = tx[1:,:].transpose()

pl.plot(T,X);
pl.show()


  