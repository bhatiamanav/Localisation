import numpy as np 
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt 

x = np.array([3,4,5,6,7,8,9,10,11,12,13,14,15])
y1=np.load('./kedm_results/S.npy')
y2=np.load('./kedm_results/S1.npy')
y3=np.load('./kedm_results/S2.npy')
y4=np.load('./kedm_results/S3.npy')

#y1=np.load('./kedm_results/S4.npy')
#y2=np.load('./kedm_results/S5.npy')
#y3=np.load('./kedm_results/S6.npy')
#y4=np.load('./kedm_results/S7.npy')


#print(len(y1))
plt.plot(x,y1,'r',label="DGP")
plt.plot(x,y2,'b',label="P=1")
plt.plot(x,y3,'g',label="P=2")
plt.plot(x,y4,'y',label="P=3")
plt.legend(loc="upper left")
plt.xlim([3,15])
plt.ylim([0,1])

plt.savefig('./kedm_results/sparsity_plot.png')

plt.show()
