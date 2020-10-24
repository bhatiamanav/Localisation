import kedm_funcs
import Helper_Funcs
import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
import cvxpy as cvx
trj_output = []
kedm_output = []

class parameters:
    def __init__(self):
        self.d = 2
        self.P = 3
        self.omega = 2*np.pi
        self.mode = 2
        self.T_sampling = np.array([0, 1])
        self.T_estim = np.array([0, 1])
        self.N_samples = Helper_Funcs.K(self)
        self.K = Helper_Funcs.K(self)
        self.N_estim= 500
        self.Nr = 5
        self.n_del = 0
        self.sampling = 2
        self.success_prob = 0.99
        self.std = 0
        self.maxIter = 5
        self.n_del_init = 0
        self.bipartite = False
        self.N0 = 3
        self.Pr = 0.9
        self.path = './kedm_results/'
param = parameters()

N=10
M=200
eo_tot = np.zeros(param.N_estim)
temp=M

for i in range(M):
    print(i)
    A = Helper_Funcs.randomAs(param,N)
    kedm_output = kedm_funcs.KEDM(param,A,N)
    if(kedm_output.status=='optimal'):
        eo_tot = np.array(eo_tot)
        eo = np.array(kedm_output.eo)
        eo_tot = np.add(eo_tot,eo)
    else:
        temp = temp - 1
        
    print(eo_tot)

eo_tot = eo_tot / temp
eo_tot = eo_tot / N
np.save(param.path+'error_rand_cheby',eo_tot)

x_new = np.linspace(-1, 1, param.N_estim)

plt.plot(x_new,eo_tot)
plt.xlim([-1,1])
plt.ylim([0,0.07])

plt.show()



