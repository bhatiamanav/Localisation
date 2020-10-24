import kedm_funcs
import Helper_Funcs
import numpy as np
import matplotlib.pyplot as plt 
trj_output = []
kedm_output = []

class parameters:
    def __init__(self):
        self.d = 2
        self.P = 1
        self.omega = 2*np.pi
        self.mode = 2
        self.T_sampling = np.array([-1, 1])
        self.T_estim = np.array([-1, 1])
        self.N_samples = Helper_Funcs.K(self)
        self.K = Helper_Funcs.K(self)
        self.N_estim= 500
        self.Nr = 5
        self.n_del = 0
        self.sampling = 1
        self.success_prob = 0.99
        self.std = 0
        self.maxIter = 5
        self.n_del_init = 0
        self.bipartite = False
        self.N0 = 3
        self.Pr = 0.9
        self.path = './kedm_results/'
param = parameters()

Sp_lvl = []
for N in range(4,16):
    print(N)
    param, S = kedm_funcs.Max_Sparsity(param,N)
    Sp_lvl.append(S)
    print('The maximum sparisty level is ', S)
#np.save(param.path+'S3',Sp_lvl)#Uncomment this line and change file name inside '' to store values, then use sparsity_plots
print(Sp_lvl)



