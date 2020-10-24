import kedm_funcs
import Helper_Funcs
import numpy as np
trj_output = []
kedm_output = []

class parameters:
    def __init__(self):
        self.d = 2
        self.P = 3
        self.omega = 2*np.pi
        self.mode = 1
        self.T_sampling = np.array([-1, 1])
        self.T_estim = np.array([-1, 1])
        self.N_samples = Helper_Funcs.K(self)
        self.K = Helper_Funcs.K(self)
        self.N_estim= 500
        self.Nr = 5
        self.n_del = 0
        self.sampling = 1
        self.success_prob = 0.99
        self.std = 1
        self.maxIter = 5
        self.n_del_init = 0
        self.bipartite = False
        self.N0 = 3
        self.Pr = 0.9
        self.path = './kedm_results/'
param = parameters()

N=6
A = Helper_Funcs.randomAs(param,N)
colors = np.random.rand(3,N)
fig_name = 'sketch.pdf'
eDi, eDo, eX = kedm_funcs.SketchX(param, A, colors,fig_name,N)
np.save(param.path+'eDi',eDi)
np.save(param.path+'eDo',eDo)
np.save(param.path+'eX',eX)

