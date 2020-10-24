#Importing the libraries
import cvxpy as cvx
import os, sys
import numpy as np
import Helper_Funcs
import math

#Defining the parameters of KEDM and the Trajectory
class KEDM_OUT:
    def __init__(self):
        self.eo         = None
        self.ei         = None
        self.G          = None
        self.sampling_times   = None
        self.estim_times   = None
        self.tau_list   = None
        self.A          = None
        self.status     = None
class TRJ_OUT:
    def __init__(self):
        self.A_         = None
        self.eX         = None

#Defining the KEDM function 
def KEDM(param, A,N):
    d = param.d
    K = Helper_Funcs.K(param)
    
    MaxIter = 500
    verbosity = True
    
    output = KEDM_OUT()
    
    #Generating the required functions of time
    tau_list = Helper_Funcs.generalsampler(param, 'required')
    sampling_times = Helper_Funcs.generalsampler(param, 'sample_interv')
    estim_times = Helper_Funcs.generalsampler(param, 'estim_interv')
    log_list = Helper_Funcs.generalsampler(param, 'basic_start')
    
    #Generating the required data
    D, W, ei = Helper_Funcs.generateData(param, sampling_times, A,N)
    
    G = []
    con = []
    for k in range(K):
        G.append(cvx.Variable((N, N), PSD=True))#Generating a random variable and transforming it into a column vector
        con.append(cvx.sum(G[k],axis = 0) == 0)
    for t in log_list:
        weights = Helper_Funcs.W(t, tau_list, param)#Defining the weights in the given times
        weights = weights/np.linalg.norm(weights)#Normalize weights
        G_tot = Helper_Funcs.G_t(G, weights, True)#Generate Gt
        con.append(G_tot>>0)#Ensure positive semidefiniteness
    cost = 0
    for i_t, t in enumerate(sampling_times):#During the sampling_times times, we sample the matrix at given times
        weights = Helper_Funcs.W(t, tau_list, param)
        G_tot = Helper_Funcs.G_t(G, weights, True)
        con.append(G_tot >> 0)
        D_G = Helper_Funcs.gram2edm(G_tot, N, True)
        W_vec = np.diag(W[i_t].flatten())
        alpha = (np.linalg.norm( np.matmul(W_vec, D[i_t].flatten()) ) )**(-2)#Significance of different weights
        cost += alpha*cvx.norm( cvx.matmul(W_vec, cvx.vec(D[i_t]-D_G) ) )**2#Finding the total error, taking into consideration significance of different weights 
    
    #Defining the CVX problem
    obj = cvx.Minimize(cost)
    prob = cvx.Problem(obj,con)

    #Solve the problem, if unable to solve show error message
    try:
        prob.solve(solver=cvx.CVXOPT,verbose=False,normalize = True)
    except Exception as message:
        print(message)

    output.status = str(prob.status)
    #If we can find an optimal solution, save the output parameters
    if str(prob.status) == 'optimal':
        G = Helper_Funcs.rankProj(G, param,N)
        eo = Helper_Funcs.testError(param, tau_list, estim_times, G, A,N)
        output.G        = G
        output.tau_list = tau_list
        output.sampling_times = sampling_times
        output.estim_times = estim_times
        output.ei       = ei
        output.eo       = eo
        output.A        = A
    return output

#Finding Max sparsity given the number of points
def Max_Sparsity(param,N):
    n_del_init = param.n_del_init
    maxIter = param.maxIter
    Pr = param.Pr
    success_prob = param.success_prob
    
    #If number of edges to be deleted is greater than the number of edges then stop
    if n_del_init >= Helper_Funcs.edgeCnt(param,N):
        print('Fix n_del_init!')
        return param
    cnt_threshold = maxIter - maxIter*Pr#The threshold at which localozation is incorrect
    cnt_threshold = math.ceil(cnt_threshold)
    param.n_del = n_del_init
    while param.n_del < Helper_Funcs.edgeCnt(param,N):
        cnt = 0
        cnt_wrong = 0
        while cnt < maxIter:
            A = Helper_Funcs.randomAs(param,N)
            output = KEDM(param, A,N)
            cvx_status = output.status
            if cvx_status == 'optimal':
                error_out = output.eo
                cnt += 1
                if Helper_Funcs.mean(error_out) > 1-success_prob:#If output error is greater than 1-success_threshold, increase the counts for wrong
                    cnt_wrong = cnt_wrong + 1
                if maxIter - cnt < cnt_threshold - cnt_wrong:#If enough number of tries are not left,skip the chance and reset
                    break
                if cnt_wrong > cnt_threshold:#If localization is incorrect,delete the edge and return the sparsity
                    param.n_del -= 1
                    S = param.n_del/(N*(N-1)/2)
                    return param, S
        print('n_del is ', param.n_del)
        param.n_del += 1#Increase number of deleted pointd by one
    param.n_del -= 1
    S = param.n_del/(N*(N-1)/2)#Return sparsity
    return param, S

#Define trajectory,use algorithm 2 
def Retrieve_Traj(param, kedm_output,N):
    d = param.d
    P = param.P
    N_samples = param.N_samples
    mode = param.mode
    omega = param.omega
    K = param.K
    
    anchor_idx = Helper_Funcs.randomAnchor(N_samples, N, d+1)#Generate the random anchors
    
    tau_list = kedm_output.tau_list
    sampling_times = kedm_output.sampling_times
    estim_times = kedm_output.estim_times
    
    A = kedm_output.A
    G = kedm_output.G
    N_estim = param.N_estim
    output = TRJ_OUT()
    
    M = np.shape(anchor_idx)[0]#If number of anchors is not equal to times of sampling, it is wrong
    if M != len(sampling_times):
        print('Fix this!')
        return output
    
    X_ = []
    T = []
    cnt_wrong = 0
    cnt = 0
    for m in range(M):
        list_m = np.nonzero(anchor_idx[m,:])[0]#Find number of non-zero anchors
        if len(list_m) < d:#If number of non-zero anchors are less than embedding dimension, then this is wrong
            cnt_wrong += 1
            continue
        t = sampling_times[m]#Take a sampling time
        Xt = Helper_Funcs.X_t(A,t,param,N)#Generate Xt
        Xtm = Xt[:,list_m]#define known anchors
        #Define weights,Gt and Xt real
        weights = Helper_Funcs.W(t, tau_list, param)
        Gt = Helper_Funcs.G_t(G, weights, False)
        Xt_ = Helper_Funcs.gram2x(Gt, d)
        Xt_m = Xt_[:,list_m]
        R = Helper_Funcs.rotationXY(Xtm, Xt_m)#Find required rotation matrix
        Nt = Xtm.shape[1]
        Xt_ = np.matmul(R,Xt_ - np.matmul(Xt_m,np.ones((Nt,N))/Nt) ) + np.matmul(Xtm,np.ones((Nt,N))/Nt)#Define the trajectory
        Tt = Helper_Funcs.T_t(t, param)#Define the time
        if cnt == 0:
            X_ = Xt_
            T = Tt
        else:
            X_ = np.concatenate((X_, Xt_), axis=0)
            T = np.concatenate((T, Tt), axis=0)
        cnt += 1
    if cnt_wrong == M:#If all trials were wrong, this is wrong
        print('Wrong!')
        return output
    T_inv = np.linalg.pinv(T)
    A_ = np.matmul(T_inv,X_)#A=T^-1*X
    A_ = Helper_Funcs.deconcat(A_, param)#Take trajectory dimension to P


    if mode == 2:#If polynomial model, take fourier transform
        A_ = Helper_Funcs.sine2fourier(A_)
    
    #Find trajectory mismatch
    eX = np.zeros(N_estim)
    for i_t, t in enumerate(estim_times):
        X_ = Helper_Funcs.X_t(A_,t,param,N)
        X = Helper_Funcs.X_t(A,t,param,N)
        eX[i_t] = np.linalg.norm(X-X_,'fro')/np.linalg.norm(X,'fro')

    output.A_ = A_
    output.eX = eX

    return output

#Define Sketch function
def SketchX(param, A,colors, figName,N):
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d
    
    #if d is not 2 or 3, it is rong
    d = param.d
    if d > 3 or d < 2:
        print('Error')
        return
    M = 300#Number of realizations
    maxIter = param.maxIter
    T_estim = param.T_estim
    t = np.linspace(T_estim[0],T_estim[1],M)
    
    
    fSize = 18
    
    #If 2D plot, plot in 2D, else plot in 3D
    if d == 2:
        plt.xlabel('x_axis',fontsize=18)
        plt.ylabel('y_axis',fontsize=18)
    else:
        ax = plt.axes(projection='3d')
        ax.set_xlabel('x_axis',fontsize=fSize)
        ax.set_ylabel('y_axis',fontsize=fSize)
        ax.set_zlabel('z_axis',fontsize=fSize)
    
    i = 0
    eX = 0
    eDo = 0
    eDi = 0
    while i < maxIter:
        kedm_output = KEDM(param, A,N)
        cvx_status = kedm_output.status
        if cvx_status == 'optimal':
            i += 1
            print(i,'out of', maxIter)
            trj_output = Retrieve_Traj(param,kedm_output,N)#Take trajectory input 
            eX = eX + Helper_Funcs.mean(trj_output.eX)/maxIter#Trajectory mismatch
            eDo = eDo + Helper_Funcs.mean(kedm_output.eo)/maxIter#Output error
            eDi = eDi + Helper_Funcs.mean(kedm_output.ei)/maxIter#Input error
            A_ = trj_output.A_
            X = np.zeros((M,d,N))
            for m in range(M):
                X[m,:,:] = Helper_Funcs.X_t(A_,t[m],param,N)#Plot the trajectory of x
            for n in range(N):
                if d == 2:
                    plt.plot(X[:,0,n], X[:,1,n],c = colors[:,n],linewidth=1/(maxIter**0.5))#If d=2,plot in 2D
                else:
                    ax.plot3D(X[:,0,n], X[:,1,n],X[:,2,n],c = colors[:,n],linewidth=1/(maxIter**0.5))#iF d=3,plot in 3D
                plt.pause(0.01)

    plt.savefig(param.path+figName)
    plt.close()
    return eDi, eDo, eX
