import numpy as np
import cvxpy as cvx

def generalsampler(param, sampler):#Generating the sampling times for the matrices and the times to calculate the error
    import math
    K = param.K
    mode = param.mode
    omega = param.omega
    sampling = param.sampling
    T0 = 2*np.pi/omega
   
    if sampler == 'required':
        N = K
        T = param.T_estim
    elif sampler == 'sample_interv':
        N = param.N_samples
        T = param.T_sampling
    elif sampler == 'estim_interv':
        N = param.N_estim
        T = param.T_estim
    else:
        N = param.Nr
        T = np.array([0,45])
   
    if mode == 2:
        T[1] = min(T[1]-T[0],T0)+T[0]
   
    if sampler == 'required' or sampler == 'estim_interv':
        if mode == 1:#Equispaced
            interval = np.linspace(T[0],T[1],N)
        else:
            interval = np.linspace(T[0],T[1],N+1)
            interval = interval[0:N]
    elif sampler == 'sample_interv':
        if mode == 1:#Equispaced
            interval = np.linspace(T[0],T[1],N)
        else:
            interval = np.linspace(T[0],T[1],N+1)
            interval = interval[0:N]
        if sampling == 2:
            if mode == 1:#Chebyshev
                for i in range(N):
                    interval[i] = (T[0]+T[1])/2 + (T[1]-T[0])/2 * math.cos((2*i+1)/(2*N) * np.pi)
            else:
                interval = np.linspace(T[0],T[1],N+1)
                for i in range(N+1):
                    interval[i] = (T[0]+T[1])/2 + (T[1]-T[0])/2 * math.cos((2*i+1)/(2*(N+1)) * np.pi)
                interval = interval[0:N]
        elif sampling == 3:#random
            interval = np.random.uniform(T[0],T[1],N)
    else:
        if mode == 1:#Logarithmic time sampling
            Nc = int(np.ceil(N/2))
            Nf = int(np.floor(N/2))
            interval_p = np.exp(np.random.uniform(T[0],T[1],Nc)-5) - np.exp(-5)
            interval_n = -np.exp(np.random.uniform(T[0],T[1],Nf)-5) + np.exp(-5)
            interval = np.concatenate((interval_n, interval_p), axis=0)
        else:
            interval = np.random.uniform(T[0],T[1],N)

    return interval


def edgeCnt(param,N):
    if param.bipartite:#If graph is partite, nodes are cleanly divided into two, and edges exist between one node and all the other nodes of the other group,so edges = (N0)*(N-N0)
        N0 = param.N0
        N1 = N-N0
        M = N0*N1
    else:#For a bipartite graph, there is a edge between every pair of vertices(for a EDM)
        M = N*(N-1)//2
    return M

def gram2edm(G, N, mode):#Generating the EDM from the gramian matrix,using the formula
    if mode:
        dG = cvx.vstack( cvx.diag(G) )
        D = cvx.matmul(dG ,np.ones((1,N)))-2*G + cvx.matmul(np.ones((N,1)),dG.T)
    else:
        dG = np.diag(G)[:,None]
        D = dG-2*G+dG.T
    return D

def symNoise(N):#Simulating the normalized gaussian random noise with variance 1
    Noise = np.random.randn(N,N)
    Noise = ( 2 ** (-1/2) ) * (Noise + Noise.T)
    return Noise

def projD(D):#Rank projection is diag(diag(D))
    D = D - np.diag(np.diag(D))
    D[D <= 0.0] = 0.0#Ensuring positive semidefinite
    return D

def randMask(param,N):#Generating random mask
    W = np.zeros((N,N))
    n_del = param.n_del
    if param.bipartite:#If the graph is bipartite
        N0 = param.N0
        N1 = N - N0
        M = edgeCnt(param,N)
        k = np.random.choice(np.arange(M),size=(n_del), replace = False)#Generate a random choice from a 1D array
        Wb = np.ones((N0,N1))
        ij = np.where(Wb)
        for i in k:
            Wb[ij[0][i], ij[1][i]] = 0#Set random points to 0
        W[0:N0,N-N1:] = Wb
        W = W + W.T#Generate the mask
    else:#For a non-bipartite graph
        M = edgeCnt(param,N)
        k = np.random.choice(np.arange(M),size=(n_del), replace = False)
        ij = np.where(np.tril(np.ones((N,N)), k=-1))#Find the lower triangular matrix equal to 1
        for i in k:
            W[ij[0][i], ij[1][i]] = 1
        W = 1 - (W+W.T)#Set random points to 0
    return W

def generateData(param,sampling_times,A,N):#Use for generating the Data of a rank projection of D,random mask, and the error between the edm and its rank projection
    N_samples = len(sampling_times)
    std = param.std
    
    D = np.zeros((N_samples,N,N))
    W = np.zeros((N_samples,N,N))
    ei = np.zeros(N_samples)
    
    for i_t, t in enumerate(sampling_times):
        Xi = X_t(A,t,param,N)#Generating the Trajectory coeffs.
        Gi = np.matmul(Xi.T,Xi)#G=XT.X
        Di = gram2edm(Gi, N, False)#The edm matrix
        D[i_t] = projD(Di + std*symNoise(N))#Rank projection of a noisy model of D
        W[i_t] = randMask(param,N)#Random mask to be deleted
        ei[i_t] = np.linalg.norm(Di-D[i_t],'fro') / (np.linalg.norm(Di,'fro') )#Error between edm and its noisy rank projection
    return D, W, ei

def W(t, tau_list, param):#Generating the weight matrix
    #Generate using propositions 1 and 2
    P = param.P
    mode = param.mode
    omega = param.omega
    K = param.K

    M = np.zeros((K,K))
    T = np.zeros((K,1))
    if mode == 1:#For polynomial models
        for i in range(K):
            M[i,:] = tau_list**i#Generating tau^i
            T[i,0] = t**i#Generating t^i
    else:#For bandlimited models
        for i in range(K):
            if i % 2 == 1:#if odd
                f = (i+1)/2
                M[i,:] = np.sin(omega*f*tau_list)#Generate using Proposition 2
                T[i,0] = np.sin(omega*f*t)
            else:
                f = i/2
                M[i,:] = np.cos(omega*f*tau_list)
                T[i,0] = np.cos(omega*f*t)
    weights = np.matmul(np.linalg.inv(M), T)#W=(M^(-1))*T
    return weights

def G_t(G, w, mode):#Generating the time-varying gramians in terms of positive semidefinite gramians
    K = np.shape(w)[0]#The x dimension of the weights(2 dimensions for the bandlimited models)
    if mode:#For polynomial models
        G_tot = w[0]*G[0]
        for k in range(K-1):
            G_tot += w[k+1]*G[k+1]#Use proposition 1
    else:#For bandlimited models
        G_tot = np.zeros(G[0].shape)
        for k in range(K):
            G_tot += np.multiply(w[k],G[k])#Use proposition 2
    return G_tot


def mean(V):#Finding the mean of an array
    N = len(V)
    Vm = 0
    for n in range(N):
        Vm += V[n]
    return Vm/N

def K(param):#defining the number of basis Gramians
    P = param.P
    mode = param.mode
    if mode == 1:#If polynomial Model,number is 2K+1
        K = 2*P+1
    else:#If bandlimited model,number is 4K+1
        K = 4*P+1
    return K

def randomAs(param,N):#Defining the A matrices from equation 7 and 10
    A = []
    d = param.d
    P = param.P
    mode = param.mode
    if mode == 1:#If polynomial trajectory model
        for p in range(P+1):
            A.append(np.random.randn(d,N))#We have random matrices from R{dxN}
    else:#For the bandlimited model
        j = (-1) ** (1/2)
        for p in range(P):
            Ar = np.random.randn(d,N)
            Aq = np.random.randn(d,N)#Defining Ap's
            A.append(Ar+j*Aq)
        A.append(np.random.randn(d,N))#B0
        for p in range(P):
            A.append(np.conj(A[P-p-1]))#Defining Bp's
    return A

def X_t(A,t,param,N):#Defining the trajectories
    d = param.d
    P = param.P
    omega = param.omega
    mode = param.mode
    j = (-1) ** (1/2)
    X = np.zeros((d,N))+0*j
    if mode == 1:#Defining Polynomial trajectories
        for p in range(P+1):
            X += np.multiply(t**p, A[p])#Using equation(7)
    else:
        for p in range(P+1):#Defininf bandlimited trajectories
            ci = np.exp(j*p*omega*t)
            X += np.multiply(ci,A[p+P])+np.multiply(np.conj(ci),A[P-p])#Using equation (10),modified into exponential realm
        X -= A[P]
    X = np.real(X)
    return X

def gram2x(G,d):#Converting gram matrix to X
    N = G.shape[0]
    [U,S,V] = np.linalg.svd(G,full_matrices=True)#Eigen value decomposition
    S = S ** (1/2)#Square root of S
    S = S[0:d]#Taking dimension as d
    X = np.matmul(np.diag(S),V[0:d])#X = diag(S)*V
    return X

def rankProj(G,param,N):#Rank proection of Gramians,using proposition 5
    d = param.d
    K = param.K
    for k in range(K):
        [U,S,V] = np.linalg.svd(G[k].value,full_matrices=True)#Eigen value decomposition of values of G
        S[d:N] = 0
        G[k] = np.matmul(np.matmul(U,np.diag(S)),V)#ans is U*diag(S)*V, where rank of S is d
    return G


def fourier2sine(A):#Fourier transform in terms of factors of sine
    L, d, N = np.shape(A)
    C = []
    P = (L-1)//2
    j = (-1)**(1/2)
    
    C.append(np.real(A[P]))
    for i in range(2*P):
        p = i//2+1
        if i % 2 == 1:
            C.append( np.real(A[p+P]+np.conj(A[p+P])) )
        else:
            C.append( np.real(j*(A[p+P]-np.conj(A[p+P]))) )
    return C

def sine2fourier(A):#Converting factors of sine into fourier transform
    L, d, N = np.shape(A)
    C = []
    P = (L-1)//2
    j = (-1)**(1/2)
    
    for p in range(-P,P+1):
        C.append(np.zeros((d,N)))
    C[P] = A[0]
    for p in range(P):
        C[P+p+1] = (-j*A[2*p+1]+A[2*p+2])/2
        C[P-p-1] = np.conj(C[P+p+1])
    return C

def T_t(t, param):#Generating T_t, generated from section Spectral factorization of Gramian
    d = param.d
    P = param.P
    omega = param.omega
    mode = param.mode
    eyed = np.eye(d)
    if mode == 1:#If polynomial model
        for p in range(P+1):#The matrix is t^p * I(d,d) concateated for all P
            Tp = (t**p)*eyed
            if p == 0:
                T = Tp
            else:
                T = np.concatenate((T, Tp), axis=1)
    else:#For bandlimited model
        for i in range(2*P+1):
            if i % 2 == 1:#Alternating sin and cos time variants, multiplied by I(d,d)
                f = (i+1)/2
                Tp = np.sin(omega*f*t)*eyed
            else:
                f = i/2
                Tp = np.cos(omega*f*t)*eyed
            if i == 0:
                T = Tp
            else:
                T = np.concatenate((T, Tp), axis=1)
    return T

def deconcat(A_, param):#Ensures Trajectory dimension is P for Polynomial model and 2*P for bandlimited models
    d = param.d
    P = param.P
    mode = param.mode
    A = []
    if mode == 1:#If polynomial model
        for i in range(P+1):
            A.append(A_[i*d:(i+1)*d,:])
    else:
        for i in range(2*P+1):#If bandlimited models
            A.append(A_[i*d:(i+1)*d,:])
    return A

def testError(param, tau_list, estim_times, G, A,N):#Calculating the error between the real EDM and the edm we generate
    P = param.P
    omega = param.omega
    mode = param.mode
    erro=0
    N_estim = len(estim_times)
    eo = np.zeros(N_estim)
    for i_t, t in enumerate(estim_times):
        weights = W(t,tau_list,param)#Generating required parameters
        G_tot = G_t(G,weights,False)
        D_G = gram2edm(G_tot,N,False)#generating our matrix
        X = X_t(A,t,param,N)#Generating the real matrix
        #X_nois = X + np.random.randn(d,N)
        eo[i_t] = np.linalg.norm(edm(X)-D_G,'fro') / np.linalg.norm(edm(X),'fro')#Output error
    return eo

def rotationXY(X,Y): # X is anchor, Y is not,Rotation Matrix
    N = X.shape[1]
    J = np.eye(N) - np.ones((N,N))/N
    XY = np.matmul(np.matmul(X,J),np.conjugate(Y.T))
    U,_,V = np.linalg.svd(XY,full_matrices = True)  #XJY' = UV' = R
    R = np.matmul(U,V)
    return R

def randomAnchor(M,N,n):#Generating random anchor indexes
    anchor_idx = np.zeros((M,N))
    idx = np.zeros(N)
    idx[0:n] = 1#Generating n random instances of 1 and then shuffling them to get the dimension same
    for m in range(M):
        np.random.shuffle(idx)
        anchor_idx[m,:] = idx
    return anchor_idx


def edm(X, Y=None):
    if Y is None:
        Y = X
    
    d = X.shape[0]
    m = X.shape[1]
    n = Y.shape[1]
    D = np.zeros((m, n))

    for i in range(d):
        D += (X[np.newaxis, i, :].T - Y[i, :])**2#EDM = (XT-Y)^2,newaxis upscales the dimension of an array, used to ensure that correct dimension gets arithmetically used
    return D
