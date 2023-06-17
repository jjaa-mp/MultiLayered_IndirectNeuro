# -*- coding: utf-8 -*-
"""
Functional HALS
@author: C. Hautecoeur

Functions to apply Hierarchical Alternating Least-Squares on functional data.

The functional data can be expressed as a linear combination of a finite number of basis elements : 
    A = Pi*B     where B are the coefficients and Pi are the basis elements considered. 
    Evaluating Pi at some points allows the discretization of these functional data.

The problem can be expressed as    min ||Y-A*X||  where elements in A are (discretization of) functions

Two possible case are distinguished: 
    the I-case (integral-case) where data Y are functions
    the S-case (sum-case) where data Y are discretization of functions
    
Polynomials and Splines are implemented, but any function with a known projection on the non-negative set can be considered (given that this function is linearly parametrizable).
"""

import numpy as np
#from PHALS import buildPoly
from SHALS import buildBS
import time
import scipy.integrate as si

"""
Tools
"""
#
# Initialization of the options
#
def initialization(Y, case, d, options):
    
    ftol = options.get("ftol",1e-10) #tolerance
        
    maxiter = options.get("maxiter",1e4) #maximum number of iterations
    
    if case=='S':
        assert isinstance(Y, np.ndarray), "In sum case, Y must be a matrix containing the data"
        m,n = Y.shape
        if d is None: 
            d=m
    if case=='I':
        assert isinstance(Y,list) and d is not None , "In integral case, Y must be a list of functions and d must be defined"
        n = len(Y)
        m=0
     
    # Normalization
    normB = options.get("normB",True)
    normX = options.get("normX",True)
    
    xp = options.get("xp",np.linspace(-1,1,m))
        
    return(ftol, maxiter, m, n, normB, normX, d, xp)

#
# Build needed matrices and projection functions
#
def buildfun(funtype,m,n,r,d,case,initB,options):
    iterX = None
    PP = None
    
    if funtype == "original":
        assert d==m, "To use discretized HALS, vectors in Y must be the same size than vectors in B=A. Expected %d, found %d" % (d,m)
        case = 'S' # make no sense to use integral case here
        iterX = 1+0.5*(1+(n+r)*(2*d-1)/(2*n*r))
        Pi,fun,optf,iterB = None,None,{},None
    elif funtype == "poly":
        (Pi,fun,optf,iterB,initB,PP) = buildPoly(m,n,r,d,case,initB,options)  
    elif funtype == "bs":
        (Pi,fun,optf,iterB,d,PP) = buildBS(m,n,r,d,case,options)   
    elif funtype == "other":
        (Pi,fun,optf,iterB) = buildOther(options)
    else:
        raise AttributeError('Function type not handled. Please provide a valid name or "other" if you define your own problem')  
        
    return(Pi,fun,optf,iterX,iterB,initB,d,PP)

#
# Initialize matrices B and X, copying initB and initX if they are not "None", putting random nonnegative numbers otherwise
# Check if initial matrices have the requested shape and if matrix X contains only nonnegative values (the check is not made for B)
# d,n,r are respectively the degree of the functions (number of basis elements),
# the number of observations and the rank of the factorization
#
def init(d,n,r,initB,initX,iterB,iterX,funtype):
    if initB is None:
        B =  np.random.rand(d,r)
    else:
        assert initB.shape==(d,r), "Initial value of matrix B is invalid : wrong shape. Expected (%d,%d), found (%d,%d)" %(d,r,initB.shape[0],initB.shape[1])
        B = np.copy(initB)
    if initX is None:
        X = np.random.rand(r,n) # positive as between 0 and 1
    else:
        assert initX.shape==(r,n), "Initial value of matrix X is invalid : wrong shape. Expected (%d,%d), found (%d,%d)" %(r,n,initX.shape[0],initX.shape[1])
        assert np.all(initX>=0), "Initial value of matrix X is invalid : negative values found."
        X = np.copy(initX)
        X[X==0] = 1e-16 # Replace 0 values by a small number to avoid division by 0
    
    assert iterB is None or iterB>=1, "iterB must be greater than 1 or None. Your iterB is %f. " %(iterB)
    if iterB is None:
        iterB = 1 +0.5*(1+((2*n-1)*(d+r)+2*d**2+2*d+n-1)/(2*d*r))
        if funtype == "original":
            iterB = 1 +0.5*(1+((2*n-1)*(d+r)+3*d+n)/(2*d*r))
    iterB = min(int(iterB),10) 
        
    assert iterX is None or iterX>=1, "iterX must be greater than 1 or None. Your iterX is %f. " %(iterB)
    if iterX is None:
        iterX = 1+ 0.5*(1+((2*d-1)*(d+r+n)+4*n+d-1)/(2*n*r))
        if funtype == "original":
            iterX = 1+ 0.5*(1+((2*d-1)*(r+n)+4*n+d-1)/(2*n*r))
    iterX = min(int(iterX),10)
    
    return(B,X,iterB,iterX)
 
    
# Build a matrix containing the integrals between each the functions in A and B
# A is a list of functions and B is a list of vectors OR functions
# M[i,j] = integral(A_i * B_j)
# xp are the discretization points (if B contains vectors) and kn the knots at which the functions have a discontinuity in their curvature (useful for splines)
def buildIntMat(A,B,xp=None,kn=None):    
    n = len(A)
    m = len(B)
        
    M = np.zeros((n,m))
    
    if type(B[0]) == np.ndarray: # The functions in B are known via discretization points
        C = np.array(B)
        mul = (C[:,1:] - C[:,:-1])/(xp[1:]-xp[:-1])
        mul = np.transpose(mul)
        
        for i in range(n):
            def fun(t):
                return A[i](t)*t
            for k in range(len(xp)-1):
                i1,_ = si.quadrature(fun,xp[k],xp[k+1],vec_func=True)
                i0,_ = si.quadrature(A[i],xp[k],xp[k+1],vec_func=True)    
                
                M[i,:] = M[i,:] + i1*mul[k,:] + i0*(C[:,k+1] - mul[k,:]*xp[k+1]) 
    else:
        for i in range(n):
            for j in range(m):
                def fun(t):
                    return A[i](t) * B[j](t)
                if kn is not None:
                    init = kn[0]
                    for k in kn[1:]:
                        M[i,j],_ = M[i,j] +  si.quadrature(fun,init,k,vec_func=True )
                        init=k
                else:
                    M[i,j],_ = si.quadrature(fun,-1,1,vec_func=True,miniter=n )
                    
    return M

#
# Initialize somme useful matrices and return the initial error
#
def initMat(Y,B,X,Pi,case,m,PP,xp,normB,optf={}):
    kn= optf.get("kn",None)
        
    if case == 'S':
        if Pi is not None:
            Pinv = np.linalg.pinv(Pi)
            Yinv = Pinv@Y
            PY = np.transpose(Pi)@Y
            if PP is None:
                PP = np.transpose(Pi)@Pi            
        else:
            Yinv = Y
            PY = Y
            PP = None
        
    elif case == 'I':
        if PP is None:
            PP = buildIntMat(Pi,Pi,xp,kn)
        PY = buildIntMat(Pi,Y,xp,kn)
        Yinv = np.linalg.inv(PP)@PY
    else:
        raise AttributeError('Case not handled. Please define a valid case ("S" or "I")')
    error = cost2(PY,B,X,m,PP)
    
    optnorm = {}
    if normB:
        optnorm['PP'] = PP
        if case == 'S':
            optnorm['scaling'] = np.sqrt(m/2)
    
    return  (Yinv,PY,PP,error,kn,optnorm)

#
# Initialize the informations contained in the log
#
def initlog(maxiter,t0,Y,B,X,Pi,xp,options,case,error,funtype):
        err = np.zeros(maxiter+1)
        times = np.zeros(maxiter+1)
        realerr = np.zeros(maxiter+1)
        
        times[0] = time.perf_counter_ns()-t0
        
        Yreal = options.get("realM",Y)
                
        err[0] = error
        if funtype != 'original':
            Pi1 = Pi.copy()
            if case == 'I':
                Pi1 = np.zeros((len(xp),len(Pi)))
                for i in range(len(Pi)):
                    Pi1[:,i] = Pi[i](xp)
            realerr[0] = np.linalg.norm(Yreal-Pi1@B@X)
        else:
            Pi1 = None
            realerr[0] = np.linalg.norm(Yreal-B@X)
        
        return (err, times, realerr, Yreal, Pi1)

#
# Update the stopping criteria
#
def StopCond(PY,m,B,X,PP,error,prevB,prevX):
    error2 = cost2(PY,B,X,m,PP)       
    de = abs(error-error2)/abs(error2)
    error = error2
        
    dx = np.linalg.norm(prevB-B)**2+np.linalg.norm(prevX-X)**2
    dx = np.sqrt(dx)
    prevB = np.copy(B)
    prevX = np.copy(X)
    
    
    if PP is not None:
        PB = PP@B
    else:
        PB = B
    grad1 = -np.transpose(B)@PY+np.transpose(B)@PB@X
    grad1[X<1e-16] = np.minimum(grad1[X<1e-16],np.zeros(grad1[X<1e-16].shape))
    grad1 = np.linalg.norm(grad1)
    ngrad = grad1
    
    return (de,dx,ngrad,error,prevB,prevX)

#
# Build the necessary matrices when the option "other" is chosen
#
def buildOther(options):
    fun = options.get("fun",None)
    optf = options.get("optf",{})
    
    if 'Pi' in options:
        Pi = options['Pi']
    else: 
        raise AttributeError("You must provide at least martix/list Pi to define your own functional HALS")

    return  (Pi,fun,optf)
    
"""
F-HALS functions
"""
#
# Compute the cost function
# PY is the multiplication of the transpose of matrix Pi with the data (Pi'@Y) it is the matrix Z in the paper
# B is the coefficient matrix and X is the mixing one
# m is the number of discretization points in Y (if m==0, we are in integral case)
# PP is equal to Pi'Pi, where Pi is the basis matrix of the choosen function (if None, original case is considered)
# PP is called "M" in paper
# 
def cost2(PY,B,X,m,PP=None):
    BX = B@X
    if PP is None:
        e1 = np.sum(BX*BX)
    else:
        e1 = np.sum(BX*(PP@BX))
    e2 = 2*np.sum(PY * BX)
    if m==0:
        return e1-e2
    else:
        return (e1-e2)*2/m
    
#
# Update of mixing matrix (matrix X). This matrix does not contain functional data but usual vectors
#
# PY contains the modified data (matrix Z in paper )
# B is the matrix of coefficients of the functional basis
# X is the mixing matrix to update
# PP is equal to Pi'Pi in S-case, or to its integral in I-case (matrix M in paper), if none the identity matrix is considered (original HALS)
# nbup are the (maximum) number of updates to perform before to switch to the update of B (default 1) 
# norm is True if the matrix must be normalized (default False)
#
def updateX(PY,B,X,PP=None,nbup=1,norm=False):
    x1 = X.shape[0]

    BT = np.transpose(B)
    BY = BT@PY
    if PP is None:
        BB = BT@B
    else:
        BB = BT@PP@B
    
    if nbup>1:
        controlN = 0 # to stop earlier iterations if necessary
        actualN = 1
        X0 = np.copy(X)
        prevX = np.copy(X)
        
    for i in range(0,nbup):
        if(nbup<=1 or actualN>0.1*controlN):
            for j in range(0,x1):
                C = BY[j,:] - BB[j,:]@X
                xj = C/BB[j,j] + X[j,:] 
                xj[xj<0]=1e-16
                
                X[j,:] = xj # update
                    
            if i<nbup-1: # not the last update
                actualN = np.linalg.norm(X-prevX)
                controlN = np.linalg.norm(X-X0)
                prevX = np.copy(X)
        else:
            break
            
    for j in range(x1):
        dX = 1e-16 + np.linalg.norm(X[j,:])
        B[:,j] = B[:,j]*dX
        X[j,:] = X[j,:]/dX
    return X,B

#
# Update of coefficient matrix (matrix B). This matrix contains coefficient of functions.
#
# Yinv contains the modified data ( M^(-1) Z in paper)
# B is the matrix of coefficients of the functional basis to update
# X is the mixing matrix 
# fun is the projection function. If None, the projection is made over the nonnegative coefficients (default)
# optf are the options to provide to function "fun"
# nbup are the number of updates to perform before to switch to the update of X (default 1) 
# norm is True if functions must be normalized (default False)
# if norm is True, optnorm must contain matrix PP ( = matrix M in paper) it can also contain a "scaling" (B is normalized by norm=(B@PP@B)^1/2 and then multiplied by "scaling")
#
def updateB(Yinv,B,X,fun=None,optf=None,nbup=1,norm=False,optnorm=None):    
    b2 = B.shape[1]
    
    XT = np.transpose(X)    
    YX = Yinv@XT
    XX = X@XT
    
    if nbup>1:
        controlN = 0 # to stop earlier iterations if necessary
        actualN = 1
        B0 = np.copy(B)
        prevB = np.copy(B)
    if norm:
        assert optnorm is not None and 'PP' in optnorm, "To normalize functions, optnorm must contain matrix PP "
        PP = optnorm['PP']
        scaling = optnorm.get("scaling",1)
            
    for i in range(0,nbup):
        if(nbup<=1 or actualN>0.1*controlN):
            for j in range(0,b2):
                C = YX[:,j]-B@XX[:,j]
                bj = C/XX[j,j] + B[:,j]
                
                if fun is None:
                    bj[bj<0]=1e-16
                else:
                    bj = fun(bj,optf)
                
                B[:,j] = bj # update
                
            if i<nbup-1: # not the last update
                actualN = np.linalg.norm(B-prevB)
                controlN = np.linalg.norm(B-B0)
                prevB = np.copy(B)
        else:
            break
    if norm:
        for j in range(0,b2):
            if PP is None:
                dB = np.linalg.norm(B[:,j])
            else:
                dB = np.sqrt(np.transpose(B[:,j])@PP@B[:,j])/scaling
            B[:,j] = B[:,j]/dB
            X[j,:] = X[j,:]*dB
    return B,X



#
# Apply HALS algorithm
# Y is the data to approximate, or a list of callable function representing the data
# r is the rank of factorization
# d is the size of the parametrization basis (degree of polynomial +1 for example), or the degree of the splines (only 3 is supported in this case)
# case is either S (sum case) or I (integral case), by default S 
# funType is the choosen function type of elements recovered in A. 
#   - original : basic HALS (non functional case)
#   - poly : (by default) nonnegative polynomials (Chebyshev basis)
#   - bs : nonnegative B-splines
#   - other : non-predifined function. User has to specify needed matrices (Pi) in options, as well as the projection
# initB is the initial matrix B. random matrix with nonnegative numbers if not specified
# initX is the initial matrix B. random matrix with nonnegative numbers if not specified
# 
# options : other needed parameters. Depends on the choosen function. Must contains matrix or list Pi and projection (at least) for non pre-builded functions in this case it can also contain:
#   - Pi is matrix of discretized basis functions (S-case) or the list of basis function (I-case)
#   - fun is the projection function (for good convergence properties, should be define such as the function cannot be zero everywhere)
#   - optf are the options for this projection
#   - iterB is the number of time the matrix B is computed before to change to X (iterX is computed automatically). This value should depend on the choosen function and must be greater than 1 or None. If None, the value for iterB for a None function is used (projection over the nonnegative coefficients)
#   - other options depends on the choosen function type => See builder for these functions 
#   For all :
#       - xp : discretization points
#       - normB : boolean to normalize or not matrix B (default: true)
#       - normX : boolean to normalize or not matrix X (default: true)
#       - compuTime : boolean to indicate if the computation time must be returned or not
#       - log: if true returns more informations about the computations as well as the "real cost" if matrix realM is given
#       - realM : original matrix used to produce the data
#       - ftol : tolerence considered for stopping criteria  (default 1e-07)
#       - maxiter : maximal number of iterations (default 5000)
def fHALS(Y, r, d=None, case='S', funtype='poly', initB=None, initX=None, options={}):
    t0 = time.perf_counter_ns()
    
    #Initialize the options
    (ftol, maxiter, m, n, normB, normX, d, xp) = initialization(Y, case, d, options)
        
    #Build needed matrices and projection functions
    (Pi,fun,optf,iterX,iterB,initB,d,PP) = buildfun(funtype,m,n,r,d,case,initB,options)
    
    #Initialization of matrices B and X, as well as the number of inner iterations (iterB/iterX)
    (B,X,iterB,iterX) = init(d,n,r,initB,initX,iterB,iterX,funtype)
    prevB = np.copy(B)
    prevX = np.copy(X)  
        
    #Other useful matrices         
    (Yinv,PY,PP,error,kn,optnorm) = initMat(Y,B,X,Pi,case,m,PP,xp,normB,optf)
    
    #Stopping variables (more than just checking error evolution, also look at updates evolution and gradient) 
    de = 1 # error difference
    dx = 1 # updates difference
    ngrad = max(m,10) # gradient value
    iter = 0 # number of iterations
    
    log = 'log' in options and options['log']
    if log:
        (err, times, realerr, Yreal, Pi1) = initlog(maxiter,t0,Y,B,X,Pi,xp,options,case,error,funtype)
        t0 = time.perf_counter_ns()
            
    vdx = np.linalg.norm(prevB)**2+np.linalg.norm(prevX)**2
    if PP is not None:
        PB = PP@B
    else:
        PB = B
    grad1 = -np.transpose(B)@PY+np.transpose(B)@PB@X
    grad1 = np.linalg.norm(grad1)
    
    t = time.perf_counter_ns()
    ind = 0
    
    #Iterations
    while(de>ftol and iter<maxiter and dx>ftol*(ftol+np.sqrt(vdx)) and ngrad>ftol*grad1):
        (X,B) = updateX(PY,B,X,PP, nbup = iterX, norm = normX)
        
        (B,X) = updateB(Yinv,B,X,fun,optf, nbup = iterB, norm = normB,optnorm = optnorm)
        
        (de,dx,ngrad,error,prevB,prevX) = StopCond(PY,m,B,X,PP,error,prevB,prevX)
        
        vdx = np.linalg.norm(prevB)**2+np.linalg.norm(prevX)**2
        
        if log and (iter%10==0 or iter==1):
            times[ind+1] = times[ind] + time.perf_counter_ns() - t0
            err[ind+1] = error
            if funtype != 'original':
                realerr[ind+1] = np.linalg.norm(Yreal-Pi1@B@X)
            else:
                realerr[ind+1] = np.linalg.norm(Yreal-B@X)
            ind = ind+1
            t0 = time.perf_counter_ns()
        iter = iter+1
        
    if "compuTime" in options and options["compuTime"]:
        compu = (time.perf_counter_ns()-t)/1e9
        return(B,X,error,iter,compu)
    if log:
        return(B,X,err[0:ind+1],realerr[0:ind+1]/np.linalg.norm(Yreal), times[0:ind+1]/1e9)
    return(B,X,error,iter)
