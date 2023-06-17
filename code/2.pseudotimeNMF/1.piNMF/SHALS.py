# -*- coding: utf-8 -*-
"""
S-HALS
@author: C. Hautecoeur

Contains tool functions to do functional HALS using B-Splines

"""

import numpy as np
from scipy.interpolate import BSpline
import mosek.fusion as mos
import scipy.integrate as si

"""
General tools
"""

# Returns a random initialisation of matrix B with rank r and B-Splines of degree d using interior knots kn
# This random initialisation is made using positive coefficients
# sparsity is the sparsity level of the initialisation. If negative, create a dense matrix


def randInitBS(d,kn,r,sparsity = -1):
    s = kn.size + d-1  ## jm: ?
    if sparsity<0 or sparsity>s:
        B = np.random.rand(s,r) #dense matrix
    else:
        index = np.random.randint(0, s, size=(sparsity,r))
        B = np.zeros((s,r))
        for i in range(0,r):
            B[index[:,i],i] = np.random.rand(sparsity)
    return B     

# Returns a random initialisation of matrix B with rank r and B-Splines of degree d using interior knots kn
# This random initialisation is not imposed to have positive coefficients
# and is found as the projection of a random spline (slower than randInitBS)

def randNonNegBS(d,kn,r):
    B = np.random.randn(len(kn)+d-1,r)
    
    S = buildS(d,kn, discrete=False)
    M = buildM(S,kn) 
    
    optf = {}
    optf["T"] = buildT(kn)

    L = np.transpose(np.linalg.cholesky(M))
    optf["L"] = L
    
    for i in range(r):
        B[:,i] = findNN(B[:,i],optf)
        
    return B

# Returns the B-Spline basis
# kn are the interior knots (whith the borders included only once)
# d is the degree of the splines
# xp are the discretization points. If not given, equals to kn
# If discrete is True, the basis elements are discretized.
def buildS(d,kn,xp=None,discrete=True):
    if discrete and xp is None:
        xp = kn
    s = kn.size + d-1
    if discrete:
        S = np.zeros((xp.size,s))
    else:
        S = []
        
    # Repeat the borders d times more (d+1 occurences in total) 
    # to have a closed interval
    kn2 = np.zeros(kn.size + 2*d)
    kn2[0:d]=kn[0]
    kn2[d:kn.size+d] = kn
    kn2[kn.size+d:] = kn[-1]

    for i in range(s):
        #coefficient of the spline (zero vector with a 1 at the ith element)
        c = np.zeros(s)
        c[i] = 1
        
        spl = BSpline(kn2,c,d, extrapolate=True)
        
        if discrete:
            S[:,i] = spl(xp)
        else:
            S.append(spl)
    return S

# Build matrix M for the integral case
# S contains functions and kn are knots at which the curvature of the function changes (interior knots of splines)
# M[i,j] is equal to the integral of S_i(t)*S_j(t) over interval [kn[0], kn[-1]]
def buildM(S,kn):  
    n = len(S)
    
    M = np.zeros((n,n))

    for i in range(n):
        for j in range(i,min(i+4,n)):
            def fun(t):
                return S[i](t) * S[j](t)
            
            deb = max(0,j-3)
            end = min(i+2,len(kn)+1)
            init = kn[deb]
            for k in kn[deb+1:end]:
                M[i,j],_ = M[i,j] +  si.quadrature(fun,init,k,vec_func=True )
                init=k
            M[j,i] = M[i,j]
    return M

"""
Tools for projection
"""
# Build the linking matrix between the coefficients of the second order polynomials and the cubic one (see paper)
# (c3x^3 c2x^2 c1x c0)' = Q (a1 b1 c1 a2 b2 c2)'
def buildQ():
    Q = np.zeros((4,6))
    
    b = np.sqrt(2)
    
    Q[0,0] = 1
    Q[0,3] = -1
    Q[1,2] = -b
    Q[1,3] = 1
    Q[1,5] = b
    Q[2,1] = 1
    Q[2,4] = -1
    Q[2,5] = -b
    Q[3,4] = 1
    
    return Q

# Construct the matrices of transformation for each considered interval
# returns a list with the appropiate matrix at index i (0 intervals are not considered)
# we consider here only cubic splines
# kn are the knots of the spline. They must be defined in increasing order
#
# [ci3x^3 ci2x^2 ci1x ci0]' = N_i [s_i-3 s_i-2 s_i-1 s_i]'
# returns T = N_i^(-1) Q (see paper)
def buildT(kn):
    Q = buildQ()
    deg=3 # only case supported
    d = np.zeros(kn.size+deg)
    d[deg-1:-deg+1] = kn[1:]-kn[0:-1] # Intervals between knots
    
    assert all(d>=0), "the knots must be defined in increasing order" 
    
    T = []
    
    for i in range(deg-1, deg-2+kn.size): # number of intervals
        if d[i] >0: # Do not consider 0 intervals (no constraint in this case)
            Ni = np.zeros((deg+1,deg+1))
            # function defined by s_i
            Ni[0,deg] = d[i]**2/((d[i] + d[i+1]) * (d[i] + d[i+1] + d[i+2]) )
            
            # function defined by s_i-1
            det1 = (d[i-1]+d[i])*(d[i-1]+d[i]+d[i+1])
            det2 = (d[i+1]+d[i])*(d[i-1]+d[i]+d[i+1])
            det3 = (d[i+1]+d[i])*(d[i+2]+d[i]+d[i+1])
            
            Ni[0,2] = -d[i]**2 * (1/det1 + 1/det2 + 1/det3 ) #x^3
            Ni[1,2] = d[i] * ( (d[i]-2*d[i-1])/det1 + (d[i+1] + d[i] - d[i-1])/det2 + (d[i] + d[i+1] + d[i+2])/det3 ) #x^2
            Ni[2,2] = d[i-1] * ( (2*d[i] - d[i-1])/det1 + (d[i+1] + d[i])/det2 ) #x
            Ni[3,2] = d[i-1]**2/det1
            
            # function defined by s_i-2
            det1 = (d[i-1]+d[i])*(d[i-1]+d[i]+d[i-2])
            det2 = (d[i-1]+d[i])*(d[i-1]+d[i]+d[i+1])
            det3 = (d[i+1]+d[i])*(d[i-1]+d[i]+d[i+1])
            
            Ni[0,1] = d[i]**2 * (1/det1 + 1/det2 + 1/det3) #x^3
            Ni[1,1] = d[i] * ( (d[i-2] + d[i-1] - 2*d[i])/det1 - (2*d[i] + d[i+1] - d[i-1])/det2 -2*(d[i] + d[i+1])/det3 ) #x^2
            Ni[2,1] = (d[i]*(d[i] - 2*d[i-2] - 2*d[i-1])/det1 + (-d[i-1]*(d[i]+d[i+1]) + d[i]*(d[i]+d[i+1]-d[i-1]))/det2 + (d[i] + d[i+1])**2/det3 ) #x
            Ni[3,1] = d[i]*(d[i-2]+d[i-1])/det1 + d[i-1]*(d[i] + d[i+1])/det2
            
            # function defined by s_i-3
            mul = d[i]**2/((d[i] + d[i-1])*(d[i] + d[i-1] + d[i-2]))
            Ni[0,0] = -mul
            Ni[1,0] = mul*3
            Ni[2,0] = mul*(-3)
            Ni[3,0] = mul
        
            Ti = np.linalg.inv(Ni) @ Q
            T.append(Ti)
        
    return T

""" Projection """

# Exact projection
# Find the closest non-negative spline to f, using matrices in T 
# only degre 3 splines are accepted
# optf contains:
#   T the matrices of transformation for each considered interval
#   L the "norm" matrix. See paper for more details
def findNN(f,optf):
    assert 'T' in optf and 'L' in optf, "You need to define matrix T and matrix L to use exact projection."
    T = optf['T']
    L = optf['L']
    
    if np.all(f>=0): # The spline is always nonnegative if all its coefficients are nonnegative
        return f
    norm = np.max(np.abs(f)); f2 = f/norm # avoid to deal with too large numbers
    
    f2[0] = f2[0]-1e-16 # to avoid obtaining 0-projection
    
    #Optimization 
    M = mos.Model()
    
    alpha = M.variable("alpha",mos.Domain.unbounded(f.size)) # the coefficients of the projection
    
    for i in range(len(T)): # constraints of positivity
        a = M.variable("a"+str(i),mos.Domain.inRotatedQCone(3))
        b = M.variable("b"+str(i),mos.Domain.inRotatedQCone(3))
        
        cons1 = alpha.pick(np.arange(i,i+4,dtype=np.int32))
        vector = mos.Expr.vstack(a,b)
        Ti = mos.Matrix.dense(T[i])
        cons2 = mos.Expr.mul(Ti,vector)
        M.constraint(mos.Expr.sub(cons1,cons2),mos.Domain.equalsTo(0.))

    L1 = mos.Matrix.sparse(L)
    l = M.variable("l",mos.Domain.inQCone(f.size+1)) #[t,u]
    t = l.pick(np.array([0],dtype=np.int32))
    u = l.pick(np.arange(1,f.size+1,dtype=np.int32))
    
    cons = mos.Expr.sub(u,mos.Expr.mul(L1,alpha))
    M.constraint(cons,mos.Domain.equalsTo(-L@f2))
    
    M.objective(mos.ObjectiveSense.Minimize,t)
    
    M.solve()
    
    if(M.getPrimalSolutionStatus() != mos.SolutionStatus.Optimal): # in case of non convergence
        print("Optimization did not work, use nonnegative coefficients instead")
        print(M.getPrimalSolutionStatus())
        print(f)
        g = f.copy()
        # If the optimization did not work, project the coefficients over the nonnegative set instead
        g[g<0] = 0
        if g[0]==0:
            g[0] = g[0]+1e-16
        M.dispose()
        return (g)
    
    alpha = alpha.level()
    
    M.dispose()
    alpha[0] = alpha[0] + 1e-16
    return alpha*norm

# Projection using nonnegative coefficients
# Find the closest nonnegative spline to f
# only degre 3 splines are accepted
# optf contains:
#   L the "norm" matrix. See paper for more details
def findCoefNN(f,optf):
    assert  'L' in optf, "You need to define matrix L to use projection with nonnegative coefficients."
    L = optf['L']
    
    if np.all(f>=0): # The spline is always nonnegative if all its coefficients are nonnegative
        return f
    norm = np.max(np.abs(f))
    f2 = f/norm # avoid to deal with too large numbers
    
    f2[0] = f2[0]-1e-16
    M = mos.Model()
    
    # Definition of the variables
    alpha = M.variable("alpha",mos.Domain.greaterThan(np.zeros(f.size)))

    L1 = mos.Matrix.sparse(L)
    l = M.variable("l",mos.Domain.inQCone(f.size+1)) #[t,u]
    u = l.pick(np.arange(1,f.size+1,dtype=np.int32))
    
    cons = mos.Expr.sub(u,mos.Expr.mul(L1,alpha))
    M.constraint(cons,mos.Domain.equalsTo(-L@f2))
    
    t = l.pick(np.array([0],dtype=np.int32))
    M.objective(mos.ObjectiveSense.Minimize,t)
    
    M.solve()
    
    if(M.getPrimalSolutionStatus() != mos.SolutionStatus.Optimal):
        print("Optimization did not work, use nonnegative coefficients instead")
        print(M.getPrimalSolutionStatus())
        print(f)
        g = f.copy()
        # If the optimization did not work, project the coefficients over the nonnegative set instead
        g[g<0] = 0
        if g[0]==0:
            g[0] = g[0]+1e-16
        M.dispose()
        return (g)
    
    alpha = alpha.level()
    
    M.dispose()
    alpha[0] = alpha[0] + 1e-16
    return alpha*norm

"""
Build the matrices and function needed by F-HALS following the options given

n is the number of observations
m the number of discretization points (S case)
r is the rank of factorization
d is the degree of the splines (only d=3 is possible when using exact projection)

Case is either 'I' for integral case or 'S' for sum case

OPTIONS:
    - kn : interior knots for splines (default : linspace(0,1,m)). Determine the basis used!
    - xp : discretization points for splines 
    - heuri : (0,1). If 0, basic heuristic projecting the coefficients on the nonnegative set. If 1 exact projection. Default 1
"""
def buildBS(m,n,r,d,case,options):
    kn = options.get("kn", np.linspace(0,1,m))
        
    xp = options.get("xp", np.linspace(kn[0],kn[-1],m))
    
    s = kn.size + d-1
        
    optf = {}
    optf['kn']=kn
    
    S = buildS(d,kn,xp)
    optf['S'] = S
    if case=='I':
        S = buildS(d,kn, discrete=False)
        PP = buildM(S,kn)    
    elif case=='S':
        PP = np.transpose(S)@S
    else:
        raise AttributeError('Case not handle. Please define a valid case ("S" or "I")')
        
    L = np.transpose(np.linalg.cholesky(PP))
    optf['L'] = L
        
    # Choice of projection
    heuri = options.get("heuri",0) #exact projection by default
        
    if heuri == 0: # Basic heuristic
        fun = None #the projection is the default projection in this case
        iterB = None
    elif heuri == 1: # Exact projection
        T = buildT(kn)
        optf['T'] = T
        fun = findNN
        iterB = 1 # the projection is time-consuming
    else:
        raise AttributeError('Case not handle. Please define a valid heuristic (0 or 1)')
        
    return S,fun,optf,iterB,s,PP
