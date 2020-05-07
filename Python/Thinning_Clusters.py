import numpy as np
import matplotlib.pyplot as plt

def PoissonByThinning(T, intensity, M):
    t = 0
    P = []
    
    while t < T:
        E = np.random.exponential(M)
        t += E
        U = np.random.uniform(high = M)
        
        if t < T and U < intensity(t):
            P.append(t)
    
    return np.asarray(P)

def HawkesByThinning(T, v, a, b, l0 = 0):
    t = 0
    P = []
    
    l = v + l0
    
    while t < T:
        
        M = l
        E = np.random.exponential(M)
        t += E
        U = np.random.uniform(high = M)
        
        l = v + (l - v) * np.exp(-b*E)
        
        if t < T and U <= l:
            P.append(t)
            l += a
    
    return np.asarray(P)

def HawkesByClusters(T, v, a, b):
    P = []
    k = np.random.poisson(v*T)
    C = np.random.uniform(high = T, size = k)
    D = np.random.poisson(a/b, size = k)
    
    for i, D_i in enumerate(D):
        P.extend( C[i]+ np.random.exponential(b, size = D_i) )
    
    P.extend(C).sort()
    
    return np.asarray(P)

def plotHawkes(P, v, a, b, l0 = 0):
    if len(P)<2:
        return 'Not enough points.'
    resolution = np.min(P[1:] - P[:-1])/3
    
    times = [0,]
    intensities = [v+l0, ]
    
    i = 0
    n = len(P)
    while i < n:
        nextTimePoint = times[-1] + resolution
        
        if nextTimePoint > P[i]:
            intensities.append( a + v+ (intensities[-1] - v) * np.exp(-b*(P[i]-times[-1])) )
            times.append(P[i])
            i += 1
            
        else:
            times.append(nextTimePoint)
            intensities.append( v+ (intensities[-1] - v) * np.exp(-b*resolution) )
            
    
    plt.plot(times, intensities)
    plt.show()
