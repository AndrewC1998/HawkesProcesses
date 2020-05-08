import numpy as np
import matplotlib.pyplot as plt

def sample_d_im(u,l,b):
    if u >= 1 - np.exp(-l/b):
        return np.inf
    return - 1/b * np.log(1 + b/l * np.log(1-u))
sample_d_im = np.vectorize(sample_d_im)

def HawkesSim(T, mu, b, Y, lamb0 = None, plot = True):
    '''
    Simulates M-dimentional Hawkes point process

    Parameters
    ----------
    T : Time period covered by the Simulator
    mu : Base intensities for each point process
    b : Cross-perseverance of the effect of an event on each point process.
            High b => Short-lived effect
    Y : random generator for jumps in intensity due to an event
    lamb0 : Initial cross-intensities

    Returns
    -------
    r : Event times for all point processes
    causes : Cause of the event at the corresponding event time in r.
                -1 => Immigrant
    effects : Label (0 to M-1) for event at the corresponding event time in r

    '''
    
    M = len(mu)
    r = [0,]
    causes = []
    effects = []
    if lamb0 is not None:
        lamb = lamb0.copy()
    else:
        lamb = np.zeros(shape = (M,M))
    
    while r[-1] < T:
        
        U = np.random.uniform(size = (M,M))
        
        # Sample waiting times for children
        d_im = sample_d_im(U, lamb, b)
        min_index_child = np.argmin(d_im)
        i = min_index_child//M
        m = min_index_child - M*i
        min_d_child = d_im[i,m]
        
        # Sample waiting times for immigrants
        d0 = np.random.exponential(scale = 1/mu)
        min_index_immigrant = np.argmin(d0)
        min_d_immigrant = d0[min_index_immigrant]
        
        # Find lowest waiting time
        if min_d_immigrant < min_d_child:
            d = min_d_immigrant
            i = -1
            m = min_index_immigrant
        else:
            d = min_d_child
        
        # Logging
        r.append(r[-1]+d)
        causes.append(i)
        effects.append(m)
        
        # Update intensities
        lamb = lamb * np.exp( -b * d )
        lamb[m,:] += Y(m, M)
    
    r = np.asarray(r)[1:-1]
    causes = np.asarray(causes)[:-1]
    effects = np.asarray(effects)[:-1]
    
    if plot:
        plotHawkes(M, r, effects)
    
    return r, causes, effects

def plotHawkes(M, r, effects):
    if len(r)<3:
        print('Not enough points.')
        return
    
    ax = plt.axes()
    for i in range(M):
        ith_process = r[effects == i]
        plt.step(ith_process, np.arange(1, len(ith_process)+1), label = i)
    if len(r)<1e3:
        plt.xticks(r)
        ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax.tick_params(length = 8, width = 1, direction = 'inout')
    plt.legend()
    plt.show()

M = 6
b = np.random.normal(loc = 0.5, size = (M,M))**2
mu = 1+np.random.normal(size = M)**2
def generateY(m, M):
    return np.min(b)

r, c, e = HawkesSim(10, mu, b, generateY)
