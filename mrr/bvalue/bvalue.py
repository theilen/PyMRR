import numpy as np
import math

import bvalue

#b : s/mm^2 = 10^6 s/m^2
#G : mT/m = 10^-3 T/m
#d : ms = 10^-3 s
#D : ms = 10^-3 s
#dD : ms = 10^-3 s

gamma = 42.57e6*2.*math.pi #[Hz/T]

Dd_flag = True

print 'Using D = Delta-delta for calculations.'


#def set_gamma(g):
#    'Set gamma to a new value (in SI-units [Hz/T]'
#    bvalue.gamma = float(g)
    
def switch_Delta(flag):
    '''
    Switch between Calculations using 'delta-delta' (space between gradients)
    and 'delta' (space between the center of the gradients)
    '''
    if flag == 'delta':
        bvalue.Dd_flag = False
        print 'Using D = Delta for calculations.'
    elif flag == 'delta-delta':
        bvalue.Dd_flag = True
        print 'Using D = Delta-delta for calculations.'
    else:
        raise ValueError('Wrong flag provided!')


def D(b,G,d):
    '''
    Calculate Delta (or Delta-delta, depending on the configuration)
    in [ms].
    
    Provide all data in common units, e.g.:    
    b : [s/mm]
    G : [mT/m]
    d : [ms]
    '''
    b = b*1e6
    G = G*1e-3
    d = d*1e-3
    if bvalue.Dd_flag:
        return calc_Dd(b,G,d)*1e3
    else:
        return calc_D(b,G,d)*1e3
    
def d(b,G,D, all_solutions=False):
    '''
    Calculate delta in [ms]. D is either the space between the gradients 
    or the space between the center of the gradients, depending on the 
    current configuration.
    
    Mathematically there are up to three possible solutions. If verbose is set
    to True (not default), all three are returned, otherwise the physically 
    meaningful solution is guessed and returned.
    
    Provide all data in common units, e.g.:    
    b : [s/mm]
    G : [mT/m]
    D : [ms]
    '''
    b = b*1e6
    G = G*1e-3
    D = D*1e-3
    if bvalue.Dd_flag:
        sol = bvalue.calc_d_Dd(b,G,D,verbose=False)
    else:
        sol = bvalue.calc_d(b,G,D)
    if all_solutions:
        return [ d_*1e3 for d_ in sol ]
    else:
        return min([ d_ for d_ in sol if d_>0 ])*1e3
    
def G(b,d,D):
    '''
    Calculate the gradient strength G in [mT/m].
    
    Provide all data in common units, e.g.:    
    b : [s/mm]
    d : [ms]
    D : [ms]
    '''
    b = b*1e6
    d = d*1e-3
    D = D*1e-3
    if bvalue.Dd_flag:
        return bvalue.calc_G_Dd(b,d,D)*1e3
    else:
        return bvalue.calc_G(b,d,D)*1e3
        
def b(G,d,D):
    '''
    Calculate the b-value in [s/mm^2].
    
    Provide all data in common units, e.g.:    
    G : [mT/m]
    d : [ms]
    D : [ms]
    '''
    G = G*1e-3
    d = d*1e-3
    D = D*1e-3
    if bvalue.Dd_flag:
        return bvalue.calc_b_Dd(G,d,D)*1e-6
    else:
        return bvalue.calc_b(G,d,D)*1e-6
        
        
        




def calc_G(b,d,D):
    'Calculate G in [T/m] using Delta'
    g = bvalue.gamma
    return math.sqrt(b/(g**2.*d**2.*(D-d/3.)))

def calc_G_Dd(b,d,Dd):
    'Calculate G in [T/m] using Delta-delta'
    return calc_G(b,d,Dd+d)

def calc_b(G,d,D):
    'Calculate b in [s/m^2] using Delta'
    g = bvalue.gamma
    return g**2.*G**2.*d**2.*(D-d/3.)

def calc_b_Dd(G,d,Dd):
    'Calculate b in [s/m^2] using Delta-delta'
    return calc_b(G,d,Dd+d)

def calc_D(b,G,d):
    'Calculate Delta in [s]'
    g = bvalue.gamma
    return b/(g**2.*G**2.*d**2.) + d/3.

def calc_Dd(b,G,d):
    'Calculate Delta-delta in [s]'
    return calc_D(b,G,d) - d

def calc_d(b,G,D,verbose=False):
    '''
    Calculate d in [s] using Delta.
    
    Returns three possible solutions if verbose == True and a list of 
    meaningful solutions otherwise.
    '''
    g = bvalue.gamma
    phi = math.acos(1.-3.*b/(2.*g**2*G**2.*D**3.))
    d1 = 2.*D*math.cos(phi/3.) + D
    d2 = D - 2.*D*math.cos((phi+math.pi)/3.)
    d3 = D - 2.*D*math.cos((phi-math.pi)/3.)
    if verbose: return d1,d2,d3
    return [d_ for d_ in [d1,d2,d3] if d_ > 0 and d_ < D]

def calc_d_Dd(b,G,Dd):
    '''
    Calculate d in [s] using Delta-delta.
    
    Returns up to three positive possible solutions.
    '''
    g = bvalue.gamma
    p = -0.25*Dd**2.
    q = Dd**3/8. - 3.*b/(4.*g**2*G**2)
    D = p**3 + q**2.
    if D > 0:
        #if verbose: print 'D =',D,'> 0'
        d1 = (-q+math.sqrt(D))**(1/3.) + (-q-math.sqrt(D))**(1/3.)-0.5*Dd
        d2 = d3 = None
    elif D < 0:
        #if verbose: print 'D =',D,'< 0'
        phi = math.acos(-q/math.sqrt(abs(p)**3))
        d1 = 2.*math.sqrt(abs(p))*math.cos(phi/3.)-0.5*Dd
        d2 = -2.*math.sqrt(abs(p))*math.cos((phi+math.pi)/3.)-0.5*Dd
        d3 = -2.*math.sqrt(abs(p))*math.cos((phi-math.pi)/3.)-0.5*Dd
    return [d_ for d_ in [d1,d2,d3] if d_ > 0]# and d_ < Dd])

def create_list(Dd, b_0, d_0, length=50, dd = 1e-3):
    data = np.empty((length,4))
    data[:,0] = np.linspace(d_0, d_0+dd*(length-1), length)
    G = calc_G_Dd(b_0,d_0,Dd)
    data[0,0] = d_0
    data[0,1] = b_0
    for i in range(1,length):
        data[i,1] = calc_b_Dd(G,data[i,0],Dd)
    for i in range(length):
        data[i,2] = round(data[i,1],-6)
        data[i,3] = calc_d_Dd(data[i,2],G,Dd)
    data[:,[0,3]] *= 1000.
    data[:,1:3] /= 1000000.
    return {'G':G,'Dd':Dd},data
