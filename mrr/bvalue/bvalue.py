# -*- coding: utf-8 -*-
"""
Calculate gradient parameters for Stejskal Tanner unipolar gradient forms
from bvalues and vice versa.

@author: Sebastian Theilenberg
"""

__version__ = '1.0'
# $Source$


import math


# b : s/mm^2 = 10^6 s/m^2
# G : mT/m = 10^-3 T/m
# d : ms = 10^-3 s
# D : ms = 10^-3 s
# dD : ms = 10^-3 s

# Gamma
gamma = 42.57e6*2.*math.pi  # [rad/s/T]
slope = 0.36e-3  # s


# Trapezoidal gradient forms

def trapezoid_G(b, d, D, s=slope):
    """
    Calculate gradient's strength for a trapezoidal gradient form.

    Parameters
    ----------
    b : float
        bvalue in s/m^2 (= 10^6 s/mm^2)
    d : float
        gradient's length in s
    D : float
        time between the two gradients' center in s
        Note: This is not DW Delta Delta as provided to the EPI!
    s : float (Optional)
        ramptime of the gradient in s.

    Returns
    -------
    G : float
        Gradient's strength in T/m
    """
    g = gamma
    res = math.sqrt(
        b/(g**2*(d**2.*(D - d/3.) + s**3./30. - d*s**2./6))
        )
    return res


def trapezoid_b(G, d, D, s=slope):
    """
    Calculate the bvalue for a trapezoidal gradient form.

    Parameters
    ----------
    G : float
        Gradient's strength in T/m
    d : float
        gradient's length in s
    D : float
        time between the two gradients' center in s
        Note: This is not DW Delta Delta as provided to the EPI!
    s : float (Optional)
        ramptime of the gradient in s.

    Returns
    -------
    b : float
        bvalue in s/mm^2 (NO SI-unit!)
    """
    return gamma**2.*G**2*(d**2.*(D - d/3.) + s**3./30. - d*s**2./6.) * 1e-6


def trapezoid_D(G, b, d, s=slope):
    """
    Calculate the distance between the centers of two trapezoidal gradient
    forms.

    Parameters
    ----------
    G : float
        Gradient's strength in T/m
    b : float
        bvalue in s/m^2 (= 10^6 s/mm^2)
    d : float
        gradient's length in s
    s : float (Optional)
        ramptime of the gradient in s.

    Returns
    -------
    D : float
        time between the two gradients' center in s
        Note: This is not DW Delta Delta as provided to the EPI!
    """
    return (30.*b/(gamma*G)**2. + 10.*d**3. + 5.*d*s**2. - s**3.)/30./d**2.


def trapezoid_d(G, b, D, s):
    """
    Calculate the duration of a trapezoidal gradient form.

    There may be more than one mathematically meaningful solution, in this
    case the minimum meaningful solution is returned.

    Parameters
    ----------
    G : float
        Gradient's strength in T/m
    b : float
        bvalue in s/m^2 (= 10^6 s/mm^2)
    D : float
        time between the two gradients' center in s
        Note: This is not DW Delta Delta as provided to the EPI!
    s : float (Optional)
        ramptime of the gradient in s.

    Returns
    -------
    d : float
        gradient's length in s
    """
    gG = gamma*G
    a_ = -3.*D
    b_ = 0.5*s**2.
    c_ = -s**3./10. + 3*b/gG**2
    p = b_ - a_**2./3.
    q = 2*a_**3./27. - a_*b_/3. + c_
    D_ = (q/2.)**2. + (p/3.)**3.
    if D_ >= 0:
        z = ((-q/2. + math.sqrt(D_))**(1./3.)
             + (-q/2. - math.sqrt(D_))**(1./3.)
             - a_/3.)
    else:
        z_ = [
            (-math.sqrt(-4./3.*p)
                * math.cos(
                    1./3.*math.acos(-q/2.*math.sqrt(-27./p**3.))
                    + const/3.)
                - a_/3.)
            for const in [+math.pi, 0.0, -math.pi]]
        z = min([d for d in z_ if d > 0 and d < D])
    return z


# Rectangular gradient forms

def rectangular_G(b, d, D):
    return math.sqrt(b/(gamma**2.*d**2.*(D-d/3.)))


def rectangular_b(G, d, D):
    return gamma**2.*G**2.*d**2.*(D-d/3.)


def rectangular_D(b, G, d):
    return b/(gamma**2.*G**2.*d**2.) + d/3.


def rectangular_d(b, G, D, verbose=False):
    '''
    Calculate d in [s] using Delta.

    Returns three possible solutions if verbose == True and the minimal
    meaningful solution otherwise.
    '''
    g = gamma
    phi = math.acos(1.-3.*b/(2.*g**2*G**2.*D**3.))
    d1 = 2.*D*math.cos(phi/3.) + D
    d2 = D - 2.*D*math.cos((phi+math.pi)/3.)
    d3 = D - 2.*D*math.cos((phi-math.pi)/3.)
    if verbose:
        return d1, d2, d3
    else:
        res = [d_ for d_ in [d1, d2, d3] if d_ > 0 and d_ < D]
        return min(res)

#def create_list(Dd, b_0, d_0, length=50, dd = 1e-3):
#    data = np.empty((length,4))
#    data[:,0] = np.linspace(d_0, d_0+dd*(length-1), length)
#    G = calc_G_Dd(b_0,d_0,Dd)
#    data[0,0] = d_0
#    data[0,1] = b_0
#    for i in range(1,length):
#        data[i,1] = calc_b_Dd(G,data[i,0],Dd)
#    for i in range(length):
#        data[i,2] = round(data[i,1],-6)
#        data[i,3] = calc_d_Dd(data[i,2],G,Dd)
#    data[:,[0,3]] *= 1000.
#    data[:,1:3] /= 1000000.
#    return {'G':G,'Dd':Dd},data
