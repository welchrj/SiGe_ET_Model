# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 18:17:03 2020
This program calculates the material properties of Si_(1-x)Ge_x based off of
equations from: http://www.ioffe.ru/SVA/NSM/Semicond/SiGe/basic.html
In addition, the center beam temperature during L-PBF is calculated using
equations specified in, Wayne, K. E., et al. "Observation of keyhole-mode
laser melting in laser powder-bed fusion additive manufacturing" 2013, Journal
of Materials Processing Technology.

@author: welch
"""
import numpy
pi = numpy.pi
def thermal_cond(X):
    '''Computes thermal conductivity of SiGe in W/mK'''
    k = 0.046+0.084*X
    k = k/1e-2 #W/mK
    return k

def specific_heat(X):
    '''Computes specific heat of SiGe
    Parameters
    ----------
    X = ratio Si(1-X)Ge(X)'''
    
    Cp = 19.6+2.9*X #J/molk
    Si = 28.09*X
    Ge = 72.63*X
    Cp = (Cp*(1/(Si+Ge)))*1000 #J/kgK
    return Cp

def density(X):
    '''Computes the density of Si_(1-x)Ge_x
    Parameters
    ----------
    X = ratio of Si and Ge'''
    rho = 2.329 + 3.493*X - 0.499*X**2 #g/cm^3
    rho1 = rho/(1e-2)**3/1000
    return rho1

def melting_temperature(X):
    '''Computes the liquidus temperature (K) of Si_(1-x)Ge_x'''
    T_l = 1412-80*X-295*X**2+273.15
    
    return T_l

def laser_intensity(P, s):
    '''Computes laser intensity from King, et al.
    Parameters
    ----------
    P: Laser power (W)
    s: spot size (diameter) (m)
    Returns
    -------
    I: Intensity (W/m^2'''
    
    I = P/(2*pi*s**2)
    
    return I
def center_temperature(A, P, s, k, D, u):
    '''Calculates the temperature at the center of a stationary beam assumes
    all of the laser power is absorbed in the material
    
    Parameters
    ----------
    A: absorptivity
    P: laser power (W)
    s: spot size
    k: thermal conductivity
    D: thermal diffusiivity
    u: scan speed
    
    Returns
    -------
    T = temperature at center of beam'''
    
    t = s/u
    I = laser_intensity(P, s)
    T = ((numpy.sqrt(2)*A*I*s)/(k*numpy.sqrt(pi)))*numpy.arctan((numpy.sqrt((2*D*t)/s)))
    
    return T

def keyhole_threshold(D, s, k, Tb, A, P):
    '''Calculates the keyhole threshold scan velocity
    
    Parameters
    ----------
    D: thermal diffusivity
    s: spot size
    k: thermal conductivity
    Tb:boiling temperature
    A: absorptivity
    P: laer power
    
    Returns
    -------
    Ut: keyhole threshold velocity'''
    
    I =laser_intensity(P,s)
    Ut = ((4*D)/s)*((numpy.sqrt(pi)*k*Tb)/(A*I*s))**(-2)
    
    return Ut
    
def main():
    X = 0.5
    cp = specific_heat(X)
    print('Specific Heat [J/kgK]: {0}'.format(cp))
    rho1 = density(X)
    print('Density [kg/m^3]: {0}\n[g/cm^3]:{1}'.format(rho1, rho))
    k = thermal_cond(X)
    print('Thermal conductivity [W/mK]: {0}'.format(k))
    Tb = 1400
    #absorptivity
    A = 0.67
    #beam spot size (diameter)
    s = 50e-6
    #thermal diffusivity
    d = k/(rho1*cp) #m^-2 s^-1
    P = numpy.array((25, 25, 25, 50, 50, 50, 75, 75, 75, 100, 100, 100))
    V = numpy.array((.1, .3, .5, .1, .3, .5, .1, .3, .5, .1, .3, .5))
    T = center_temperature(A, P, s, k, d, V)
    print(T)
    V = keyhole_threshold(d, s, k, Tb, A, P)
    print(V)
    
if __name__=='__main__':
    main()

