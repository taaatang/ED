#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 11:30:41 2021

@author: tatang
"""
#%% import
import numpy as np
from scipy.special import sph_harm
from scipy.special import eval_genlaguerre
from scipy.special import factorial
import scipy.integrate as integrate
from numba import njit

import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab
#%% wavefunction
# @jit
def WF(n, l, m, x, y, z):
    a0 = 0.8
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z/r)
    phi = np.arctan2(y,x)
    r = 2*r/n/a0
    return np.sqrt((2/n/a0)**3*factorial(n-l-1)/2/n/factorial(n+l))*np.exp(-r/2)*r**l*eval_genlaguerre(n-l-1,2*l+1,r)*sph_harm(m,l,phi,theta)
# @jit
def WF_2px(x, y, z):
#     return ((WF(2, 1, 1, x, y, z) + WF(2, 1, -1, x, y, z))/np.sqrt(2)).real
    return -np.real(WF(2, 1, 1, x, y, z))
# @jit
def WF_2py(x, y, z):
#     return ((WF(2, 1, 1, x, y, z) - WF(2, 1, -1, x, y, z))/np.sqrt(2)/1j).real
    return -np.imag(WF(2, 1, 1, x, y, z))
# @jit
def WF_2pz(x, y, z):
    return np.real(WF(2, 1, 0, x, y, z))
# @jit
def WF_3d(x, y, z):
#     return ((WF(3, 2, 1, x, y, z) - WF(3, 2, -1, x, y, z))/np.sqrt(2)/1j).real
    return np.real(WF(3, 2, 2, x, y, z))
# @jit
def WForb(x, y, z, d, orb):
    if (orb == 0):
        return WF_3d(x, y, z)
    elif (orb == 1):
        return WF_2px(x - d, y, z)
    elif (orb == 2):
        return WF_2py(x, y - d, z)
    elif (orb == 3):
        return WF_2pz(x, y, z - d)
    elif (orb == -1):
        return WF_2px(x + d, y, z)
    elif (orb == -2):
        return WF_2py(x, y + d, z)
    elif (orb == -3):
        return WF_2pz(x, y, z + d)

def drawOrbWave(coords, orbs, L = 15, d = 40, contNum = 4, cmap = 'seismic'):
    x, y, z = np.mgrid[-L:L:300j, -L:L:300j, -L:L:300j]
    coords *= d
    dhf = d * 0.5
    for orb in orbs:
        data = WForb(x , y , z, dhf, orb) 
        for coord in coords:
            mlab.contour3d(x + coord[0], y + coord[1], z + coord[2], data, colormap = cmap, contours = contNum, transparent = True)
            mlab.show()

def plotSphere(ax, center, r, color, alpha):
    # Make data
    u = np.linspace(0, 2 * np.pi, 25)
    v = np.linspace(0, np.pi, 25)
    x = r * np.outer(np.cos(u), np.sin(v)) + center[0]
    y = r * np.outer(np.sin(u), np.sin(v)) + center[1]
    z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]
    # Plot the surface
    ax.plot_surface(x, y, z, color = color, alpha = alpha)
#%% plot
def add_subplot_axes(ax, rect, axisbg = 'w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]
    subax = fig.add_axes([x,y,width,height])
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax


def text_y(ax, height):
    lim = ax.get_ylim()
    return lim[0] + height * (lim[1] - lim[0])
def text_x(ax, width):
    lim = ax.get_xlim()
    return lim[0] + width * (lim[1] - lim[0])
#%% spectra
@njit
def Continue_frac(z, an, bn, n, s = 0.0):
    """continuous franction used for calculating dynamic structure with ED

    Args:
        z (complex64): frequency
        an (n-vector): alpha coefficients from Lanczos
        bn ((n-1)-vector): beta coefficients from Lanczos
        n (integer): length of an
        s (complex64, optional): previous Continue_frac value. Defaults to 0.

    Returns:
        complex64: dynamic structure value at z
    """
    if n == 0:
        return 0
    if n > 1:
        s = bn[n-2]**2/(z-an[n-1]-s)
        return Continue_frac(z,an,bn,n-1,s)
    if n == 1:
        return 1/(z-an[n-1]-s)

@njit
def Conductivity(w, w0, an, bn, norm, epsilon = 0.02): 
    '''
    regular part of the conductivity
    '''
    n = len(an)
    z = w+w0+1j*epsilon
    sig = - (1/2/(1j*w)*norm*Continue_frac(z, an, bn, n))
    return sig

@njit
def Spectra(w, w0, an, bn, norm, epsilon = 0.02):
    '''
    spectra
    '''
    n = len(an)
    z = w+w0+1j*epsilon
    sig = - 1/np.pi*np.imag(norm*Continue_frac(z,an,bn,n))
    return sig

@njit
def Spectra_inv(w, w0, an, bn, norm, epsilon = 0.02):
    '''
    spectra
    '''
    n = len(an)
    z = w0-w-1j*epsilon
    sig = 1/np.pi*np.imag(norm*Continue_frac(z, an, bn, n))
    return sig

@njit
def findmax(sig, wlist):
    for idx,w in enumerate(wlist):
        if(idx>0 and idx<(len(wlist)-1)):
            if(sig[idx]>sig[idx-1] and sig[idx]>sig[idx+1]):
                return w

@njit
def findmax_rev(sig, wlist):
    for idx,w in enumerate(wlist[-1::-1]):
        if(idx>0 and idx<(len(wlist)-1)):
            if(sig[idx]>sig[idx-1] and sig[idx]>sig[idx+1]):
                return w

def read_cond(condPath, wlist):
    nums = 400
    w0 = np.fromfile(condPath + "/w0", dtype = np.complex128)
    vecNorm = np.fromfile(condPath + "/vecNorm",dtype=np.float64)
    alpha = np.fromfile(condPath + "/alpha",dtype=np.float64)
    beta = np.fromfile(condPath + "/beta",dtype=np.float64)
    w0 = np.real(w0[0])
    norm = abs(vecNorm[0])
    sig = []
    for w in wlist:
        sig.append(Conductivity(w,w0,alpha[:nums],beta,norm,epsilon=0.02)+Conductivity(-w,w0,alpha[:nums],beta,norm,epsilon=0.02))
    return w0, np.real(np.array(sig))

def read_spec(specPath, wlist, state = 0, delta = 0.02, nums = -1):
    """read spectra from files

    Args:
        specPath ([type]): [description]
        wlist ([type]): [description]
        state (int, optional): [description]. Defaults to 0.
        delta (float, optional): [description]. Defaults to 0.02.
        nums (int, optional): [description]. Defaults to -1.

    Returns:
        [type]: [description]
    """
    label = ""
    if state >= 0:
        label = "_" + str(state)
    w0 = np.fromfile(specPath + "/w0" + label, dtype = np.complex128)
    vecNorm = np.fromfile(specPath + "/vecNorm" + label,dtype=np.float64)
    alpha = np.fromfile(specPath + "/alpha" + label,dtype=np.float64)
    beta = np.fromfile(specPath + "/beta" + label,dtype=np.float64)
    w0 = np.real(w0[0])
    norm = abs(vecNorm[0])
    sig = []
    for w in wlist:
        sig.append(Spectra(w,w0,alpha[:nums],beta,norm,epsilon=delta))
    return np.array(sig)

def read_spec_inv(specPath, wlist, state = 0, delta = 0.02, nums = -1):
    label = ""
    if state >= 0:
        label = "_" + str(state)
    w0 = np.fromfile(specPath + "/w0" + label, dtype = np.complex128)
    vecNorm = np.fromfile(specPath + "/vecNorm" + label,dtype=np.float64)
    alpha = np.fromfile(specPath + "/alpha" + label,dtype=np.float64)
    beta = np.fromfile(specPath + "/beta" + label,dtype=np.float64)
    w0 = np.real(w0[0])
    norm = abs(vecNorm[0])
    sig = []
    for w in wlist:
        sig.append(Spectra_inv(w,w0,alpha[:nums],beta,norm,epsilon=delta))
    return np.array(sig)
#%% time evolution
w0 = 1.0
t0 = 1.05457/1.60218 # fs, unit time
L0 = 300*t0
E_field0 = 16.0218/3/1.055*1e6 #V/m
lam =1241.5267/w0# pulse central wavelength: nm
width = 10 # pulse width: fs
t1 = 50/t0 #fs/t0
sigma0 = width*np.sqrt(2)/t0
a0 = 0.2429/L0 # Lattice constant: nm/L0
epsilon = 8.85*1e-12
c = 3*1e8
def E(t,psi=np.pi/2,tc=t1,sigma=sigma0,w=w0,E0=E_field0):
    return (abs(t-tc)<5*sigma).astype(int)*E0*np.exp(-((t-tc)/sigma)**2)*np.sin(w*(t-tc)+psi) 
def A(t,tc=t1,sigma=sigma0):
    T1 = tc+5*sigma
    if t<T1:
        return -integrate.quad(E,0,t,limit=250)[0]
        return -integrate.quad(E,0,T1,limit=250)[0]
def E2(t):
    return E(t)**2
def Energy(t1,t2,t0=t0,E0=E_field0,c=3e8,epsilon=8.85e-12):
    # retunerrn fluence mJ/cm^2
    return c*epsilon*t0*1e-15*E0**2*integrate.quad(E2,t1,t2,limit=500)[0]/10
# %%
def getRamanLabel(channel):
    name = channel[0]
    sub = "" if len(channel) == 1 else "_{" + channel[1] + "}"
    sup = "" if len(channel) == 2 else "^{(" + channel[2:] + ")}"
    label = "$" + name + sub  + sup + "$"
    return label