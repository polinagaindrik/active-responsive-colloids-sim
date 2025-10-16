# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate


def potential_V(r, sigma_ij):
    phi = np.zeros((np.size(r)))
    for i in range (np.size(r)):
        if r[i] < sigma_ij:
            phi[i] = (1 - r[i] / sigma_ij)**2.5
    return phi

# Calculate Radial distribution function (RDF) g(r)
def RDF_LDL0(x, T_sigma):
    sigma1 = np.linspace(0.01, 2, 100)
    sigma2 = np.linspace(0.01, 2, 100)
    intedr_func1 = np.zeros((np.size(x), np.size(sigma1), np.size(sigma2)))
    intedr_func2 = np.zeros((np.size(x), np.size(sigma1)))
    phi = Herz_pot(x, sigma1, sigma2)
    p01 = init_prob2(sigma1, T_sigma)
    p02 = init_prob2(sigma2, T_sigma)
    k_B = 1.0
    T_eff = 1.0
    beta = 1 / k_B * T_eff
    for j in range (np.size(sigma1)):
        for k in range (np.size(sigma2)):
            intedr_func1[:, j, k] = p01[j] * np.exp(- beta * phi[:, j, k])
            
    for k in range (np.size(sigma2)):              
        intedr_func2[:, k] = integrate.trapz(intedr_func1, sigma1, axis=1)[:, k] * p02[k]
    g = integrate.trapz(intedr_func2, sigma2, axis=1)
    return g  

# Herzian potential
def Herz_pot(x, sigma1, sigma2):
    eps = 500.0
    sigma12 = np.zeros((np.size(sigma1), np.size(sigma2)))
    phi = np.zeros((np.size(x), np.size(sigma1), np.size(sigma2)))
    for i in range (np.size(x)):
        for j in range (np.size(sigma1)):
            for k in range (np.size(sigma2)):
                sigma12[j, k] = 0.5 * (sigma1[j] + sigma2[k])
                if x[i] <= sigma12[j, k]:
                    phi[i, j, k] = eps * (1 - x[i] / sigma12[j, k])**2.5
    return phi


def init_prob2(sigma, T_sigma):
    Teff = T0
    var = 0.04 * T_sigma / Teff
    sigma0 = 1.0
    p = np.exp(- (sigma - sigma0)**2 / (2 * var)) / np.sqrt(2 * np.pi * var)
    return p


def B2coef(r, g):
    integrand = r**2 * (g - 1)
    B2 = - 2 * np.pi * integrate.trapz(integrand, r)
    return B2


if __name__ == "__main__":
    T0 = 1.0
    r = np.linspace(0, 10, 100)
    # Radial distribution function (RDF) 
    g = RDF_LDL0(r, 1.0)
    B2 = B2coef(r, g)