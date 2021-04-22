#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 08:59:17 2021

@author: jorgenorena
"""

from Vectors import Vector, VectSum, VectScalMul, Dot
from Vectors import Scalar, ScalPow, ScalMul, ScalSum
from Vectors import fract, sympy, translate_sympy, vectors, scalars
from itertools import permutations
from scipy.special import factorial
import sympy as sp

def alpha(k1, k2):
    return Dot((k1 + k2), k1)/Dot(k1,k1)

def beta(k1, k2):
    return (Dot(k1 + k2, k1 + k2)*Dot(k1, k2))/Dot(k1,k1)/Dot(k2,k2)/2

def G(n, ks):
    if n == 1:
        return 1
    else:
        return sum(Gs(m, ks[:m])*\
            (3*alpha(sum(ks[:m]),sum(ks[m:]))*Fs(n-m, ks[m:]) +\
            2*n*beta(sum(ks[:m]),sum(ks[m:]))*Gs(n-m, ks[m:]))/(2*n+3)/(n-1)
            for m in range(1, n))
            
def F(n, ks):
    if n == 1:
        return 1
    else:
        return sum(Gs(m, ks[:m])*\
            ((2*n+1)*alpha(sum(ks[:m]),sum(ks[m:]))*Fs(n-m, ks[m:]) +\
            2*beta(sum(ks[:m]),sum(ks[m:]))*Gs(n-m, ks[m:]))/(2*n+3)/(n-1)
            for m in range(1, n))
            
def Fs(n, ks):
    return ScalPow(int(factorial(n)),-1)*\
        sum(map(lambda x: F(n, x), permutations(ks)))
    
def Gs(n, ks):
    return ScalPow(int(factorial(n)),-1)\
        *sum(map(lambda x: G(n, x), permutations(ks)))

def F2(k1, k2):
    return Fs(2, [k1, k2])

if __name__=='__main__':
    k1, k2 = vectors('k_1', 'k_2')
    ks = [k1, k2]
    m = 1
    n = 2
    
    f = Fs(2, ks)
    sf = sympy(f).simplify()
    print(sp.latex(sf))