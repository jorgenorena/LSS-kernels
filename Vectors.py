#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 07:23:49 2021

@author: jorgenorena
based on code writen by Joaquin Rohland
"""

import sympy as sp

symbol_dict = {}

def vectors(names):
    return tuple(Vector(name) for name in names)

def scalars(names):
    return tuple(Scalar(name) for name in names)

def srepr(a):
    try:
        return a.srepr()
    except AttributeError:
        return repr(a)

def _gcd(a, b):
    while b:
        a, b = b, a%b
    return a

def prod(args):
    res = 1
    for a in args:
        res *= a
    return res

def is_number(n):
    return isinstance(n, int) or isinstance(n, float)

def is_not_number(n):
    return not (isinstance(n, int) or isinstance(n, float))

def is_vector(expr):
    try:
        return expr.vector
    except AttributeError:
        return False

def is_scalar(expr):
    try:
        return expr.scalar
    except AttributeError:
        if isinstance(expr, int) or isinstance(expr, float):
            return True
        else:
            return False

class VExpr():
    
    def __init__(self):
        self._mhash = None
    
    
    def __hash__(self):
        h = self._mhash 
        if h is None:
            h = hash((type (self).__name__,) + \
                     tuple (self.args))
            self._mhash = h
        return h
        
    
    def srepr(self):
        return type(self).__name__ + '(' +\
            ', '.join(srepr(a) for a in self.getargs()) + ')'
            
    def __repr__(self):
        return self.srepr()
    
    def __add__(self, other):
        
        if other == 0:
            return self
        
        if is_vector(self):
            if is_vector(other):
                return VectSum(self, other)
            else:
                raise TypeError("unsupported operand type(s) for +: " +\
                                type(self).__name__ + ' and ' +\
                                    type(other).__name__)
            
        elif is_scalar(self):
            if is_scalar(other):
                return ScalSum(self, other)
            else:
                raise TypeError("unsupported operand type(s) for +: " +\
                                type(self).__name__ + ' and ' +\
                                    type(other).__name__)

        
    def __radd__(self, other):
        return self + other
    
    
    def __mul__(self, other):
        if is_vector(self):
            if is_scalar(other):
                return VectScalMul(other, self)
            elif is_vector(other):
                return Dot(self, other)
            else:
                raise TypeError("unsupported operand type(s) for *: " +\
                                type(self).__name__ + ' and ' +\
                                    type(other).__name__)
                
        elif is_scalar(self):
            if is_scalar(other):
                return ScalMul(self, other)
            elif is_vector(other):
                return VectScalMul(self, other)
            else:
                raise TypeError("unsupported operand type(s) for *: " +\
                                type(self).__name__ + ' and ' +\
                                    type(other).__name__)
        
        
    def __rmul__(self, other):
        return self * other
    
    def __eq__(self, other):
        return hash(self) == hash(other)
    
    def __truediv__(self, other):
        if is_vector(self):
            if is_scalar(other):
                return VectScalMul(ScalPow(other, -1), self)
            else:
                return TypeError("3 unsupported operand type(s) for /: " +\
                                type(self).__name__ + ' and ' +\
                                    type(other).__name__)
        
        elif is_scalar(self):
            if is_scalar(other):
                return ScalMul(ScalPow(other, -1), self)
            else:
                raise TypeError("2 unsupported operand type(s) for /: " +\
                                type(self).__name__ + ' and ' +\
                                    type(other).__name__)
        
        else:
            raise TypeError("1 unsupported operand type(s) for /: " +\
                                type(self).__name__ + ' and ' +\
                                    type(other).__name__)
        
    def __rtruediv__(self, other):
        if is_scalar(self):
            if is_vector(other):
                return VectScalMul(ScalPow(self, -1), other)
            elif is_scalar(other):
                return ScalMul(ScalPow(self, -1), other)
            else:
                raise TypeError("4 unsupported operand type(s) for /: " +\
                                type(self).__name__ + ' and ' +\
                                    type(other).__name__)
        else:
            raise TypeError("5 unsupported operand type(s) for /: " +\
                                type(self).__name__ + ' and ' +\
                                    type(other).__name__)
                    
                    
    
    def getargs(self):
        return self.args
 
class Associative():   
    
    def make_associative(self):    
        new_args = []
        for a in self.args:
            if type(a) == type(self):
                new_args.extend(a.args)
            else:
                new_args.append(a)
        self.args = tuple(new_args)
        
    
class Commutative():  

    def make_commutative(self):
        constlist = list(filter(is_number, self.args))
        arglist = sorted(list(filter(is_not_number, self.args)), key=hash)
        if len(constlist) > 0:
            arglist = [self._number_version(*constlist)] + arglist
        self.args = tuple(arglist)


class Identity():
    
    def ignore_identity(self):
        self.args = tuple(filter(self._identity_.__ne__, self.args))
    

class NullElement():
    
    def is_null(self):
        return self._null_ in self.args
    
    
class Cummulative():
    
    def simplify(self, repeated, operate, separate):
        previous = None
        c = None
        terms = []
        def key(term):
            ci, t  = separate(term)
            return hash(t)
        args = sorted(self.args, key=key)
        for term in args:
            ci, current = separate(term)
            if current != previous:
                if previous != None:
                    terms.append(repeated(previous, c))
                c = ci
                previous = current
            else:
                c = operate(c, ci)
        terms.append(repeated(previous, c))
        self.args = tuple(terms)
    
    
class Vector(VExpr):
    
    def __init__(self, name, components = None):
        if not isinstance(name, str):
            raise ValueError('Name of vector must be a string.')
        self.name = name
        self.components = components
        self.vector = True
        self.args = (name,)
        self._mhash = None
    
    def getargs(self):
        return self.args
    
    def latex(self):
        return '\\vec{' + self.name + '}'
    
    def __repr__(self):
        return 'v('+self.name+')'
    
    def srepr(self):
        return repr(self)
    
    def set_components(self, components):
        self.components = components
        
    def val(self):
        if self.components != None:
            return self.components
        else:
            return self
    
class VectSum(VExpr, Associative, Commutative, Identity, Cummulative):
    
    def __new__(cls, *args):
        if not all(map(is_vector, args)):
            raise TypeError('VectSum should only involve vector objects.')
        instance = super(VectSum, cls).__new__(cls)
        instance._identity_ = Vector('0')
        instance.args = args
        instance.make_associative()
        instance.make_commutative()
        instance.ignore_identity()
        instance.simplify(lambda a, b: VectScalMul(b, a), ScalSum, 
                          instance._separate_scal)
        if len(instance.args) == 1:
            return instance.args[0]
        else:
            return instance
    
    def __init__(self, *args):
        self.vector = True
        self._mhash = None 

    def _separate_scal(self, term):
        if isinstance(term, VectScalMul) and is_scalar(term.args[0]):
            return term.args[0], term.args[1]
        else:
            return 1, term

    def __repr__(self):
        return '(' + ' + '.join(repr(a) for a in self.args) + ')'
        
    def latex(self):
        l = [a.latex() for a in self.args]
        return '(' + ' + '.join(l) + ')'
    
    def val(self):
        return sum(a.val() for a in self.args)
    

class VectScalMul(VExpr, Identity):
    
    def __new__(cls, *args):
        if len(args) != 2:
            raise TypeError('VectScalMul takes 1 argument, ' + len(args) +\
                            'were given')
        elif not is_scalar(args[0]) or not is_vector(args[1]):
            raise TypeError('VectScalMul takes a scalar and a vector.')
        
        instance = super(VectScalMul, cls).__new__(cls)
        instance.args = args
        instance._identity_ = 1
        instance.ignore_identity()
        if len(instance.args) == 1:
            return instance.args[0]
        elif instance.args[0] == 0 or instance.args[1] == Vector('0'):
            return Vector('0')
        elif isinstance(instance.args[1], VectSum):
            return VectSum(*[instance.args[0]*v 
                             for v in instance.args[1].args])
        else:
            return instance
        
    def __init__(self, *args):
        self.vector = True
        self._mhash = None
        
    def __repr__(self):
        return repr(self.args[0]) + '*' + repr(self.args[1])
        
    def latex(self):
        return self.args[0].latex() + ' ' + self.args[1].latex()
    
    def val(self):
        return self.args[0].val()*self.args[1].val()


class Dot(VExpr, Commutative, NullElement):
    
    def __new__(cls, *args):
        if len(args) != 2:
            raise TypeError('Dot takes 1 argument, ' + len(args) +\
                            'were given')
        if not all(map(is_vector, args)):
            raise TypeError('Dot should only involve vector objects.')
        instance = super(Dot, cls).__new__(cls)
        instance._null_ = Vector('0')
        instance.args = args
        instance.make_commutative()
        if instance.is_null():
            return 0
        # if there are sums or products with scalars, expand
        if isinstance(instance.args[0], VectSum) \
            or isinstance(instance.args[1], VectSum):
            return instance.expand_sum()
        elif isinstance(instance.args[0], VectScalMul) \
            or isinstance(instance.args[1], VectScalMul):
            return instance.expand_mul()
        else:
            return instance
     
    def expand_sum(self):
        if isinstance(self.args[0], VectSum):
            terms0 = self.args[0].args
        else:
            terms0 = [self.args[0]]
        
        if isinstance(self.args[1], VectSum):
            terms1 = self.args[1].args
        else:
            terms1 = [self.args[1]]
            
        sumargs = [Dot(t1, t2) for t1 in terms0 for t2 in terms1]
        return ScalSum(*sumargs)
    
    
    def expand_mul(self):
        if isinstance(self.args[0], VectScalMul):
            s1, t1 = self.args[0].args
        else:
            s1 = 1
            t1 = self.args[0]
        
        if isinstance(self.args[1], VectScalMul):
            s2, t2 = self.args[1].args
        else:
            s2 = 1
            t2 = self.args[1]
        
        return ScalMul((s1*s2),Dot(t1, t2))
    
        
    def __init__(self, *args):
        self.scalar = True
        self._mhash = None
        self.symbol = None
        
    def __repr__(self):
        return '(' + repr(self.args[0]) + '.' + repr(self.args[1]) + ')'
    
    def latex(self):
        return self.args[0].latex() + ' \cdot ' + self.args[1].latex()
    
    def val(self):
        return sum(i[0]*i[1] 
                   for i in zip(self.args[0].val(), self.args[1].val()))
    
    def sympy(self):
        if self.symbol == None:
            self.make_symbol()
        return self.symbol
    
    def make_symbol(self):
        self.symbol = sp.Symbol(self.latex())
        symbol_dict[self.symbol] = self
        
    def clean_symbol(self):
        self.symbol = None
        del symbol_dict[self.symbol]
    

    
class Scalar(VExpr):
    
    def __init__(self, name, value = None):
        if not isinstance(name, str):
            raise ValueError('Name of scalar must be a string.')
        self.name = name
        self.value = value
        self.scalar = True
        self.symbol = None
        self.args = (name,)
        self._mhash = None
        
    def getargs(self):
        return (self.name,)
    
    def __repr__(self):
        return self.name
        
    def latex(self):
       return self.name
   
    def srepr(self):
        return repr(self)
   
    def set_value(self, value):
        self.value = value
   
    def val(self):
        if self.value != None:
            return self.value
        else:
            return self
    
    def sympy(self):
        if self.symbol == None:
            self.make_symbol()
        return self.symbol
    
    def make_symbol(self):
        self.symbol = sp.symbols(self.name)
        symbol_dict[self.symbol] = self
        
    def clean_symbol(self):
        self.symbol = None
        del symbol_dict[self.symbol]


class ScalPow(VExpr):
    
    def __new__(cls, *args):
        if not all(map(is_scalar, args)):
            raise TypeError('ScalPow should only involve Scalar objects.')
        instance = super(ScalPow, cls).__new__(cls)
        instance.args = args
        if len(instance.args) > 2:
            raise TypeError('ScalPow takes 2 arguments but ' + \
                            len(instance.args)  + ' were given.')
        elif len(instance.args) == 1:
            return instance.args[0]
        elif instance.args[0] == 1:
            return 1
        elif instance.args[1] == 0:
            return 1
        elif instance.args[1] == 1:
            return instance.args[0]
        elif isinstance(instance.args[0], float) and \
            is_number(instance.args[1]):
            return instance.args[0]**instance.args[1]
        elif is_number(instance.args[0]) and \
            isinstance(instance.args[1], float):
            return instance.args[0]**instance.args[1]
        else:
            return instance
        
    def __init__(self, *args):
        self._mhash = None
        self.scalar = True
        self.base = self.args[0]
        self.exp = self.args[1]
        
    def __repr__(self):
        if is_number(self.base):
            base_string = str(self.base)
        else:
            base_string = repr(self.base)
        
        if is_number(self.exp) and self.exp < 0:
            if self.exp == -1:
                return '(1/' + base_string + ')'
            exp_string = str(-self.exp)
            return '(1/(' + base_string + '^' + exp_string + '))'
        else:
            exp_string = self.exp if is_number(self.exp) else repr(self.exp)
            return '(' + base_string + ')^' + exp_string
        
    def latex(self):
        if is_number(self.base):
            base_string = str(self.base)
        else:
            base_string = self.base.latex()
        
        if is_number(self.exp) and self.exp < 0:
            exp_string = str(-self.exp)
            return '\\frac{1}{' + base_string + '^' + exp_string + '}'
        else:
            exp_string = self.exp if is_number(self.exp) else self.exp.latex()
            return '(' + base_string + ')^' + exp_string
    
    def val(self):
        return self.base.val()**self.exp.val()
    
    def sympy(self):
        return sp.Pow(sympy(self.base), sympy(self.exp))


class ScalMul(VExpr, Associative, Commutative, Identity, Cummulative, 
              NullElement):
    
    def __new__(cls, *args):
        if not all(map(is_scalar, args)):
            raise TypeError('ScalMul should only involve Scalar objects.')
        instance = super(ScalMul, cls).__new__(cls)
        instance._identity_ = 1
        instance._null_ = 0
        instance.args = args
        instance.make_associative()
        instance.make_commutative()
        if instance.is_null():
            return 0
        instance.simplify(ScalPow, lambda a, b: a + b,
                          instance._separate_exp)
        instance.ignore_identity()
        if len(instance.args) == 1:
            return instance.args[0]
        elif all([is_number(a) for a in instance.args]):
            return prod(instance.args)
        else:
            return instance
        
    def __init__(self, *args):
        self.scalar = True
        self._mhash = None
        
    def __repr__(self):
        s = [self._separate_exp(a) for a in self.args]
        numer = ''
        denom = ''
        for p, b in s:
            b_str = str(b) if is_number(b) else repr(b)
            if is_number(p) and p < 0:
                if p == -1:
                    denom += b_str
                else:
                    p_str = str(-p)
                    denom += ' ' + b_str + '^' + p_str
            else:
                if p == 1:
                    numer += ' ' + b_str
                else:
                    p_str = str(p) if is_number(p) else repr(p)
                    numer += ' ' + b_str + '^' + p_str
        if len(numer) == 0:
            numer = str(1)
        else:
            numer = numer[1:]
        if len(denom) == 0:
            return numer
        else:
            return '(' + numer + '/(' + denom + '))'
    
    def latex(self):
        s = [self._separate_exp(a) for a in self.args]
        numer = ''
        denom = ''
        for p, b in s:
            b_str = str(b) if is_number(b) else b.latex()
            if is_number(p) and p < 0:
                if p == -1:
                    denom += b_str
                else:
                    p_str = str(-p)
                    denom += ' ' + b_str + '^' + p_str
            else:
                if p == 1:
                    numer += ' ' + b_str
                else:
                    p_str = str(p) if is_number(p) else p.latex()
                    numer += ' ' + b_str + '^' + p_str
        if len(numer) == 0:
            numer = str(1)
        else:
            numer = numer[1:]
        if len(denom) == 0:
            return numer
        else:
            return '\\frac{' + numer + '}{' + denom + '}'
    
    def val(self):
        return prod([a.val() for a in self.args])
    
    def sympy(self):
        return sp.Mul(*[sympy(a) for a in self.args])
    
    def _number_version(self, *args):
        return prod(args)
    
    def _separate_exp(self, term):
        if isinstance(term, ScalPow) and is_number(term.args[1]):
            return term.args[1], term.args[0]
        else:
            return 1, term

class ScalSum(VExpr, Associative, Commutative, Identity, Cummulative):
    
    def __new__(cls, *args):
        if not all(map(is_scalar, args)):
            raise TypeError('ScalSum should only involve Scalar objects.')
        instance = super(ScalSum, cls).__new__(cls)
        instance._identity_ = 0
        instance.args = args
        instance.make_associative()
        instance.make_commutative()
        instance.ignore_identity()
        instance.simplify(ScalMul, instance._number_version, 
                          instance._separate_num)
        if len(instance.args) == 1:
            return instance.args[0]
        if all([is_number(a) for a in instance.args]):
            return sum(args)
        else:
            return instance
    
    def __init__(self, *args):
        self.scalar = True
        self._mhash = None 

    def _separate_num(self, term):
        if isinstance(term, ScalMul) and is_number(term.args[0]):
            return term.args[0], ScalMul(*term.args[1:])
        else:
            return 1, term
        
    def __repr__(self):
        l = [(str(a) if is_number(a) else repr(a)) for a in self.args]
        return '(' + ' + '.join(l) + ')'
        
    def latex(self):
        l = [(str(a) if is_number(a) else a.latex()) for a in self.args]
        return '(' + ' + '.join(l) + ')'
    
    def val(self):
        return sum(a.val() for a in self.args)
    
    def sympy(self):
        return sp.Add(*[sympy(a) for a in self.args])
    
    def _number_version(self, *args):
        return sum(args)

def fract(numer, denom):
    if (isinstance(numer, float) and is_number(denom)) or\
        (isinstance(denom, float) and is_number(numer)):
        return numer/denom
    elif not is_scalar(numer) or not is_scalar(denom):
        raise TypeError('Can only divide scalars (or numbers).')
    else:
        return ScalMul(numer, ScalPow(denom, -1))

def sympy(expr):
    try:
        return expr.sympy()
    except AttributeError:
        return expr

def is_sympy_number(tree):
    return isinstance(tree, sp.Float) or isinstance(tree, sp.Integer)

def translate_sympy(tree):
    if isinstance(tree, sp.Add):
        return(ScalSum(*[translate_sympy(a) for a in tree.args]))
    elif isinstance(tree, sp.Mul):
        return(ScalMul(*[translate_sympy(a) for a in tree.args]))
    elif isinstance(tree, sp.Pow):
        return ScalPow(*[translate_sympy(a) for a in tree.args])
    elif isinstance(tree, sp.Rational):
        return fract(translate_sympy(tree.p), translate_sympy(tree.q))
    elif isinstance(tree, sp.Symbol) or is_sympy_number(tree):
        return symbol_dict[tree]
    elif is_number(tree):
        return tree
    else:
        print('Could not translate:', tree)
    
if __name__ == '__main__':
    
    n=3
    k1, k2 = vectors(['k_1', 'k_2'])
    print(sympy(Dot(k1, k1) + Dot(k2, k1)))