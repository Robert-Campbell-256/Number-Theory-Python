################################################################################
# Elliptic Curve
# Elliptic curves in reduced Weierstrass form over prime order fields
# Author: Robert Campbell, <campbell@math.umbc.edu>
# Date: 18 May, 2014
# Version 0.2
# License: Simplified BSD (see details at bottom)
######################################################################################
"""Elliptic Curve
	EllipticCurve(p,[a,b]) is the elliptic curve in affine restricted Weierstrass form
	y^2 = x^3 +ax + b (mod p)
	Usage: 
		>>> from ellipticcurve import *
		>>> ec29 = EllipticCurve(29,[4,20])
		>>> print ec29
			y^2 = x^3 + 4x + 20 (mod 29)
		>>> pt = EllipticCurveElt(ec29,[2,6])
		>>> print 7*pt
			(3, 28)
		>>> print (9*pt) - pt
			(15, 27)
		>>> ec29
			EllipticCurve(29, [4,20])
		>>> pt
			EllipticCurveElt(EllipticCurve(29, [4,20]), [2,6])
"""   
   
Version = 'ELLIPTICCURVE.PY, version 0.1, 9 May, 2014, by Robert Campbell, <campbell@math.umbc.edu>'
import numbthy # For xgcd (for modinv) and sqrtmod

# Assumptions: Affine (later Projective?) Reduced Weierstrass form
#	over prime field.  Thus identity is point at infinity,
#	and -P is gotten from P by negating the y value.
#	y^2 = x^3 + ax + b (mod p)
#	Refs:
#		[HMV04] Guide to Elliptic Curve Cryptography by Hankerson, Menezes & Vanstone, 2004
#			(see in particular, Sects 3.1.1 & 3.1.2)
#		[Wash03] Elliptic Curves, L. Washington, 2003 (see Sect 2.2)
#		[Many, many other decent references]
#	Addition Rules:
#		i) 0 + P = P
#		ii) P + (-P) = 0 [i.e. x1==x2 but y1==-y2]
#		iii) P + P [i.e. x1==x2 and y1==y2]
#			lambda = (3x^2+a)/2y   # "tangent slope"
#			x3 = lambda^2 - 2x
#			y3 = lambda*(x-x3) - y
#		iv) P1 + P2 [i.e. x1!=x2]
#			lambda = (y1-y2)/(x1-x2) # "slope"
#			x3 = lambda^2 - x1 - x2
#			y3 = lambda*(x1-x3) - y1
#	Zero point (aka point at infinity) represented as ["Infinity","Infinity"]

class EllipticCurve(object):
	"""Elliptic Curve
		EllipticCurve(p,[a,b]) is the elliptic curve in affine restricted Weierstrass form
		y^2 = x^3 +ax + b (mod p)
		Usage: 
			>>> from ellipticcurve import *
			>>> ec29 = EllipticCurve(29,[4,20]) # Define the elliptic curve over GF(29)
			>>> print ec29
				y^2 = x^3 + 4x + 20 (mod 29)
			>>> pt = EllipticCurveElt(ec29,[2,6])  # Define a point on the curve
			>>> print pt+pt, pt-pt  # Add the point to itself, then subtract (yielding zero pt on curve)
				(1, 5) (Infinity, Infinity)
			>>> print 5*pt  # Compute 5*pt = pt+pt+pt+pt+pt
				(0, 7)
			>>> repr(ec29)  # Format usable as input
				EllipticCurve(29, [4,20])
			>>> repr(pt)    # Format usable as input
				EllipticCurveElt(EllipticCurve(29, [4,20]), [2,6])
	"""   
	def __init__(self,prime,coeffs,verbose=False,fmtspec="d"):
		self.prime = prime
		if(not(numbthy.isprime(self.prime))): raise ValueError("***** Error *****: Characteristic of base field {0} must be prime".format(self.prime))
		self.a = coeffs[0]
		self.b = coeffs[1]
		self.discriminant = -16*(4*(self.a**3)+27*(self.b**2)) % self.prime
		if(self.discriminant == 0): raise ValueError("***** Error *****: Not an elliptic curve - Zero discriminant (-16*(4*({0}^3)+27*({1}^2)))".format(self.a,self.b))
		self.verbose = verbose
		self.fmtspec = fmtspec
		if(self.verbose): print "Elliptic Curve {0} of discriminant {1}".format(self,self.discriminant)
	def __format__(self,fmtspec):  # Over-ride format conversion
		return "y^2 = x^3 + {1:{0}}x + {2:{0}} (mod {3:{0}})".format(self.fmtspec,self.a,self.b,self.prime)
	def __str__(self):   # Over-ride string conversion used by print (?maybe?) and str()
		return "y^2 = x^3 + "+str(self.a)+"x + "+str(self.b)+" (mod "+str(self.prime)+")"
	def __repr__(self):  # Over-ride string conversion used repr() for output which can be used as input
		return "EllipticCurve("+str(self.prime)+", ["+str(self.a)+","+str(self.b)+"])"


class EllipticCurveElt(object):
	"""EllipticCurveElt(ec,[x,y]) is an element of the elliptic curve ec, with coordinates (x,y) in
		affine Weierstrass form.
		Usage: 
			>>> from ellipticcurve import *
			>>> ec29 = EllipticCurve(29,[4,20]) # Define the elliptic curve over GF(29)
			>>> print ec29
				y^2 = x^3 + 4x + 20 (mod 29)
			>>> pt = EllipticCurveElt(ec29,[2,6])  # Define a point on the curve
			>>> print pt+pt, pt-pt  # Add the point to itself, then subtract (yielding zero pt on curve)
				(1, 5) (Infinity, Infinity)
			>>> print 5*pt  # Compute 5*pt = pt+pt+pt+pt+pt
				(0, 7)
			>>> repr(ec29)  # Format usable as input
				EllipticCurve(29, [4,20])
			>>> repr(pt)    # Format usable as input
				EllipticCurveElt(EllipticCurve(29, [4,20]), [2,6])
	"""
	def __init__(self, ellipticcurve, coords):
		self.ec = ellipticcurve
		self.x = coords[0]
		self.y = coords[1]
		if(self.ec.verbose): print "Create point {0} on elliptic Curve {1}".format(self,self.ec)
	def __format__(self,fmtspec):  # Over-ride format conversion
		if(self.x == "Infinity"):
			return "(Infinity,Infinity)"
		else:
			return "({1:{0}}, {2:{0}})".format(self.ec.fmtspec,self.x,self.y)
	def __str__(self):  # Over-ride string conversion used by print (?maybe?) and str()
		return "("+str(self.x)+", "+str(self.y)+")"
	def __repr__(self): # Over-ride string conversion used repr() for output which can be used as input
		return "EllipticCurveElt("+repr(self.ec)+", ["+repr(self.x)+","+repr(self.y)+"])"
	def add(self,summand):
		"""add elements of elliptic curves"""
		if (self.x == "Infinity"):  # Add to zero (i.e. point at infinity)
			if(self.ec.verbose): print "Add: 0 + ({0},{1}) = ({0},{1})".format(summand.x,summand.y)
			return summand
		elif (summand.x == "Infinity"):  # Add zero (i.e. point at infinity)
			if(self.ec.verbose): print "Add: ({0},{1}) + 0 = ({0},{1})".format(self.x,self.y)
			return self
		elif ((summand.x == self.x) and ((summand.y + self.y) % self.ec.prime == 0)): # P + (-P) = infty
			if(self.ec.verbose): print "Add: ({0},{1}) + ({2},{3}) = 0".format(self.x,self.y,summand.x,summand.y)
			return EllipticCurveElt(self.ec, ["Infinity","Infinity"])
		else:  # Usual addition and doubling (what a nuisance: lambda is a keyword - shorten to lamb)
			if (self.x == summand.x):  # Point doubling
				if(self.ec.verbose): print "Add (Double): ",
				lamb = (3*(self.x**2)+self.ec.a)*numbthy.xgcd(2*self.y,self.ec.prime)[1] % self.ec.prime
			else:  # Point addition
				if(self.ec.verbose): print "Add: ",
				lamb = (self.y - summand.y) * numbthy.xgcd((self.x - summand.x), self.ec.prime)[1]  % self.ec.prime
			if(self.ec.verbose): print "lambda = {0}; ".format(lamb),
			x3 = (lamb*lamb - self.x - summand.x) % self.ec.prime
			y3 = (lamb*(self.x-x3) - self.y) % self.ec.prime
			if(self.ec.verbose): print "({0},{1}) + ({2},{3}) = ({4},{5})".format(self.x,self.y,summand.x,summand.y,x3,y3)
			return EllipticCurveElt(self.ec, [x3,y3])
	def __add__(self,summand):   # Overload the "+" operator
		return self.add(summand)
	def __iadd__(self,summand): # Overload the "+=" operator
		self = self + summand
		return self
	def __neg__(self):  # Overload the "-" unary operator
		return EllipticCurveElt(self.ec, [self.x, ((-self.y) % self.ec.prime)])
	def __sub__(self,summand):  # Overload the "-" binary operator
		return self.__add__(-summand)
	def __isub__(self,summand): # Overload the "-=" operator
		self = self - summand
		return self
	def mult(self,multand):  # Multiply EC point by integer (repeated addition in EC)
		"""multiply elliptic curve point by integer (repeated addition in the elliptic curve)"""
		accum = EllipticCurveElt(self.ec, ["Infinity","Infinity"]) # start with identity
		i = 0
		bpow2 = self
		while ((multand>>i) > 0):
			if((multand>>i) & 1):
				accum = (accum + bpow2)
			bpow2 = (bpow2 + bpow2)
			i+=1
		return accum
	def __rmul__(self,multip):  # Overload the "*" operator
		return self.mult(multip)
	def __imul__(self,multip): # Overload the "*=" operator
		self = self.mult(multip)
		return self
################################################################################
# Example Curves: (e.g. [HMV04]
# y^2 = x^3 + x + 1 (mod 5), order 9, disc=4
#    P = [0,1] has full order
# y^2 = x^3 + 4x + 20 (mod 29), order 37, disc=4
#    P = [2,6] has full order
# Simple minded code:
# >>> p=5; a=1; b=1
# >>> disc = -16*(4*a**3+27*b**2) % p
# >>> disc
# 4
# >>> pt = [0,1]
# >>> def plus(p1,p2):
# ...     lamb=(p1[1]-p2[1])*xgcd(p1[0]-p2[0],p)[1] % p
# ...     x3 = (lamb**2 - p1[0] - p2[0]) % p
# ...     return [x3, (lamb*(p1[0]-x3) - p1[1]) % p]
# >>> def doub(p1):
# ...     lamb=(3*p1[0]**2+a)*xgcd(2*p1[1],p)[1] % p
# ...     x3 = (lamb**2 - 2*p1[0]) % p
# ...     return [x3, (lamb*(p1[0]-x3) - p1[1]) % p]
# >>> pt2 = doub(pt)
# >>> pt2 # 2*pt
# [4, 2]
# >>> plus(pt, pt2) # 3*pt
# [2, 1]
# >>> plus(pt, _) # 4*pt
# [3, 4]
# >>> plus(pt, _) # 5*pt (note: 4*pt = -5*pt)
# [3, 1]
############################################################################
# License: Freely available for use, abuse and modification
# (this is the Simplified BSD License, aka FreeBSD license)
# Copyright 2001-2014 Robert Campbell. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#    1. Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#
#    2. Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in 
#       the documentation and/or other materials provided with the distribution.
############################################################################
