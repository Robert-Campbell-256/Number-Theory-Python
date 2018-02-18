######################################################################################
# Elliptic Curve
# Elliptic curves in reduced Weierstrass form over prime order fields
# Author: Robert Campbell, <r.campbel.256@gmail.com>
# Date: 17 Feb, 2018
# Version 0.26
# License: Simplified BSD (see details at bottom)
######################################################################################
"""Elliptic Curve
	EllipticCurve(p,[a,b]) is the elliptic curve in affine restricted Weierstrass form
	y^2 = x^3 +ax + b (mod p)
	Usage: 
		>>> from ellipticcurve import *
		>>> ec29 = EllipticCurve(29,[4,20]); ec29
			y^2 = x^3 + 4x + 20 (mod 29)
		>>> pt = EllipticCurveElt(ec29,[2,6]); pt
			(2, 6)
		>>> 7*pt
			(3, 28)
		>>> 37*pt
			(Infinity, Infinity)
		>>> pt1 = 9*pt; pt1 - pt
			(15, 27)
		>>> '{0:f}'.format(ec29) # Full format
			'EllipticCurve(29, (4,20))'
"""   
   
__version__ = '0.26' # Format specified in Python PEP 396
Version = 'ELLIPTICCURVE.PY, version ' + __version__ + ', 17 Feb, 2018, by Robert Campbell, <r.campbel.256@gmail.com>'

import numbthy # For xgcd (for modinv) and sqrtmod
import random   # Generate random elements
import sys      # Check Python2 or Python3
import math  # For sqrt

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
		>>> ec29 = EllipticCurve(29,[4,20]); ec29
			y^2 = x^3 + 4x + 20 (mod 29)
		>>> pt = EllipticCurveElt(ec29,[2,6]); pt
			(2, 6)
		>>> 7*pt
			(3, 28)
		>>> 37*pt
			(Infinity, Infinity)
		>>> pt1 = 9*pt; pt1 - pt
			(15, 27)
		>>> '{0:f}'.format(ec29) # Full format
			'EllipticCurve(29, (4,20))'
	"""   
	def __init__(self,prime,coeffs,fmtspec="s"):
		self.prime = prime
		if(not(numbthy.isprime(self.prime))): raise ValueError("***** Error *****: Characteristic of base field {0} must be prime".format(self.prime))
		self.a = coeffs[0]
		self.b = coeffs[1]
		self.discriminant = -16*(4*(self.a**3)+27*(self.b**2)) % self.prime
		if(self.discriminant == 0): raise ValueError("***** Error *****: Not an elliptic curve - Zero discriminant (-16*(4*({0}^3)+27*({1}^2)))".format(self.a,self.b))
		self.fmtspec = fmtspec

	def isIntType(self,x):
		if sys.version_info < (3,): return isinstance(x,(int, long,))
		else: return isinstance(x,(int,))
                
	def __call__(self,pt):  # Coerce constant or array of coeffs as elt of field
		"""Create a point on the curve from a tuple or list of integers. (not [Infinity,Infinity])"""
		if not (isinstance(pt,(list,tuple,)) and len(pt)==2 and self.isIntType(pt[0]) and self.isIntType(pt[1])):
			raise ValueError('{0} should be a list or tuple of two integers'.format(pt))
		if not ((pow(pt[0],3,self.prime) + self.a*pt[0] + self.b - pt[1]*pt[1]) % self.prime == 0):
			raise ValueError('{0} is not a point on the curve {1}'.format(pt,self))
		return EllipticCurveElt(self,pt)

	def __iter__(self):
		"""Generator producing all points on the elliptic curve."""
		yield EllipticCurveElt(self, ("Infinity","Infinity"))
		x = 0
		for x in range(self.prime):
			ysq = (pow(x,3,self.prime) + self.a*x + self.b) % self.prime
			if((ysq == 0) or (pow(ysq,(self.prime-1)//2,self.prime)==1)):
				if (ysq == 0): y = 0
				else: y = numbthy.sqrtmod(ysq,self.prime)
				if((y % 2)==1): y = self.prime - y # Always even y first (consistent order)
				yield EllipticCurveElt(self, (x,y))
				if (y != 0): yield EllipticCurveElt(self, (x,self.prime - y)) # Distinct unless y==0
		raise StopIteration
                
	def random_element(self):
		"""A random element of the elliptic curve."""
		# Currently, choosing point at infinity (group identity) and point
		# with y=0 is twice as likely as any other point.
		# Find a random x such that y^2 = x^3 + ax + b has a solution (mod p)
		xrand = random.randint(-1,self.prime-1)
		if(xrand == -1): return EllipticCurveElt(self, ("Infinity","Infinity"))
		ysq = (pow(xrand,3,self.prime) + self.a*xrand + self.b) % self.prime
		while((ysq != 0) and (pow(ysq,(self.prime-1)//2,self.prime)!=1)):
			xrand = random.randint(-1,self.prime-1)
			if(xrand == -1): return EllipticCurveElt(self, ("Infinity","Infinity"))
			ysq = (pow(xrand,3,self.prime) + self.a*xrand + self.b) % self.prime
		# Given x, find a y solving y^2 = x^3 + ax + b (mod p)
		if (ysq == 0): yrand = 0
		else: yrand = numbthy.sqrtmod(ysq,self.prime)
		if(random.randint(0,1)==1): yrand = self.prime - yrand # Choose between pt and -pt
		return EllipticCurveElt(self,(xrand,yrand))
		
	def __format__(self,fmtspec):  # Over-ride format conversion
		"""Override the format when outputting an elliptic curve.
		A default can be set when the curve is defined or it can be specified for each output.
		Possible formats are:
		        s - short format (default)
		        f - full format, can be used as input
		        t - LaTeX format"""
		if(fmtspec == ''): fmtspec = self.fmtspec
		if(fmtspec == 's'): # Short format
			return "y^2 = x^3 + {0}x + {1} (mod {2})".format(self.a,self.b,self.prime)
		if(fmtspec == 'f'): # Full format
			return "EllipticCurve({0},({1},{2}))".format(self.prime,self.a,self.b)
		if(fmtspec == 't'): # LaTeX format
			return "\mathbb{{E}}_{{y^2 = x^3 + {0}x + {1} \pmod{{{2}}}}}".format(self.a,self.b,self.prime)
                        
	def __str__(self):   # Over-ride string conversion used by print (?maybe?) and str()
		return format(self)
	def __repr__(self):  # Over-ride string conversion for output
		return format(self)


class EllipticCurveElt(object):
	"""EllipticCurveElt(ec,[x,y]) is an element of the elliptic curve ec, with coordinates (x,y) in
		affine Weierstrass form.
	Usage: 
		>>> from ellipticcurve import *
		>>> ec29 = EllipticCurve(29,[4,20]); ec29
			y^2 = x^3 + 4x + 20 (mod 29)
		>>> pt = EllipticCurveElt(ec29,[2,6]); pt
			(2, 6)
		>>> 7*pt
			(3, 28)
		>>> 37*pt
			(Infinity, Infinity)
		>>> pt1 = 9*pt; pt1 - pt
			(15, 27)
		>>> '{0:f}'.format(pt) # Full format
			'EllipticCurveElt(EllipticCurve(29,(4,20)), (2,6))'
	"""
	def __init__(self, ellipticcurve, coords):
		self.ec = ellipticcurve
		self.x = coords[0]
		self.y = coords[1]
	def __format__(self,fmtspec):  # Over-ride format conversion
		"""Override the format when outputting a point on an elliptic curve.
		A default can be set when the curve is defined or it can be specified for each output.
		Possible formats are:
		        s - short format (default)
		        f - full format, can be used as input
		        t - LaTeX format"""
		if(fmtspec == 'f'):
			if(self.x == "Infinity"):
				return "EllipticCurveElt("+'{0:f}'.format(self.ec)+", (Infinity,Infinity))"
				return "(Infinity,Infinity)"
			else:
				return "EllipticCurveElt("+'{0:f}'.format(self.ec)+", ("+format(self.x)+","+format(self.y)+"))"
		else:  # Both short and LaTeX formats
			if(self.x == "Infinity"):
				if(fmtspec == 't'): # LaTeX format
					return "(\infty,\infty)"
				else: # short format
					return "(Infinity,Infinity)"
			else:
				return "({0}, {1})".format(self.x,self.y)
	def __str__(self):  # Over-ride string conversion used by str()
		return format(self)
	def __repr__(self): # Over-ride string conversion for output
		return format(self)
	def __cmpec__(self,other): # Implement cmp for both Python2 and Python3
		"""compare two points for equality and (possibly in future allow sorting)
		overloaded to allow comparisons to lists of integers"""
		# Coerce if comparing list (x,y) and point
		if (isinstance(other,(list,tuple,)) and len(other)==2 and self.ec.isIntType(pt[0]) and self.ec.isIntType(pt[1])):
			if (other[0]==self.x) and (other[1]==self.y): return 0
			else: return 1
		elif(self.ec != other.ec):
			raise ValueError("Cannot compare elements of different elliptic curves: <{0}> and <{1}>".format(self.ec,other.ec))
		else:
			if (other.x==self.x) and (other.y==self.y): return 0
			else: return 1
	def __eq__(self,other): return (self.__cmpec__(other) == 0)
	def __ne__(self,other): return (self.__cmpec__(other) != 0)
	def add(self,summand):
		"""add elements of elliptic curves"""
		if (self.x == "Infinity"):  # Add to zero (i.e. point at infinity)
			return summand
		elif (summand.x == "Infinity"):  # Add zero (i.e. point at infinity)
			return self
		elif ((summand.x == self.x) and ((summand.y + self.y) % self.ec.prime == 0)): # P + (-P) = infty
			return EllipticCurveElt(self.ec, ("Infinity","Infinity"))
		else:  # Usual addition and doubling (what a nuisance: lambda is a keyword - shorten to lamb)
			if (self.x == summand.x):  # Point doubling
				lamb = (3*(self.x**2)+self.ec.a)*numbthy.xgcd(2*self.y,self.ec.prime)[1] % self.ec.prime
			else:  # Point addition
				lamb = (self.y - summand.y) * numbthy.xgcd((self.x - summand.x), self.ec.prime)[1]  % self.ec.prime
			x3 = (lamb*lamb - self.x - summand.x) % self.ec.prime
			y3 = (lamb*(self.x-x3) - self.y) % self.ec.prime
			return EllipticCurveElt(self.ec, (x3,y3))
	def __add__(self,summand):   # Overload the "+" operator
		return self.add(summand)
	def __iadd__(self,summand): # Overload the "+=" operator
		self = self + summand
		return self
	def __neg__(self):  # Overload the "-" unary operator
		return EllipticCurveElt(self.ec, (self.x, ((-self.y) % self.ec.prime)))
	def __sub__(self,summand):  # Overload the "-" binary operator
		return self.__add__(-summand)
	def __isub__(self,summand): # Overload the "-=" operator
		self = self - summand
		return self
	def mult(self,multand):  # Multiply EC point by integer (repeated addition in EC)
		"""multiply elliptic curve point by integer (repeated addition in the elliptic curve)"""
		accum = EllipticCurveElt(self.ec, ("Infinity","Infinity")) # start with identity
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

############################################################################
# License: Freely available for use, abuse and modification
# (this is the Simplified BSD License, aka FreeBSD license)
# Copyright 2001-2018 Robert Campbell. All rights reserved.
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
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
# SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################
# 4 Feb 2018: ver 0.25
#   Formatting - for EllipticCurve and EllipticCurveElt
#   Added document strings
#   Remove verbose mode
#   Fixed support for Python 3
# 17 Feb 2018: ver 0.26
#   Added random_element
#   Added iterator
#   Added call (coerce list as point on curve)
#   Changed list to tuple for (x,y) - fix comparison bug
