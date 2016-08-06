######################################################################################
# GAUSSINT.PY
# Basic Number Theory functions implemented in Python
# Note: Currently requires Python 2.x (uses +=, %= and other 2.x-isms)
# Author: Robert Campbell, <r.campbel.256@gmail.com>
# Date: 8 June, 2013
# Version 1.1
# License: Simplified BSD (see details at bottom)
# Requirements:
#	Requires at least Python 2.x (runs fine on Python 2.2)
# Bugs:
#	Division relies on floating point (complex) computations - fails for large values
#		(Need to reimplement, even at the cost of speed)
#   Need to clean up input/output routines, so a-bi does not print a+-bi and so
#       one can easily cut and paste to input previous output
######################################################################################
__version__ = '1.1' # Format specified in Python PEP 396
Version = 'GAUSSINT.PY, version ' + __version__ + ', 8 June, 2013, by Robert Campbell, <r.campbel.256@gmail.com>'

import math  # Use floor
import types # Use IntType, LongType

class GaussInt:
	"""Gaussian Integer functions.
	Functions implemented are:
	     Arithmetic functions: +,*,/,%,**(exponentiation)
	     a.gcd(b) - Compute the greatest common divisor of a and b.
	     a.xgcd(b) - Extended gcd - return gcd(a,b) and x,y such that gcd(a,b)=xa+yb.
	     n.isprime() - Is n prime (pseudoprime tests)
	     n.factor() - Return a factor of n.
	     n.factors() - Return list of the factors of n.
	Gaussian Integers can be created by:
	     n = GaussInt(5,7)  # Create (5 + 7i)
	     n = GaussInt(13)  # Create (5 + 0i)
	     z = complex(2,3); n = GaussInt(z) # Round the complex number to integer form
	A list of the functions implemented in GaussInt is printed by the command help(GaussInt).
	Usage: from gaussint import * """
   	def __init__(self, a=0, b=0):
		if (type(a) == type(complex(1,0))) and (b == 0):
		      b = int(a.imag + 0.5); a = int(a.real + 0.5)
		self.r = int(a)
		self.i = int(b)
	def __str__(self):  # Overload string conversion used by print
		return "(" + str(self.r) + " + " + str(self.i) + "i)"
	def __repr__(self):  # Overload conversion used for output
		return "GaussInt(" + str(self.r) + ", " + str(self.i) + ")"
	def __complex__(self):  # Allow conversion to complex type
		return complex(self.r, self.i)
	def __eq__(self,other):  # Overload the "==" test operator
		if (type(other) != type(GaussInt(1,0))):
		     other = GaussInt(other)  # Coerce if base not GaussInt (works for int or complex)
		return (self.r == other.r) and (self.i == other.i)
	def __ne__(self,other):  # Overload the "!=" test operator
		return not (self == other)
	def norm(self):
		return self.r*self.r + self.i*self.i
	def add(self,summand):
		sum = GaussInt()
		sum.r = self.r + summand.r
		sum.i = self.i + summand.i
		return sum
	def __add__(self,summand):  # Overload the "+" operator
		if ((type(summand) == types.IntType) or (type(summand) == types.LongType)):
			summand = GaussInt(summand)  # Coerce if adding integer and GaussInt
		return self.add(summand)
	def __radd__(self,summand):  # Overload the "+" operator
		if ((type(summand) == types.IntType) or (type(summand) == types.LongType)):
			summand = GaussInt(summand)  # Coerce if adding integer and GaussInt
		return self.add(summand)
	def __iadd__(self,summand): # Overload the "+=" operator
		self = self + summand
		return self
	def __neg__(self):  # Overload the "-" unary operator 
		return GaussInt(-self.r,-self.i)
	def __sub__(self,summand):  # Overload the "-" binary operator 
		return self.__add__(-summand)
	def __isub__(self,summand): # Overload the "-=" operator
		self = self - summand
		return self
	def mult(self,multip):
		prod = GaussInt()
		prod.r = (self.r * multip.r) - (self.i * multip.i)
		prod.i = (self.i * multip.r) + (self.r * multip.i)
		return prod
	def __mul__(self,multip):  # Overload the "*" operator
		if ((type(multip) == types.IntType) or (type(multip) == types.LongType)):
			multip = GaussInt(multip)  # Coerce if multiplying integer and GaussInt
		return self.mult(multip)
	def __rmul__(self,multip):  # Overload the "*" operator
		if ((type(multip) == types.IntType) or (type(multip) == types.LongType)):
			multip = GaussInt(multip)  # Coerce if multiplying integer and GaussInt
		return self.mult(multip)
	def __imul__(self,multip): # Overload the "*=" operator
		self = self * multip
		return self
	def div(self,divisor):
		if (type(divisor) in [types.IntType, types.LongType]):
			divisor = GaussInt(divisor)  # Coerce if dividing GaussInt by integer
		q = complex(self.r,self.i)/complex(divisor.r,divisor.i)
		return GaussInt(int(round(q.real)),int(round(q.imag)))
	def __div__(self,divisor):  # Overload the "/" operator 
		return self.div(divisor)
	def __idiv__(self,divisor): # Overload the "/=" operator
		self = self/divisor
		return self
	def mod(self,divisor):
		if (type(divisor) in [types.IntType, types.LongType]):
			divisor = GaussInt(divisor)  # Coerce if dividing GaussInt by integer
		return self - divisor * (self/divisor)
	def __mod__(self,divisor):  # Overload the "%" operator 
		return self.mod(divisor)
	def __imod__(self,divisor): # Overload the "%=" operator
		self = self%divisor
		return self
	def divmod(self,divisor):
		if (type(divisor) in [types.IntType, types.LongType]):
			divisor = GaussInt(divisor)  # Coerce if dividing GaussInt by integer
		q = self/divisor
		return q, self - divisor * q
	def xgcd(self,other):
		quot=GaussInt(); a1=GaussInt(1,0); b1=GaussInt(0,0); a2=GaussInt(0,0); 
		b2=GaussInt(1,0); a = self; b = other;
		if(b.norm() > a.norm()):
			a,b = b,a  # Swap a and b - need to start with a>b
			a1,b1,a2,b2 = a2,b2,a1,b1 # Swap (a1,b1) with (a2,b2)
		while (1):
			quot = a / b
			a %= b
			a1 -= quot*a2; b1 -= quot*b2
#			print a,b,a1,b1,a2,b2
			if (a == GaussInt(0,0)):
				return b, a2, b2
			quot = b / a
			b %= a
			a2 -= quot*a1; b2 -= quot*b1
#			print a,b,a1,b1,a2,b2
			if (b == GaussInt()):
				return a, a1, b1
	def gcd(self,other):
		return self.xgcd(other)[0]
	def powmod(self,exp,mod):
		accum = GaussInt(1,0); basepow2 = self; i=0
		while ((exp>>i)>0):
			if (((exp>>i) & 1) == 1):
				accum = (accum*basepow2) % mod
			basepow2 = (basepow2*basepow2) % mod
			i+=1
		return accum
	def pow(self,exp):
		accum = GaussInt(1,0); basepow2 = self; i=0
		while ((exp>>i)>0):
			if (((exp>>i) & 1) == 1):
				accum = (accum*basepow2)
			basepow2 = (basepow2*basepow2)
			i+=1
		return accum
	def __pow__(self,exponent):  # Overload the "**" operator 
		return self.pow(exponent)
	def isprime(self):
		"""n.isprime() - Test whether the GaussInt n is prime using a variety of pseudoprime tests."""
		# Multiply by (1,i,-1,-i) to rotate to first quadrant (similar to abs)
		if (self.r < 0): self *= (-1)
		if (self.i < 0): self *= GaussInt(0,1)
		# Check some small non-primes
		if (self in [GaussInt(0,0), GaussInt(1,0), GaussInt(0,1)]): return False
		# Check some small primes
		if (self in [GaussInt(1,1), GaussInt(2,1), GaussInt(1,2), GaussInt(3,0), GaussInt(0,3), GaussInt(3,2), GaussInt(2,3)]): return True
		return self.isprimeF(2) and self.isprimeF(3) and self.isprimeF(5)
	def isprimeF(self,base):
		"""n.isprimeF(base) - Test whether the GaussInt n is prime using the Gaussian Integer analogue of the Fermat pseudoprime test."""
		if (type(base) != type(GaussInt(1,0))):
		     base = GaussInt(base)  # Coerce if base not GaussInt (works for int or complex)
		return base.powmod(self.norm()-1,self) == GaussInt(1,0)
	# Note: Possibly more effective would be to use the characterization of primes
	# in the Gaussian Integers based on the primality of their norm and reducing mod 4.
	# This depends on the characterization of the ideal class group, and only works for
	# simple rings of algebraic integers.
	def factor(self):
		"""n.factor() - Find a prime factor of Gaussian Integer n using a variety of methods."""
		if (self.isprime()): return n
		for fact in [GaussInt(1,1), GaussInt(2,1), GaussInt(1,2), 
			     GaussInt(3,0), GaussInt(3,2), GaussInt(2,3)]:
			if self%fact == 0: return fact
		return self.factorPR()  # Needs work - no guarantee that a prime factor will be returned
	def factors(self):
		"""n.factors() - Return a sorted list of the prime factors of Gaussian Integer n."""
		if (self.isprime()):
		      return [self]
		fact = self.factor()
		if (fact == 1): return "Unable to factor "+str(n)
		facts = (self/fact).factors() + fact.factors()
		return facts

	def factorPR(self):
		"""n.factorPR() - Find a factor of Gaussian Integer n using the analogue of the Pollard Rho method.  
		Note: This method will occasionally fail."""
		for slow in [2,3,4,6]:
			numsteps=2*math.floor(math.sqrt(math.sqrt(self.norm()))); fast=slow; i=1
		        while i<numsteps:
				slow = (slow*slow + 1) % self
				i = i + 1
				fast = (fast*fast + 1) % self
				fast = (fast*fast + 1) % self
				g = gcd(fast-slow,self)
				if (g != 1):
					if (g == self):
						break
					else:
						return g
			return 1
	# Note: Possibly more effective would be to factor the norm and then use the known
	# splitting properties of primes over the Gaussian Integers. This depends on the 
	# characterization of the ideal class group, and only works for simple rings of algebraic integers.



# >>> from GaussInt import *
# >>> a = GaussInt(1,0)

############################################################################
# License: Freely available for use, abuse and modification
# (this is the Simplified BSD License, aka FreeBSD license)
# Copyright 2001-2013 Robert Campbell. All rights reserved.
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

