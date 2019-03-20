######################################################################################
# GAUSSINT.PY
# Basic Number Theory functions implemented in Python
# Note: Currently requires Python 3.x (uses floordiv, changes to the "types" module…)
# Author: Robert Campbell, <r.campbel.256@gmail.com>
# Modified: Hubert Holin, <Hubert.Holin.1982@Polytechnique.org>
# Date: 17 March, 2019
# Version 1.2
# License: Simplified BSD (see details at bottom)
# Requirements:
#	Requires at least Python 3.x (runs fine on Python 3.6)
# Bugs:
#	None currently known.
######################################################################################
__version__ = '1.2' # Format specified in Python PEP 396
Version = 'GAUSSINT.PY, version ' + __version__ +\
	', 8 June, 2013, by Robert Campbell, <r.campbel.256@gmail.com>'+\
	', modified 17 March 2019 by Hubert Holin, <Hubert.Holin.1982@Polytechnique.org>'


import math	# For tools used in primality testing


class GaussInt:
	"""Gaussian Integer functions.
	Functions implemented are:
	     Arithmetic functions: +,*,//,%,**(exponentiation)
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
	
	
	def __init__(self, a = 0, b = 0):
		
		if (type(a) is complex):
			
			if b != 0:
				
				raise TypeError("Attempting to ceate a Gauss Integer from a complex "+
					"number and another input ({0:s} and {1:s})!".format(a, b))
			
			self.r = round(a.real)
			self.i = round(a.imag)
		
		else:
			
			self.r = int(a)
			self.i = int(b)
	
	
	def __str__(self):  # Overload string conversion used by print
		
		return "(" + str(self.r) + ((" + "+str(self.i)) if (self.i >= 0) else (" - "+str(-self.i))) + " i)"
	
	
	def __format__(self, spec):# Overload string conversion used by format
		
		return "(" + str(self.r) + ((" + "+str(self.i)) if (self.i >= 0) else (" - "+str(-self.i))) + " i)"
	
	
	def __repr__(self):  # Overload conversion used for output
		
		return "GaussInt(" + str(self.r) + ", " + str(self.i) + ")"
	
	
	def __complex__(self):  # Allow conversion to complex type
		return complex(self.r, self.i)
	
	
	def __eq__(self,other):  # Overload the "==" test operator - NOTE: differs from version 1.1
		
		if (type(other) is not GaussInt):
			
			return False
		
		else:
			
			return (self.r == other.r) and (self.i == other.i)
	
	
	def __ne__(self,other):  # Overload the "!=" test operator
		
		return not (self == other)
	
	
	def neutral_element_for_multiplication():
		
		return __class__(1)
	
	
	def conjugate(self):
		
		return GaussInt(self.r, -self.i)
	
	
	def norm(self):
		
		return self.r*self.r + self.i*self.i
	
	
	def __pos__(self):	  # Overload the "+" unary operator 
		
		return self
	
	
	def add(self,summand):
		
		sum_r = self.r + summand.r
		sum_i = self.i + summand.i
		
		return GaussInt(sum_r, sum_i)
	
	
	def __add__(self,summand):  # Overload the "+" binary operator
		
		if type(summand) is int:
			
			return GaussInt(self.r+summand, self.i)
		
		else:
			
			return self.add(summand)
	
	
	def __radd__(self,summand):  # Overload the "+" binary operator
		
		if type(summand) is int:
			
			return GaussInt(self.r+summand, self.i)
		
		else:
			
			return self.add(summand)
	
	
	def __iadd__(self,summand): # Overload the "+=" operator
		
		self = self + summand
		
		return self
	
	
	def __neg__(self):  # Overload the "-" unary operator 
		
		return GaussInt(-self.r,-self.i)
	
	
	def __sub__(self,summand):  # Overload the "-" binary operator
		
		return self.__add__(-summand)
	
	
	def __rsub__(self,summand):  # Overload the "-" binary operator
		
		if type(summand) is int:
			
			return GaussInt(summand-self.r, -self.i)
		
		else:
			
			return summand-self
	
	
	def __isub__(self,summand): # Overload the "-=" operator
		
		self = self - summand
		
		return self
	
	
	def mult(self,multip):
		
		prod_r = (self.r * multip.r) - (self.i * multip.i)
		prod_i = (self.i * multip.r) + (self.r * multip.i)
		
		return GaussInt(prod_r, prod_i)
	
	
	def __mul__(self,multip):  # Overload the "*" operator
		
		if type(multip) is int:
			
			return GaussInt(self.r*multip, self.i*multip)
		
		else:
			
			return self.mult(multip)
	
	
	def __rmul__(self,multip):  # Overload the "*" operator
		
		if type(multip) is int:
			
			return GaussInt(self.r*multip, self.i*multip)
		
		else:
			
			return self.mult(multip)
	
	
	def __imul__(self,multip): # Overload the "*=" operator
		
		self = self * multip
		
		return self
	
	
	def floordiv(self,divisor):
		
		if type(divisor) is int:
			
			numerator = (-self if (divisor < 0) else self)
			
			denominator = (-divisor if (divisor < 0) else divisor)
			
			if denominator == 0:
				
				raise ZeroDivisionError("{0:s} is null!".format(divisor))
		
		else:
			
			numerator = self*divisor.conjugate()
			
			denominator = divisor.norm()	# Recall that denominator >= 0
			
			if denominator == 0:
				
				raise ZeroDivisionError("{0:s} is null!".format(divisor))
		
		candidate_r = numerator.r//denominator
		candidate_i = numerator.i//denominator
		
		# i.e. (candidate_r+1)*denominator-numerator.r < numerator.r-candidate_r*denominator
		if (2*candidate_r+1)*denominator < 2*numerator.r:
			
			candidate_r += 1
		
		# i.e. (candidate_i+1)*denominator-numerator.i < numerator.i-candidate_i*denominator
		if (2*candidate_i+1)*denominator < 2*numerator.i:
			
			candidate_i += 1
		
		return GaussInt(candidate_r,candidate_i)
	
	
	def __floordiv__(self,divisor):  # Overload the "//" operator
		
		return self.floordiv(divisor)
	
	
	def __ifloordiv__(self,divisor): # Overload the "//=" operator
		
		self = self//divisor
		
		return self
	
	
	def mod(self,divisor):
		
		return self - divisor * (self//divisor)
	
	
	def __mod__(self,divisor):  # Overload the "%" operator
		
		return self.mod(divisor)
	
	
	def __imod__(self,divisor): # Overload the "%=" operator
		
		self = self % divisor
		
		return self
	
	
	def divmod(self,divisor):
		
		q = self//divisor
		
		return q, self - divisor * q
	
	
	def xgcd(self,other):
		
		quot = GaussInt()
		
		a1 = GaussInt(1,0)
		b1 = GaussInt(0,0)
		
		a2 = GaussInt(0,0)
		b2 = GaussInt(1,0)
		
		a = self
		b = other
		
		if(b.norm() > a.norm()):	# Need to start with a>b
			
			a,b = b,a					# Swap a and b
			a1,b1,a2,b2 = a2,b2,a1,b1	# Swap (a1,b1) with (a2,b2)
		
		while (True):
			
			quot = a // b
			
			a %= b
			
			a1 -= quot*a2
			b1 -= quot*b2
			
			if (a == GaussInt(0,0)):
				
				return b, a2, b2
			
			quot = b // a
			
			b %= a
			
			a2 -= quot*a1
			b2 -= quot*b1
			
			if (b == GaussInt()):
				
				return a, a1, b1
	
	
	def Bézout(self, other):
		
		a = self
		b = other
		
		if a.norm() < b.norm():
			
			(u, v, pgcd) = b.Bézout(a)
			
			return (v, u, pgcd)
		
		if b == 0:
			
			return (1, 0, a)
		
		u_n, u_n_moins_1, v_n, v_n_moins_1 = 0, 1, 1, 0
		
		while b.norm() > 0:
		
			q,r = a.divmod(b)
			
			u_n_plus_1 = u_n_moins_1 - q*u_n
			v_n_plus_1 = v_n_moins_1 - q*v_n
			
			a, b = b, r
			u_n_moins_1, u_n, v_n_moins_1, v_n = u_n, u_n_plus_1, v_n, v_n_plus_1
		
		return (u_n_moins_1, v_n_moins_1, a)
	
	
	def gcd(self,other):
		
		a = self
		b = other
		
		if a.norm() < b.norm():
			
			return b.gcd(a)
	
		while b.norm() > 0:
		
			q,r = a.divmod(b)
			a,b = b,r
	
		return a
	
	
	def powmod(self, a_power, a_modulus):
	# We adapt the Binary Exponentiation algorithm with modulo
		
		result = GaussInt(1)
		
		auxilliary = GaussInt(self.r, self.i)
		
		while a_power:
			
			if a_power % 2:	# If power is odd
				
				result = (result * auxilliary) % a_modulus
			
			# Divide the power by 2
			a_power >>= 1
			
			# Multiply base to itself
			auxilliary = (auxilliary * auxilliary) % a_modulus
		
		return result
	
	
	def __pow__(self, a_power):  # Overload the "**" operator
	# We adapt the Binary Exponentiation algorithm (without modulo!)
		
		result = GaussInt(1)
		
		auxilliary = GaussInt(self.r, self.i)
		
		while a_power:
			
			if a_power % 2:	# If power is odd
				
				result = result * auxilliary
			
			# Divide the power by 2
			a_power >>= 1
			
			# Multiply base to itself
			auxilliary = auxilliary * auxilliary
		
		return result
	
	
	def isprime(self):
		"""n.isprime() - Test whether the GaussInt n is prime using a variety of pseudoprime tests."""
		# Multiply by (1,i,-1,-i) to rotate to first quadrant (similar to abs)
		if (self.r < 0): self *= (-1)
		if (self.i < 0): self *= GaussInt(0,1)
		# Check some small non-primes
		if (self in [GaussInt(0,0), GaussInt(1,0), GaussInt(0,1)]): return False
		# Check some small primes
		if (self in [GaussInt(1,1), GaussInt(2,1), GaussInt(1,2), GaussInt(3,0), GaussInt(0,3), GaussInt(3,2), GaussInt(2,3)]):
			return True
		return self.isprimeF(2) and self.isprimeF(3) and self.isprimeF(5)
	
	
	def isprimeF(self,base):
		"""n.isprimeF(base) - Test whether the GaussInt n is prime using the
		Gaussian Integer analogue of the Fermat pseudoprime test."""
		if type(base) is not GaussInt:
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
	
	
	def factorPR(self):	# TODO: learn and test
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

