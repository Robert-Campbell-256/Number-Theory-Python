######################################################################################
# FINITEFIELD.PY
# A finite field of prime power order
# Note: Version 0.9 - major changes to improve SAGE compatibility
#    Version 0.8 is compatible with both Python 2.5+ and 3.x 
#       and changes some names to improve SAGE compatibility
# Author: Robert Campbell, <r.campbel.256@gmail.com>
# Date: 16 Sept, 2017
# Version 0.93
# License: Simplified BSD (see details at bottom)1
######################################################################################
"""Finite fields.
	Functions implemented are:
		+, *, /, **, norm, trace, gen, random, enumerate
	Simple Examples:
	    >>> from finitefield import *     # names do not need path
	    >>> GF81 = GF(3**4)               # create finite field
	    >>> a = GF81([0,1])               # create element [0,1] = (0 + 1*x) = x
	    >>> '{0}^10 has value {1}'.format(a,a**10) # format as polynomial
	    >>> GF81.random().verbstr()       # a random value, formatted as polynomial
	    >>> a+2*(a**12)+2*(a**50)         # some arithmetic
	    >>> for i,x in enumerate(GF81): print i,x.verbstr()
	    >>> GF13 = GF(13); GF13(10)+GF13(10) # Compute mod 13

	Examples:
	    >>> from finitefield import *     # names do not need path
	    >>> GF9 = FiniteField(3,[2,1])    # Define GF(3^2), polys w/ GF3 coeffs, mod 2+x+x^2 (x^2 is assumed)
	    >>> a = FiniteFieldElt(GF9,[1,2]) # Define 1+2x in GF(9)
	    >>> a**12                         # Compute (2x+1)**12 in GF(9)
	  Define GF(5^8), defined mod z^8+3z^5+z^4+z^2+3z+4
	  Providing the factored order as (5^8-1) = (2^5)(3)(13)(313)
	  Output is coefficients only ('c'), not full polynomials ('z')
	    >>> GF5e8 = FiniteField(5,[4,3,1,0,1,3,0,0],'z',[[2,5],[3,1],[13,1],[313,1]],'p')
	    >>> a = FiniteFieldElt(GF5e8,[0,1])
	    >>> '{0}'.format(a**20)
	        '(2 + 3z^2 + 4z^3 + z^5 + z^6)'
	    >>> print a**20
	        [2, 0, 3, 4, 0, 1, 1, 0]
	    >>> a**20
	        FiniteFieldElt(FiniteField(5, [4, 3, 1, 0, 1, 3, 0, 0]),[2, 0, 3, 4, 0, 1, 1, 0])"""


__version__ = '0.93' # Format specified in Python PEP 396
Version = 'finitefield.py, version ' + __version__ + ', 16 Sept, 2017, by Robert Campbell, <r.campbel.256@gmail.com>'

import numbthy  # Use factor
import random   # Generate random elements
import types # Use IntType, LongType
from operator import add,mul,mod # To allow reduce(add,list) construct for sum

class FiniteField(object):
	"""Finite fields of prime power order.
	Driving polynomial must be monic and top coeff (i.e. 1) is implicit.
	Usage: 
	    >>> from finitefield import *
	    >>> GF9 = FiniteField(3,[2,1])    # Define GF(3^2), polys w/ GF3 coeffs, mod x^2+x+2
	    >>> a = FiniteFieldElt(GF9,[1,2]) # Define 2x+1 in GF(9)
	    >>> a**12                         # Compute (2x+1)**12 in GF(9)"""

	def __init__(self, prime, poly, var='x', orderfacts=None, fmtspec="p"):
		"""FiniteField(prime, poly, var='x', orderfacts=None, fmtspec="p")
		Create a finite field of order p**d, where d is the degree of the polynomial.
		Driving polynomial must be monic and top coeff (i.e. 1) is implicit.
		Example: 
		    >>> from finitefield import *
		    >>> GF9 = FiniteField(3,[2,1])    # Define GF(3^2), polys w/ GF3 coeffs, mod x^2+x+2
		    >>> a = FiniteFieldElt(GF9,[1,2]) # Define 2x+1 in GF(9)
		  Define GF(5^8), defined mod z^8+3z^5+z^4+z^2+3z+4
		  Providing the factored order as (5^8-1) = (2^5)(3)(13)(313)
		  Default output is coefficients only ('c'), not full polynomials ('z')
		    >>> GF5e8 = FiniteField(5,[4,3,1,0,1,3,0,0],'z',((2,5),(3,1),(13,1),(313,1)),'c')
		    >>> '{0}'.format(c**20)   # Default format
		       '20340110'
		    >>> '{0:p}'.format(c**20) # Force polynomial format
		       '(2 + 3z^2 + 4z^3 + z^5 + z^6)'"""

		self.char = prime
		self.degree = len(poly)
		self.order = self.char**self.degree
		self.modpoly = poly
		self.var = var
		self.fmtspec = fmtspec # p=polynomial; c=coeffsonly
		if(orderfacts == None):
			self.facts_order_gpunits = numbthy.factor(self.order - 1)
		else:
			self.facts_order_gpunits = orderfacts
			if((self.char**self.degree-1) != reduce(lambda theprod,primepow:theprod*primepow, [prime**thepow for [prime,thepow] in orderfacts])):
				   raise ValueError('{0} is not a factorization of ({1}^{2}-1)'.format(orderfacts,self.char,self.degree))
		self.reduc_table = [[0 for j in range(self.degree)] for i in range(2*self.degree-1)]
		for i in range(self.degree): self.reduc_table[i][i] = 1
		if(self.degree > 1):
			self.reduc_table[self.degree] = [(-self.modpoly[j])%self.char for j in range(self.degree)]
			for i in range(self.degree+1,2*self.degree-1):
				for j in range(self.degree):
					self.reduc_table[i][j] = sum(map(lambda k: (-self.modpoly[k]*self.reduc_table[i-self.degree+k][j]), range(self.degree))) % self.char

	def verbstr(self): # Requires feature from python 2.5.2 or better
		if(self.degree > 1):
			return "Z_"+str(self.char)+"["+self.var+"]/<"+self.polyprint(self.modpoly+[1],self.var)+">"
		else:
			return "Z_"+str(self.char)
			
	def __format__(self,fmtspec):  # Over-ride format conversion
		return str(self)  # Get to this later

	def __str__(self):  # Over-ride string conversion used by print
		if(self.degree > 1):
			return "GF("+str(self.char)+"^"+str(self.degree)+")"
		else:
			return "GF("+str(self.char)+")"			

	def __repr__(self):  # Over-ride string conversion used by print
		return "FiniteField("+str(self.char)+", "+str(self.modpoly)+")"

	def polyprint(self,coeffs,var):  # Coefficients and polynomial variable, e.g. [1,2,2,0,3], "x" yields "1 + 2x + 2x^2 + 3x^4"
		thestr = ""
		firstnon = -1
		for i in range(len(coeffs)): # Find first non-zero coeff
			if(coeffs[i] != 0):
				firstnon = i
				break
		if(firstnon == -1): return "0" # All coeffs are zero
		for i in range(firstnon,len(coeffs)):
			if coeffs[i] == 0: continue
			if(i == 0):  # First non-zero coeff is the constant term 
				thestr = "{0}".format(coeffs[0])
			else: # First non-zero coeff is not the constant term
				if (coeffs[i] < 0): # Negative coeff
					if(i == firstnon): # Different spacing for first
						thestr += "-"
					else:
						thestr += " - "
				elif(i != firstnon):  # Positive coeff (not the first non-zero coeff)
					thestr += " + "
				if(abs(coeffs[i]) != 1): # Suppress printing a coeff of '1'
					thestr += str(abs(coeffs[i]))
				if(i == 1): # x^1 is just 'x'
					thestr += var
				else:
					thestr += var+"^"+str(i)
		return thestr

	def __call__(self,elts=0):  # Coerce constant or array of coeffs as elt of field
		return FiniteFieldElt(self,elts)

	def __iter__(self):
		"""Generator producing all elements of the finite field."""
		digits = [0]*self.degree
		while True:
			yield FiniteFieldElt(self,digits)
			digits[0] += 1
			i = 0
			while digits[i] >= self.char:
				if ((i+1)>= self.degree): raise StopIteration
				digits[i] = 0
				digits[i+1] += 1
				i += 1

	def random(self):
		"""A random element of the finite field."""
		therand = random.randint(0,(self.char)**(self.degree)-1)
		return FiniteFieldElt(self,[(therand//(self.char)**i)%(self.char) for i in range(self.degree)])

	def gen(self):  # Generator element of the field
		"""Generating element of the finite field."""
		return FiniteFieldElt(self,[0,1])

class FiniteFieldElt(object):
	"""An element of a prime power order finite fields"
	Driving polynomial must be monic and top coeff (i.e. 1) is implicit.
	Usage:
	    >>> from finitefield import *
	    >>> GF9 = FiniteField(3,[2,1])    # Define GF(3^2), polys w/ GF3 coeffs, mod x^2+x+2
	    >>> a = FiniteFieldElt(GF9,[1,2]) # Define 2x+1 in GF(9)
	    >>> a**12                         # Compute (2x+1)**12 in GF(9)"""

	def __init__(self, field, elts=0):
		self.field = field
		if (type(elts) == type(0)): # Allow simplified form
			self.coeffs = [elts] + [0 for i in range(self.field.degree-1)]
		else:
			self.coeffs = elts + [0 for i in range(self.field.degree - len(elts))]

	def verbstr(self): # Requires feature from python 2.5.2 or better
		return "("+self.field.polyprint(self.coeffs,self.field.var)+")"

	def __format__(self,fmtspec):  # Over-ride format conversion
		if(fmtspec == ''): fmtspec = self.field.fmtspec
		if(fmtspec == 'p'): # Polynomial format
			return "("+self.field.polyprint(self.coeffs,self.field.var)+")"
		elif(fmtspec == 'c'): # Coeffs only format
			return ''.join([str(self.coeffs[i]) for i in range(len(self.coeffs))])
		else: raise ValueError("***** Error *****: FiniteFieldElt has valid fmtspec values p and c, not <{0}>".format(fmtspec))

	def __str__(self):
		"""over-ride string conversion used by print"""
		return str(self.coeffs)

	def __repr__(self):
		"""over-ride string conversion used by print"""
		return "FiniteFieldElt("+self.field.__repr__()+","+str(self.coeffs)+")"

	def __cmp__(self,other):
		"""compare two elements for equality and allow sorting
		overloaded to allow comparisons to integers and lists of integers"""
		if((type(other) == types.IntType) or (type(other) == types.LongType)):
			return cmp(self.coeffs,[other]+[0 for i in range(self.field.degree-1)])
		elif(type(other) == types.ListType):
			return cmp(self.coeffs,other+[0 for i in range(self.field.degree-len(other))])
		elif(self.field != other.field):
			return -1
		else:
			return cmp(self.coeffs,other.coeffs)

	def norm(self):
		"""The norm of an element over the base field is the product of its conjugates, 
		i.e. norm(a) = a * (a**p) * (a**(p**2)) * ... * (a**(p**k))
		where p is the characteristic of the field and k is the extension degree."""
		return reduce(mul,[self**(self.field.char**k) for k in range(self.field.degree)],1)

	def trace(self):
		"""The trace of an element over the base field is the sum of its conjugates, 
		i.e. trace(a) = a + (a**p) + (a**(p**2)) + ... + (a**(p**k))
		where p is the characteristic of the field and k is the extension degree."""
		return sum([self**(self.field.char**k) for k in range(self.field.degree)])

	def add(self,summand):
		"""add elements of finite fields (overloaded to allow adding integers and lists of integers)"""
		return FiniteFieldElt(self.field, map(lambda x,y: (x+y)%self.field.char, self.coeffs, summand.coeffs))

	def __add__(self,summand):   # Overload the "+" operator
		if ((type(summand) == types.IntType) or (type(summand) == types.LongType)):
			# Coerce if adding integer and FiniteFieldElt
			return self.add(FiniteFieldElt(self.field,[summand]))
		elif(type(summand) == types.ListType):
			return self.add(FiniteFieldElt(self.field,summand))
		else:
			return self.add(summand)

	def __radd__(self,summand):  # Overload the "+" operator when first addend can be coerced to finfld
		if ((type(summand) == types.IntType) or (type(summand) == types.LongType)): # Coerce if adding int and finfld
			return self.add(FiniteFieldElt(self.field,[summand]))
		elif(type(summand) == types.ListType): # Coerce if adding list and finfld
			return self.add(FiniteFieldElt(self.field,summand))
		else:
			return self.add(summand)

	def __iadd__(self,summand): # Overload the "+=" operator
		self = self + summand
		return self

	def __neg__(self):  # Overload the "-" unary operator 
		return FiniteFieldElt(self.field, map(lambda x: self.field.char-x, self.coeffs))

	def __sub__(self,summand):  # Overload the "-" binary operator 
		return self.__add__(-summand)

	def __isub__(self,summand): # Overload the "-=" operator
		self = self - summand
		return self

# 	def mult(self,multand):  # Elementary multiplication in finite fields
# 		thelist = [0 for i in range(self.field.degree)]
# 		for d in range(2*self.field.degree-2):
# 			thelist = map(add, thelist, [sum(self.coeffs[j]*multand.coeffs[d-j] for j in range(max(0,d-self.field.degree),min(d,self.field.degree-1)+1))*i for i in self.field.reduc_table[d]])
# 		return FiniteFieldElt(self.field,map(lambda x: x%self.field.char, thelist))

	def mult(self,multand):  # Elementary multiplication in finite fields
		"""multiply elements of finite fields (overloaded to allow integers and lists of integers)"""
		thelist = [0 for i in range(self.field.degree)]
		for d in range(2*self.field.degree-1):
			for j in range(max(0,d-(self.field.degree-1)),min(d+1,self.field.degree)):
				list2 = [(self.coeffs[j]*multand.coeffs[d-j])*i for i in self.field.reduc_table[d]]
				thelist = map(add, thelist, list2)
		return FiniteFieldElt(self.field,map(lambda x: x%self.field.char, thelist))

	def __mul__(self,multip):  # Overload the "*" operator
		if ((type(multip) == types.IntType) or (type(multip) == types.LongType)): # Coerce if multiply int and finfld
			return self.mult(FiniteFieldElt(self.field,[multip]))
		elif (type(multip) == types.ListType): # Coerce if multiply list and finfld
			return self.mult(FiniteFieldElt(self.field,multip))
		else:
			return self.mult(multip)

	def __rmul__(self,multip):  # Overload the "*" operator
		if ((type(multip) == types.IntType) or (type(multip) == types.LongType)): # Coerce if mult int and and finfld
			return self.mult(FiniteFieldElt(self.field,[multip]))
		elif (type(multip) == types.ListType): # Coerce if mult list and and finfld
			return self.mult(FiniteFieldElt(self.field,multip))
		return self.mult(multip)

	def __imul__(self,multip): # Overload the "*=" operator
		self = self * multip
		return self

	def inv(self):
		"""inverse of element in a finite field"""
		# A better implementation would be xgcd over polynomials
		return self.pow(self.field.order-2)

	def div(self,divisor):
		"""divide elements of a finite field"""
		return self * divisor.inv()

	def __div__(self,divisor):
		if ((type(divisor) == types.IntType) or (type(divisor) == types.LongType)):
			divisor = FiniteFieldElt(self.field,[divisor])  # Coerce if dividing integer and FiniteFieldElt
		return self * divisor.inv()

	def __rdiv__(self,dividend):
		if ((type(dividend) == types.IntType) or (type(dividend) == types.LongType)):
			dividend = FiniteFieldElt(self.field,[dividend])  # Coerce if dividing integer and FiniteFieldElt
		return dividend * self.inv()

	def pow(self,exponent):
		"""pow(b,e) computes the eth power of finite field element b."""
		exponent = (exponent % (self.field.order-1)) # Handle negative and large exponents
		accum = FiniteFieldElt(self.field,[1])
		i = 0
		bpow2 = self
		while ((exponent>>i)>0):
			if((exponent>>i) & 1):
				accum = (accum*bpow2)
			bpow2 = (bpow2*bpow2)
			i+=1
		return accum

	def __pow__(self,exponent): # Overload the "**" operator
		return self.pow(exponent)

	def is_primitive(self):
		if (self**(self.field.order-1) != 1): return False # Not a unit
		for [theprime,thepow] in self.field.facts_order_gpunits:
			if (self**((self.field.order-1)//theprime) == 1): return False
		return True

	def order(self):
		if (self**(self.field.order-1) != 1): 
			raise ValueError('{0} is not a unit in GF({1}^{2})'.format(self,self.field.char,self.field.degree))
		orderaccum = 1
		for [theprime,maxpow] in self.field.facts_order_gpunits:
			theval = self**((self.field.order-1)//(theprime**maxpow))
			for thepow in range(maxpow+1):
				if (theval == 1):
					orderaccum *= (theprime**thepow)
					break
				theval = theval**theprime
		return orderaccum

# Read a file of Conway Polynomials, formatted as per Frank Luebeck
# [http://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/CPimport.txt]
def readconway(filepath,p,e):
        import sys
	try:
		cpfile = open(filepath)
	except IOError as e:
	    print("I/O error({0}): {1}\n   Probably couldn't find file 'CPimport.txt' of Conway Polynomials\n   Try [http://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/CPimport.txt]".format(e.errno, e.strerror))
            pass # Pass the exception up the calling stack
	cpfile.readline()   # Discard the first line
	for theline in cpfile:
		if (theline[0] == '0'):
			cpfile.close()
			raise ValueError('GF({0}^{1}) is beyond the table of Conway Polynomials in {2}'.format(p,e,filepath))
		polyspec = eval(theline[:-2]) # Strip off comma and line end char
		theprime = polyspec[0]; thepow = polyspec[1]
		if (p == polyspec[0] and e == polyspec[1]):
			cpfile.close()
			return polyspec[2]

def findprimpoly(p,e):
	import sys
	if sys.version_info[0] > 2: 
		iterrange = range  # Version 3 patch
	else:
		iterrange = xrange  # Version 2 patch
	ordfacts = numbthy.factor((p**e)-1)
	print (p**e)-1, " = ", ordfacts
	for thepolynum in iterrange(p+1,p**e):
		# Note: Skip [0,d1,d2,...] as not irred
		# Note: Skip [d0,0,0,...,0,1] as not primitive and probably not irred
		thepoly = [(thepolynum//(p**i))%p for i in range(e)]  # length e list of p digits
		a = FiniteField(p,thepoly).gen()  # create poly 'x' mod thepoly, check its order
		print thepoly, str(a)
		if ((a**((p**e)-1)) == 1):    # check for irreducibility, else try next poly
			isprim = True
			for (thefact,thepow) in ordfacts:
				print thefact
				if ((a**(((p**e)-1)/thefact)) == 1): # thepoly not prim
					isprim = False
					break # Don't need to check any more cofactors
			if isprim: return thepoly
	print "Oops" # Should throw exception - should never reach this point, all prime/exp have primitive polys

def GF(n,name='x',modulus=[]):
	if numbthy.is_prime(n):
		if len(modulus)<2:
			return FiniteField(n,[1],name)
		else:
			return FiniteField(n,modulus,name) # Explicit characteristic and polynomial modulus
	else:  # n is not prime - hope it's a prime power
		nfactors = numbthy.factor(n)
		if (len(nfactors) != 1): raise ValueError('GF({0}) only makes sense for {0} a prime power'.format(n))
		p = nfactors[0][0]; e = nfactors[0][1]  # Prime p, exponent e
		if (len(modulus) == 0):
			try:
				modulus = readconway('CPimport.txt',p,e)[:-1] # Assume monic
			except(IOError,ValueError):
				print("      Look for non-Conway primitive polynomial")
				modulus = findprimpoly(p,e)
		return FiniteField(p,modulus,name)

##################################################################################
# More Examples:
# (ref: Primitive Polynomials over Finite Fields, T. Hansen & G. L. Mullins, Math Comp, Oct 1992)
# GF32 = FiniteField(2,[1,0,1,0,0])        # Define GF(2^5) = Z_2[x]/<x^5 + x^2 + 1>
# GF256 = FiniteField(2,[1,0,1,1,1,0,0,0]) # Define GF(2^8) = Z_2[x]/<x^8 + x^4 + x^3 + x^2 + 1>
# GF27 = FiniteField(3,[1,2,0])            # Define GF(3^3) = Z_3[x]/<x^3 + 2x + 2>
# GF6561 = FiniteField(3,[2,0,0,1,0,0,0,0])# Define GF(3^8) = Z_3[x]/<x^8 + x^3 + 2>
# GF25 = FiniteField(5,[2,1])              # Define GF(5^2) = Z_5[x]/<x^2 + x + 2>
# GF125 = FiniteField(5,[2,3,0])           # Define GF(5^3) = Z_5[x]/<x^3 + 3x + 2>
# GF2197 = FiniteField(13,[6,1,0])         # Define GF(13^3) = Z_13[x]/<x^3 + x + 6>
# Irreducible Polynomials:
#    [http://theory.cs.uvic.ca/gen/poly.html]
#    [http://en.wikipedia.org/wiki/Conway_polynomial_(finite_fields)]
# Factorization of Group Orders: (Cunningham Tables)
#    [http://mersennewiki.org/index.php/Cunningham_Tables]
##################################################################################

############################################################################
# License: Freely available for use, abuse and modification
# (this is the Simplified BSD License, aka FreeBSD license)
# Copyright 2001-2017 Robert Campbell. All rights reserved.
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

# 1 June 2014: ver 0.6
#  norm, trace
#  minpoly (difficult - requires linear algebra)
# 19 July 2014: ver 0.7
#  allow direct input of factored group order
#  changed from factors to factor (changed in numbthy.py)
# 21 July 2014: ver 0.71
#  is_primitive()
#  order()
# 18 Oct 2014: ver 0.8
#   [http://www.sagemath.org/doc/reference/rings_standard/sage/rings/finite_rings/constructor.html]
#   [http://www.sagemath.org/doc/reference/rings_standard/sage/rings/finite_rings/element_base.html]
#   [http://www.sagemath.org/doc/reference/rings_standard/sage/rings/finite_rings/integer_mod.html]
# 13 Dec 2014: ver 0.9
#   GF() constructor
#   Conway polynomials
#   (still need to address simple case of prime field)
#   Bug: 1/eltfinfld - fixed (inv and __div__ each referenced the other - need polynomial xgcd)
# 15 Dec 2014: ver 0.91
#   Conway polynomials - fixed bug in readconway routine
#   findprimpoly() - added for brute force primitive polynomial search
# 7 Feb 2017: ver 0.92
#   fix to findprimpoly()
# 16 Sept 2017: ver 0.93
#   Prime field
#   Fixed bug in order() - always returned 1
#   Need to improve formatted output
