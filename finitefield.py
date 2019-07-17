######################################################################################
# FINITEFIELD.PY
# A finite field of prime power order
# Note: Version 0.9 - major changes to improve SAGE compatibility
#    Version 0.8 is compatible with both Python 2.5+ and 3.x 
#       and changes some names to improve SAGE compatibility
# Author: Robert Campbell, <r.campbel.256@gmail.com>
# Date: 23 Sept, 2018
# Version 0.971
# License: Simplified BSD (see details at bottom)1
######################################################################################
"""Finite fields.
	Functions implemented are:
		+, *, /, **, norm, trace, gen, random, enumerate
	Simple Examples:
	    >>> from finitefield import *     # names do not need path
	    >>> GF81 = GF(3**4)               # create finite field
	    >>> a = GF81([0,1])               # create element [0,1] = (0 + 1*x) = x
	    >>> '{0:l}^10 has value {1:l}'.format(a,a**10) # format as coeff list
	    >>> GF81.random_element().verbstr()  # a random value, formatted as polynomial
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
	  Output is coefficients only ('c'), not full polynomials ('p')
	    >>> GF5e8 = FiniteField(5,[4,3,1,0,1,3,0,0],'z',[[2,5],[3,1],[13,1],[313,1]],'c')
	    >>> a = FiniteFieldElt(GF5e8,[0,1])
	    >>> format(a**20)  # Default format, set when defining FiniteField
	        '20340110'
	    >>> '{0:p}'.format(a**20)  # Force polynomial format
	        '(2 + 3*z**2 + 4*z**3 + z**5 + z**6)'
	    >>> '{0:l}'.format(a**20)  # Force coeff list format
	        '[2, 0, 3, 4, 0, 1, 1, 0]'
	    >>> '{0:c}'.format(a**20)  # Force packed coeffs format
	        '20340110'
	    >>> print a**20
	        [2, 0, 3, 4, 0, 1, 1, 0]
	    >>> a**20
	        FiniteFieldElt(FiniteField(5, [4, 3, 1, 0, 1, 3, 0, 0]),[2, 0, 3, 4, 0, 1, 1, 0])"""


__version__ = '0.971' # Format specified in Python PEP 396
Version = 'finitefield.py, version ' + __version__ + ', 23 Sept, 2018, by Robert Campbell, <r.campbel.256@gmail.com>'

import numbthy  # Use factor
import random   # Generate random elements
import sys      # Check Python2 or Python3
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
		       '(2 + 3*z**2 + 4*z**3 + z**5 + z**6)'"""

		self.char = prime
		self.degree = len(poly)
		self.order = self.char**self.degree
		self.modpoly = poly
		self.var = var
		self.fmtspec = fmtspec
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
			return "Z_"+str(self.char)+"["+self.var+"]/<"+self.polyprint(self.modpoly+[1],var=self.var)+">"
		else:
			return "Z_"+str(self.char)
			
	def __format__(self,fmtspec):  # Over-ride format conversion
		"""Override the format when outputting a finite field.
		A default can be set when the field is defined or it can be specified for each output.
		Possible formats are:
			p - polynomial format
			l - list of coefficients
			c - coefficients in packed format
			f - full format, can be used as input
			s - short format
			t - LaTeX format
			tl - LaTeX format (long)
			ts - LaTeX format (short)
"""
		if(fmtspec == ''): fmtspec = self.fmtspec
		if(fmtspec == 'p'): # Polynomial format
			return "GF("+str(self.char)+"**"+str(self.degree)+","+self.polyprint(self.modpoly+[1],var=self.var)+")"
		elif(fmtspec == 'l'): # Coefficient list format
			return "GF("+str(self.char)+"**"+str(self.degree)+","+format(self.modpoly+[1])+")"
		elif(fmtspec == 'c'): # Coeffs only format
			return "GF("+str(self.char)+"**"+str(self.degree)+","+''.join([str(self.modpoly[i]) for i in range(len(self.modpoly))])+"1"+")"
		elif(fmtspec == 'f'): # Full form format - can be input
			return "FiniteField("+str(self.char)+","+str(self.modpoly)+")"
		if(fmtspec == 't' or fmtspec == 'tl'): # long LaTeX format
			return "{GF("+str(self.char)+")["+str(self.var)+r']/\left\langle{'+self.polyprint(self.modpoly+[1],var=self.var,fmtspec='t')+r'}\right\rangle}'
		if(fmtspec == 'ts'): # short LaTeX format
			return "{GF("+str(self.char)+"^{"+str(self.degree)+"})}"
		elif(fmtspec == 's'): # Short format - field size only (can be input though)
			return "GF("+str(self.char)+"**"+str(self.degree)+")"
		else: raise ValueError("***** Error *****: FiniteField has valid fmtspec values p, l, c, f, s, t, tl and ts, not <{0}>".format(fmtspec))
		return str(self)  # Get to this later

	def __str__(self):  # Over-ride string conversion used by print
		if(self.degree > 1):
			return "GF("+str(self.char)+"^"+str(self.degree)+")"
		else:
			return "GF("+str(self.char)+")"			

	def __repr__(self):  # Over-ride format conversion
		return '{0:f}'.format(self)

	def polyprint(self,coeffs,var='X',fmtspec='p'):  # Coefficients and polynomial variable, e.g. [1,2,2,0,3], "x" yields "1 + 2x + 2x^2 + 3x^4"
		"""polyprint(coeffs,var=thevar,fmtspec=thefmtspec) prints the coefficients in polynomial form, 
		using the variable thevar. Formats can be Python/SAGE, ie 1 - 2*x + 2*x**2 + 3*x**4, 
		or LaTeX, ie 1 - 2x + 2x^{2} + 3x^{4}
		Possible formats are:
			p - polynomial format
			f - full format, can be used as input
			t - LaTeX format"""
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
					if (fmtspec=='t'):
						thestr += str(abs(coeffs[i]))
					else:
						thestr += (str(abs(coeffs[i])) + "*")
				if(i == 1): # x^1 is just 'x'
					thestr += var
				else:
					if (fmtspec=='t'):
						thestr += var+"^{"+str(i)+"}"
					else:
						thestr += var+"**"+str(i)
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

	def random_element(self):
		"""A random element of the finite field."""
		therand = random.randint(0,(self.char)**(self.degree)-1)
		return FiniteFieldElt(self,[(therand//(self.char)**i)%(self.char) for i in range(self.degree)])

	def random(self):
		"""A random element of the finite field. (Renamed random_element in ver 0.96)"""
		return self.random_element()
                
	def gen(self):  # Generator element of the field
		"""Usual generating element of the finite field."""
		return FiniteFieldElt(self,[0,1])

class FiniteFieldElt(object):
	"""An element of a prime power order finite fields"
	Driving polynomial must be monic and top coeff (i.e. 1) is implicit.
	Usage:
	    >>> from finitefield import *
	    >>> GF9 = FiniteField(3,[2,1])    # Define GF(3^2), polys w/ GF3 coeffs, mod x^2+x+2
	    >>> a = FiniteFieldElt(GF9,[1,2]) # Define 2x+1 in GF(9)
	    >>> a**12                         # Compute (2x+1)**12 in GF(9)"""

	def isIntType(self,x):
		if sys.version_info < (3,): return isinstance(x,(int, long,))
		else: return isinstance(x,(int,))
                
	def __init__(self, field, elts=0):
		self.field = field
		if self.isIntType(elts): # Allow coercion from integer
			self.coeffs = [mod(elts,self.field.char)] + [0 for i in range(self.field.degree-1)]
		else:
			self.coeffs = [mod(theelt,self.field.char) for theelt in elts] + [0 for i in range(self.field.degree - len(elts))]

	def verbstr(self): # Requires feature from python 2.5.2 or better
		return "("+self.field.polyprint(self.coeffs,var=self.field.var)+")"

	def __format__(self,fmtspec):  # Over-ride format conversion
		"""Override the format when outputting a finite field element.
		A default can be set when the field is defined or it can be specified for each output.
		Possible formats are:
			p - polynomial format
			l - list of coefficients
			c - coefficients in packed format
			t - LaTeX format
			f - full format - can be used as input
		Example:
			>>> GF64 = GF(2**6,var='z',fmtspec='c')
			>>> a = GF64([1,0,1,1,0,1])
			>>> format(a)
			'101101'
			>>> '{:p}'.format(a)
			'(1 + z^2 + z^3 + z^5)'
			>>> '{:l}'.format(a)
			'[1, 0, 1, 1, 0, 1]'
			>>> '{:c}'.format(a)
			'101101'  """
		if(fmtspec == ''): fmtspec = self.field.fmtspec
		if(fmtspec == 'p' or fmtspec == 's'): # Polynomial format
			return "("+self.field.polyprint(self.coeffs,var=self.field.var)+")"
		elif(fmtspec == 'l'): # Coefficient list format
			return format(self.coeffs)
		elif(fmtspec == 'c'): # Coeffs only format
			return ''.join([str(self.coeffs[i]) for i in range(len(self.coeffs))])
		elif(fmtspec == 't'): # LaTeX format
			return "("+self.field.polyprint(self.coeffs,var=self.field.var,fmtspec='t')+")"
		elif(fmtspec == 'f'): # Full form format - can be input
			return "FiniteFieldElt("+'{0:f}'.format(self.field)+","+str(self.coeffs)+")"
		else: raise ValueError("***** Error *****: FiniteFieldElt has valid fmtspec values p, l, c and f, not <{0}>".format(fmtspec))

	def __str__(self):
		"""over-ride string conversion used by print"""
		return '{0:p}'.format(self)

	def __repr__(self):
		"""over-ride format conversion"""
		return '{0:f}'.format(self)

	def __cmpfinfld__(self,other): # Implement cmp for both Python2 and Python3
		"""compare two elements for equality and allow sorting
		overloaded to allow comparisons to integers and lists of integers"""
		if self.isIntType(other):  # Coerce if comparing int and and finfld
			return self.listcmp(tuple(reversed(self.coeffs)),tuple(reversed([other]+[0 for i in range(self.field.degree-1)])))
		elif isinstance(other,(list,tuple,)): # Coerce if comparing list (of coeffs) and finfld
			return self.listcmp(tuple(reversed(self.coeffs)),tuple(reversed(other+[0 for i in range(self.field.degree-len(other))])))
		elif(self.field != other.field):
			raise ValueError("Cannot compare elements of different FiniteFields: <{0}> and <{1}>".format(self.field,other.field))
		else:
			return self.listcmp(tuple(reversed(self.coeffs)),tuple(reversed(other.coeffs)))

	def listcmp(self,list1,list2): # Implement list cmp for Python3
		for ptr in range(len(list1)):
			if (list1[ptr] < list2[ptr]): return -1
			if (list1[ptr] > list2[ptr]): return 1
		else: return 0

	def __cmp__(self,other): return self.__cmpfinfld__(other) # Used by Python2 for sorting
	def __lt__(self,other): return (self.__cmpfinfld__(other) < 0)
	def __gt__(self,other): return (self.__cmpfinfld__(other) > 0)
	def __eq__(self,other): return (self.__cmpfinfld__(other) == 0)
	def __le__(self,other): return (self.__cmpfinfld__(other) <= 0)
	def __ge__(self,other): return (self.__cmpfinfld__(other) >= 0)
	def __ne__(self,other): return (self.__cmpfinfld__(other) != 0)

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
		return FiniteFieldElt(self.field, tuple(map(lambda x,y: (x+y)%self.field.char, self.coeffs, summand.coeffs)))

	def __add__(self,summand):   # Overload the "+" operator
		if self.isIntType(summand): # Coerce if adding integer and FiniteFieldElt
			return self.add(FiniteFieldElt(self.field,[summand]))
		elif isinstance(summand,(list,tuple,)): # Coerce if adding list (of coeffs) and FiniteFieldElt
			return self.add(FiniteFieldElt(self.field,summand))
		else:
			return self.add(summand)

	def __radd__(self,summand):  # Overload the "+" operator when first addend can be coerced to finfld
		if self.isIntType(summand): # Coerce if adding int and finfld
			return self.add(FiniteFieldElt(self.field,[summand]))
		elif isinstance(summand,(list,tuple,)): # Coerce if adding list and finfld
			return self.add(FiniteFieldElt(self.field,summand))
		else:
			return self.add(summand)

	def __iadd__(self,summand): # Overload the "+=" operator
		self = self + summand
		return self

	def __neg__(self):  # Overload the "-" unary operator 
		return FiniteFieldElt(self.field, tuple(map(lambda x: self.field.char-x, self.coeffs)))

	def __sub__(self,summand):  # Overload the "-" binary operator 
		return self.__add__(-summand)

	def __isub__(self,summand): # Overload the "-=" operator
		self = self - summand
		return self

	def mult(self,multand):  # Elementary multiplication in finite fields
		"""multiply elements of finite fields (overloaded to allow integers and lists of integers)"""
		thelist = [0 for i in range(self.field.degree)]
		for d in range(2*self.field.degree-1):
			for j in range(max(0,d-(self.field.degree-1)),min(d+1,self.field.degree)):
				list2 = [(self.coeffs[j]*multand.coeffs[d-j])*i for i in self.field.reduc_table[d]]
				thelist = map(add, thelist, list2)
		return FiniteFieldElt(self.field, tuple(map(lambda x: x%self.field.char, thelist)))

	def __mul__(self,multip):  # Overload the "*" operator
		if self.isIntType(multip): # Coerce if multiply int and finfld
			return self.mult(FiniteFieldElt(self.field,[multip]))
		elif isinstance(multip,(list,tuple,)): # Coerce if multiply list and finfld
			return self.mult(FiniteFieldElt(self.field,multip))
		else:
			return self.mult(multip)

	def __rmul__(self,multip):  # Overload the "*" operator
		if self.isIntType(multip): # Coerce if mult int and and finfld
			return self.mult(FiniteFieldElt(self.field,[multip]))
		elif isinstance(multip,(list,tuple,)): # Coerce if mult list and and finfld
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
		if self.isIntType(divisor):
			divisor = FiniteFieldElt(self.field,[divisor])  # Coerce if dividing integer and FiniteFieldElt
		return self * divisor.inv()

	def __rdiv__(self,dividend):
		if self.isIntType(dividend):
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
		"""is_primitive(e) returns True if the finite field element e has full order,
		i.e. the powers of e exhaust all non-zero elements of the field."""
		if (self**(self.field.order-1) != 1): return False # Not a unit
		for [theprime,thepow] in self.field.facts_order_gpunits:
			if (self**((self.field.order-1)//theprime) == 1): return False
		return True

	def multiplicative_order(self):
		"""multiplicative_order(b) returns the smallest positive non-zero exponent e such that b**e == 1."""
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

	def order(self):
		"""order(b) returns the smallest positive non-zero exponent e such that b**e == 1.  (Renamed multiplicative_order in ver 0.96)"""
		return self.multiplicative_order()

	def minimal_polynomial(self):  # Find first element of nullspace in matrix (1, a, a^2, a^3, ..., a^deg)
		"""minimal_polynomial(a) is the least degree monic polynomial which has
		a value of zero when a is substituted.  (Currently returns a list of 
		coefficients - will revisit when a polynomial package has been written.)
		To see in conventional polynomial format use 
		FiniteField.polyprint(a.minimal_polynomial(),var='X')"""
		thedegree = self.field.degree
		themod = self.field.char
		numrows = self.field.degree+1
		numcols = 2*self.field.degree+1
		themat = [(self**i).coeffs + [(0 if (i!=j) else 1) for j in range(numrows)] for i in range(numrows)]
		###### Nested function
		def reducemat(submat):
			pivotrow=0
			for rj in range(thedegree):
				# Find the pivot element in this column
				for ri in range(pivotrow,len(submat)):
					if (submat[ri][rj] != 0):
						# Swap rows
						temp = submat[pivotrow]; submat[pivotrow] = submat[ri]; submat[ri] = temp
						break
				else:
					continue # No pivot in this column - move to next column (for rj)
				# Subtract multiples of pivot row from lower rows (zeroing out column elts)
				invpivot = numbthy.inverse_mod(submat[pivotrow][rj],themod)
				for ri in range(pivotrow+1,len(submat)):
					rowmult = mod(invpivot*submat[ri][rj],themod)
					for rjj in range(rj,numcols):
						submat[ri][rjj] = mod(submat[ri][rjj]-rowmult*submat[pivotrow][rjj],themod)
					#if not any(themat[i][:thedegree]): return themat[i][thedegree:]
				pivotrow += 1  # Move down to next row
		###### end Nested function reducemat
		for i in range(2,numrows+1): # Reduce each submatrix - look for first null
			reducemat(themat[:i])
			if not any(map(lambda x:x!=0, themat[i-1][:thedegree])): return themat[i-1][thedegree:]
			
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
	"""A brute-force search for a primitive polynomomial mod p of degree e.
	Needs to be fixed, as it relies on a bug in the FiniteField code, allowing
	a 'field' to be defined mod a reducible polynomial.  Also needs better check
	for reducibility - often fails for non-prime degrees.  Need to implement Rabin test in polynomial class."""
	raise NotImplementedError('findpripoly() is not implemented.  Previous implementation was buggy and removed.')

def GF(n, poly=[], var='x', fmtspec="p"):
	"""A shorthand for generating finite fields.  If poly is not specified then one will be chosen from a list of Conway polynomials."""
	if numbthy.is_prime(n):
		if len(poly)<2:
			return FiniteField(n,[1],var=var,fmtspec=fmtspec)
		else:
			return FiniteField(n,poly,var=var,fmtspec=fmtspec) # Explicit characteristic and polynomial modulus
	else:  # n is not prime - hope it's a prime power
		nfactors = numbthy.factor(n)
		if (len(nfactors) != 1): raise ValueError('GF({0}) only makes sense for {0} a prime power'.format(n))
		p = nfactors[0][0]; e = nfactors[0][1]  # Prime p, exponent e
		if (len(poly) == 0):
			try:
				poly = readconway('CPimport.txt',p,e)[:-1] # Assume monic
			except(IOError,ValueError):
				print("      Look for non-Conway primitive polynomial")
				poly = findprimpoly(p,e)
		else:
			if nfactors[0][1] != len(poly):
				raise ValueError('Polynomial {0} does not have degree {1}'.format(poly+[1],nfactors[0][1]))
		return FiniteField(p,poly,var=var,fmtspec=fmtspec)

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
# 22 Jan 2018: ver 0.94
#   Formatting - both for FiniteField and FiniteFieldElt
#   Added document strings
# 28 Jan 2018: ver 0.95
#   Fixed support for Python 3
#   Removed buggy findprimpoly (often returned reducible poly if degree composite)
# 10 Feb 2018: ver 0.96
#   Add minimal_polynomial() method
#   Minor renaming for SAGE compatibility (multiplicative_order and random_element)
#   Minor formatting changes
# 11 Mar 2018: ver 0.97
#   Fix various bugs in GF
#   polyprint: Let default var be 'X'
# 11 Mar 2018: ver 0.971
#   Minor format fixes to repr and str

