################################################################################
# ECDSA Signature
# Python implementation (needs Python EllipticCurve module)
# Author: Robert Campbell, <r.campbel.256@gmail.com>
# Date: 1 June, 2014
# Version 0.3
# License: Simplified BSD (see details at bottom)
######################################################################################
# Refs:
#       NIST FIPS 186-4, March 2013 (and earlier versions to 186-2)
#       ECDSA Test Vectors [http://csrc.nist.gov/groups/ST/toolkit/documents/Examples/ECDSA_Prime.pdf]
######################################################################################

import hashlib, binascii, random

class ECDSA(object):
	"""ECDSA Signatures
	Usage:
	        >>> from ellipticcurve import *
		>>> from ellipticcurve import *
		>>> theECDSA = ECDSA(256,d=0xDEADBEEF)  # d = 0xDEADBEEF private key
		      Q = d*G = [0xB487...8394, 0x2A12...CE5E]
		>>> themsg = "Blah, blah"
		>>> [r,s] = theECDSA.sign(themsg,k=0xF00F00D00F)  # Explicit signing nonce k
		      (r,s) = [0xFF99...4ECB1, 0xDA3F...FB81]
		>>> theECDSA.verify(themsg,[r,s])
	"""
	
	def __init__(self,curve=384,d=None,verbose=True):  # Generate signing key
		self.verbose = verbose
		if (curve == 384):
			[E,G,n,thehash,numbits] = self.define384()
		elif (curve == 256):
			[E,G,n,thehash,numbits] = self.define256()
		elif (curve == 192):
			[E,G,n,thehash,numbits] = self.define192()
		else:
			raise ValueError("*****  Error  *****: curve parameter must be 192, 256 or 384, not {0}".format(curve))
		self.numbits = numbits  # For verbose output formatting
		self.E = E  # The elliptic curve - E:y^2 = x^3 - 3x + b
		# SAGE: self.p = E.base().characteristic()  # Mod p (SAGE specific)
		# SAGE: self.b = self.p - E.defining_polynomial().coefficients()[-1]  # (SAGE specific)
		# SAGE: self.n = E.cardinality()  # Order of the cyclic elliptic curve group (a prime)
		self.p = E.prime  # Mod p (Python class)
		self.b = self.p - E.b  # (Python class specific)
		self.n = E.order  # Order of the cyclic elliptic curve group (a prime) (Python class specific)
		self.G = G  # Basepoint on E with order n
		self.thehash = thehash # Hash function used
		if (d is None):  # Generate a random secret key if none is provided
			d = random.randint(1,n)
		self.d = d
		if (verbose): print "   Private Key: d = 0x{0:0{1}X}".format(self.d,numbits/4)
		self.Q = (self.d)*(self.G)
		if (verbose): print "   Public Key: Q = d*G = [0x{0:0{1}X}, 0x{2:0{3}X}]".format(
			Integer(self.Q[0]),numbits/4,Integer(self.Q[1]),numbits/4)
		
	def sign(self,msg,k=None):
		"""Generate an ECDSA signature for the message provided.
		Usage:
			sage: themsg = "Example of ECDSA with P-384"
			sage: [r,s] = theECDSA.sign(themsg)  # ECDSA signature
			sage: [r,s] = theECDSA.sign(themsg,r=0xA9876543210)  # ECDSA signature w/ explicit nonce
		"""
		verbose = self.verbose
		if (k is None):
			k = random.randint(1,self.n)
		if (verbose): print "   k = 0x{0:0{1}X}".format(k,self.numbits/4)
		r = int((k*(self.G))[0])  # x-coord of point k*G
		if (verbose): print "   r = (k*G)_x = 0x{0:0{1}X}".format(r,self.numbits/4)
		z = int(self.thehash(msg).hexdigest(),16)
		if (verbose): print "   z = hash(msg) = 0x{0:0{1}X}".format(z,self.numbits/4)
		kinv = int(inverse_mod(k,self.n))
		if (verbose): print "   kinv = 1/k (mod n) = 0x{0:0{1}X}".format(kinv,self.numbits/4)
		s = int(mod((z+r*self.d)*kinv,self.n))
		if (verbose): print "   s = kinv*(z+r*d) (mod n) = 0x{0:0{1}X}".format(s,self.numbits/4)
		if (verbose): print " Signature: (r,s) = [0x{0:0{1}X}, 0x{2:0{3}X}]".format(r,self.numbits/4,s,self.numbits/4)
		return [r,s]
		
	def verify(self,msg,thesign):
		"""Verify an ECDSA signature for the message.
		Usage:
			sage: themsg = ""
			sage: theECDSA.verify(themsg,[r,s])
		"""
		verbose = self.verbose
		[r,s] = thesign
		if (verbose): print "   Verify Signature: (r,s) = [0x{0:0{1}X}, 0x{2:0{3}X}]".format(r,self.numbits/4,s,self.numbits/4)
		z = int(self.thehash(msg).hexdigest(),16)
		if (verbose): print "   z = hash(msg) = 0x{0:0{1}X}".format(z,self.numbits/4)
		w = int(inverse_mod(s,self.n))
		if (verbose): print "   w = 1/s (mod n) = 0x{0:0{1}X}".format(w,self.numbits/4)
		u1 = int(mod(z*w,self.n))
		if (verbose): print "   u1 = z*w (mod n) = 0x{0:0{1}X}".format(u1,self.numbits/4)
		u2 = int(mod(r*w,self.n))
		if (verbose): print "   u2 = r*w (mod n) = 0x{0:0{1}X}".format(u2,self.numbits/4)
		# SAGE: rprime = int(mod((u1*(self.G)+u2*(self.Q))[0],self.n)) # x-coord of point (u1*G+u2*Q)		
		rprime = int(mod((u1*(self.G)+u2*(self.Q)).x,self.n)) # x-coord of point (u1*G+u2*Q)
		if (verbose): print "   r' = (u1*G + u2*Q)_x = 0x{0:0{1}X}".format(rprime,self.numbits/4)
		if (r == rprime):
			if (verbose): print "   Signature verifies"
			return True
		else:
			if (verbose): print "   Signature FAILS to verify"
			return False

	##############################################################################################  
	# NIST P-192/256/384 Elliptic Curves
	# [http://csrc.nist.gov/groups/ST/toolkit/documents/dss/NISTReCur.pdf]
	# [http://www.secg.org/collateral/sec2_final.pdf]
	# FIPS 186-2, pg 29 [http://csrc.nist.gov/publications/fips/archive/fips186-2/fips186-2.pdf]
	# FIPS 186-4, Sect D.1.2.1, pg 90 [http://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-4.pdf]
	##############################################################################################
	def define192(self):
		p = 2**192 - 2**64 - 1 # 6277101735386680763835789423207666416083908700390324961279
		b = 0x64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1
		n = 6277101735386680763835789423176059013767194773182842284081
		E = EllipticCurve(GF(p),[-3,b])
		#SAGE: E.set_order(n)  # Set the pre-computed curve order
		E.order = n  # Set the pre-computed curve order
		gx = 0x188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012
		gy = 0x07192b95ffc8da78631011ed6b24cdd573f977a11e794811
		g = E.point([gx,gy])  # The basepoint
		return [E,g,n,hashlib.sha1,192]
		
	def define256(self):
		p = 2**256 - 2**224 + 2**192 + 2**96 - 1 # 115792089210356248762697446949407573530086143415290314195533631308867097853951
		b = 41058363725152142129326129780047268409114441015993725554835256314039467401291
		n = 115792089210356248762697446949407573529996955224135760342422259061068512044369
		E = EllipticCurve(GF(p),[-3,b])
		E.set_order(n)  # Set the pre-computed curve order
		gx = 48439561293906451759052585252797914202762949526041747995844080717082404635286
		gy = 36134250956749795798585127919587881956611106672985015071877198253568414405109
		g = E.point([gx,gy])  # The basepoint
		return [E,g,n,hashlib.sha256,256]

	def define384(self):
		p = 2**384 - 2**128 - 2**96 + 2**32 - 1 # 3940...2319
		b = 27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575
		n = 39402006196394479212279040100143613805079739270465446667946905279627659399113263569398956308152294913554433653942643
		E = EllipticCurve(GF(p),[-3,b])
		E.set_order(n)  # Set the pre-computed curve order
		gx = 26247035095799689268623156744566981891852923491109213387815615900925518854738050089022388053975719786650872476732087
		gy = 8325710961489029985546751289520108179287853048861315594709205902480503199884419224438643760392947333078086511627871
		g = E.point([gx,gy])  # The basepoint
		return [E,g,n,hashlib.sha384,384]

############################################################################
# License: Freely available for use, abuse and modification
# (this is the Simplified BSD License, aka FreeBSD license)
# Copyright 2014 Robert Campbell. All rights reserved.
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
