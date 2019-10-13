# Number Theory Python
A collection of Python modules implementing number theory stuff.

The Number-Theory-Python package currently includes:
* NumbThy.pdf: A short course in number theory
* numbthy.py: Basic number theory functions
* gaussint.py: Basic operations over the Gaussian integers
* finitefield.py: Finite fields of prime power order
* CPimport.txt: Data file of finite field defining polynomials
* ellipticcurve.py: Elliptic curves in affine reduced Weierstrass form over prime order fields
* ECDSA.py: Elliptic curve signatures

Functions implemented in numbthy.py are:
* gcd(a,b) - Compute the greatest common divisor of a and b.
* xgcd(a,b) - Find [g,x,y] such that g=gcd(a,b) and g = ax + by.
* power_mod(b,e,n) - Compute b^e mod n efficiently.
* inverse_mod(b,n) - Compute 1/b mod n.
* is_prime(n) - Test whether n is prime using a variety of pseudoprime tests.
* euler_criterion(a, p) - Test whether a is a quadratic residue mod p.
* euler_phi(n) - Compute Euler's Phi function of n - the number of integers strictly less than n which are coprime to n.
* carmichael_lambda(n) - Compute Carmichael's Lambda function of n - the smallest exponent e such that b\*\*e = 1 for all b coprime to n.
* factor(n) - Return a sorted list of the prime factors of n with exponents.
* prime_divisors(n) - Returns a sorted list of the prime divisors of n.
* is_primitive_root(g,n) - Test whether g is primitive - generates the group of units mod n.
* sqrtmod(a,n) - Compute sqrt(a) mod n using various algorithms.
* TSRsqrtmod(a,grpord,n) - Compute sqrt(a) mod n using Tonelli-Shanks-RESSOL algorithm.

Classes implemented in finitefield.py are:
* FiniteField(p,polycoeffs) is the finite field of characteristic p and given defining polynomial.  The methods defined for this class are:
    * str(ff), format(ff) and repr(ff) - also implicitly used by print and display functions
    * Coerce integer or array of integers to element
    * Iterator over all elements of the finite field
    * ff.random() - Random element of the finite field
    * If a defining polynomial is not specified, the Conway polynomial is used
* FiniteFieldElt(ff,polycoeffs) is the element in the specified finite field.  The methods defined for this class are:
    * str(elt), format(elt) and repr(elt) - also implicitly used by print and display functions
    * elt1.add(elt2) - Add two elements - also overloads the + operator, so elt1 + elt2 (and elt1 += elt2)
    * elt1.mul(elt2) - Multiply two elements - also overloads the * operator, so elt1 * elt2 (and elt1 *= elt2)
    * elt1.neg() - Additive inverse (aka negative) of an element - also overloads the unary - operator, so -elt1
    * elt1.inv() - Multiplicative inverse of an element
    * elt1.div(elt2) - Divide two elements - also overloads the / operator, so elt1 / elt2 (and elt1 /= elt2)
    * elt1.pow(n) - The nth power of an element - also overloads the ** operator, so elt1**n
    * elt1.is_primitive() - Does elt1 generate the multiplicative group
    * elt1.order() - Compute the multiplicative order of elt1

Classes implemented in ellipticcurve.py are:
* EllipticCurve(p,[a,b]) is the elliptic curve in Weierstrass form y^2 = x^3 +ax + b (mod p) over the finite field of prime order p.  The methods defined for this class are:
    * str(ec), format(ec) and repr(ec) - also implicitly used by print and display functions
* EllipticCurveElt(ec,[x,y]) is a point on the specified elliptic curve.  The methods defined for this class are:
    * str(pt), format(pt) and repr(pt) - also implicitly used by print and display functions
    * pt1.add(pt2) - Add two points - also overloads the + operator, so pt1 + pt2 (and pt1 += pt2)
    * pt1.neg() - Additive inverse (aka negative) of a point - also overloads the unary - operator, so -pt1
    * pt1.mult(n) - The nth multiple of a point - also overloads the * operator, so n*pt1

Classes implemented in ECDSA.py are:
* ECDSA(size,d=key) defines signatures over the NIST curve of size bits, using the optional value d as a private key (otherwise it is randomly generated).  The methods defined for this class are:
    * ecdsa.sign(msg,k=nonce) - Signs the message string, using the optional value k as a signing nonce (othewise it is randomly generated).
    * ecdsa.verify(msg,thesign) - Verifies that thesign is a valid signature for the string msg.

Classes implemented in gaussint.py are:
* GaussInt(a,b) is the Gaussian integer a + bI, where a and b are integers.  The methods defined for this class are:
    * str(ff) and repr(ff) - also implicitly used by print and display functions
    * Coerce integer or complex to element
    * elt1.add(elt2) - Add two elements - also overloads the + operator, so elt1 + elt2 (and elt1 += elt2)
    * elt1.mult(elt2) - Multiply two elements - also overloads the * operator, so elt1 * elt2 (and elt1 *= elt2)
    * elt1.neg() - Additive inverse (aka negative) of an element - also overloads the unary - operator, so -elt1
    * elt1.div(elt2) - Divide two elements - also overloads the / operator, so elt1 / elt2 (and elt1 /= elt2)
    * elt1.mod(elt2) - Reduce one element mod another - also overloads the % operator, so elt1 % elt2 (and elt1 %= elt2)
    * elt1.divmod(elt2) - Divide two elements, returning both quotient and remainder
    * elt1.gcd(elt2) - The gcd of two elements
    * elt1.xgcd(elt2) - The extended gcd of two elements
    * elt1.pow(n) - The nth power of an element - also overloads the ** operator, so elt1**n
    * elt1.isprime() - Tests whether elt1 is prime (as a Gaussian integer)
    * elt1.factors() - Return a list of the prime factors of elt1
    * elt1.factor() - Return a single prime factor of elt1

