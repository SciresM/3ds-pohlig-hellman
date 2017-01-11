from xgcd import xgcd  # http://anh.cs.luc.edu/331/code/xgcd.py
from mod import Mod # http://anh.cs.luc.edu/331/code/mod.py

def ChineseRemainder(pairs):
    '''Return the solution to the Chinese Remainder Theorem, (x, M)
    Pairs contains tuples (a, m) with all m's positive and coprime.
    Return the smallest nonnegative integer x so 
    x mod m  = a mod m for each (a, m) in pairs.
    M is the product ofthe m's.

    >>> pairs = [(2, 3), (3, 4), (1, 5)]
    >>> x, M = ChineseRemainder(pairs)
    >>> (x, M)
    (11, 60)
    >>> for (a, m) in pairs:
    ...     print a % m, x % m
    2 2
    3 3
    1 1
    '''
    (a, m)=list(pairs)[0]
    for (b,p) in list(pairs)[1:]:
        k=((b-a)*xgcd(m,p)[1]) % p #moduli coprime so inverse exists for m mod p
        a=(a+m*k) % (m*p)# joining a mod m and b mod p gives a mod(mp)
        m *= p # mod mp
    return (a,m)

def CRTM(mods):
    '''Return the solution to the Chinese Remainder Theorem as an element
    of a Mod class, x mod N.  The parameter is a nonempty sequence of
    Mod class instances, where all the parameter moduli are coprime.
    N is the product of the parameter moduli.
    For each element a mod m of mods, x satisfies: a mod m equals x mod m.

    >>> vals = [Mod(2, 3), Mod(3, 4), Mod(1, 5)]
    >>> x = CRTM(vals)
    >>> print x
    11 mod 60
    >>> for val in vals:
    ...     print val, int(val - x.value)
    2 mod 3 0
    3 mod 4 0
    1 mod 5 0
    '''
    N = 1 # will be product of moduli
    for a in mods:
        N *= a.modulus()
    x = Mod(0, N)
    for a in mods:
        m = a.modulus()
        Nd = N//m
        gcd, s, t = xgcd(m, Nd) # t*Nd = 1 mod m
        assert gcd == 1, 'moduli not coprime!' # not required, but good idea
        x += a.value*t*Nd
    return x  

if __name__ == '__main__': 
    import doctest
    doctest.testmod() #verbose=True) 

