import sys
from os.path import basename, splitext
from crt import ChineseRemainder # http://anh.cs.luc.edu/331/code/crt.py

def getPrivExponent(signature):
    ''' Finds DiscreteLog(7, sig, X), where X is a predetermined modulus.
        Uses Pohlig-Hellman to do so.
    '''
    primes = [2,3,5]
    pows = [1960,13,29]
    total = sum(pows)
    
    primepows = [pow(p, n) for p, n in zip(primes, pows)]
    
    # precomputed primitive root/modular inverse for our modulus.
    # Sorry for the magic number :(
    root = 7
    invroot = 0x8c5e4273435fc1f5b3ad2d2492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492492493

    # compute the modulus on the fly to minimize magic numbers.
    modulus = 1
    for x in primepows:
        modulus *= x
    modulus += 1
    
    smooth = modulus - 1 # order of the modulus
    
    invlookup = [pow(invroot, c, modulus) for c in range(max(primes))]
    
    coefficients = [0] * len(primes) # initialize all coefficients to 0
    cnt = 0
    for i in range(len(primes)):
        prime = primes[i]
        val = signature # val starts at signature
        base = smooth # base starts at smooth
        mult = 1
        lookup = [pow(root, j * smooth / prime, modulus) for j in range(prime)] # create lookup table for modular reductions
        for n in range(1,pows[i]+1):
            base /= prime # base = smooth / pow(prime, n), this is faster.
            c = lookup.index(pow(val, base, modulus)) # find coefficient for pow(prime, n) via lookup table
            val = val * pow(invlookup[c], mult, modulus) % modulus # reduce signature by relevant amount
            coefficients[i] += c * mult # coefficient is Cn * pow(prime, n).
            mult *= prime # update multiplier to avoid annoying exponentiation every loop
            cnt += 1
            if (cnt % (total / 100) == 0): # update progress if relevant
                sys.stdout.write('%d%%\r' % (cnt / (total / 100)))
                sys.stdout.flush()
    
    print '\n'
    print 'all done!'
    print '---'

    congruenceList = zip(coefficients, primepows)
    (x,m) = ChineseRemainder(congruenceList)
    
    if x < (pow(0x100, 0x100) - smooth): # Return both possible Ds if more than one is possible
        return [x, x + smooth]
    return [x]
    
def save_bin(fn, dat):
    try:
        with open(fn, 'wb') as f:
            f.write(dat)
        print 'Wrote data to %s!' % fn
    except IOError:
        print 'Failed to write data to %s!' % fn

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: python %s [signature.bin]' % sys.argv[0]
    
    try:
        with open(sys.argv[1], 'rb') as f:
            sig = f.read()
    except IOError:
        print >> sys.stderr, 'Failed to read signature from %s!' % sys.argv[1]
        sys.exit(-1)
    if len(sig) != 0x100:
        print >> sys.stderr, 'Error: signature must be 0x100 bytes for RSA-2048!'
        sys.exit(-1)
    
    sig_num = int(sig.encode('hex'), 16)
    print 'Calculating discrete log for %X...' % sig_num
    privs = getPrivExponent(sig_num)
    privs_hex = ['%0256X' % x for x in privs]
    print 'Private Exponent: '
    print ', or '.join(privs_hex)
    
    bn = splitext(basename(sys.argv[1]))[0]
    if len(privs) == 1:
        save_bin('%s_exp.bin' % bn, privs_hex[0].decode('hex'))
    else:
        for i in range(len(privs)):
            save_bin('%s_exp_%d.bin' % (bn, i), privs_hex[i].decode('hex'))
    
    print 'All done!'
    