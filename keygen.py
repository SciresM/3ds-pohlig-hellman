# Python 3 support
from __future__ import print_function

import binascii
import hashlib
import sys

# These need to be filled in.
# I suggest https://github.com/SciresM/3ds-pohlig-hellman,
# although you can also find them in boot9_prot.
retail_exponent = ''
dev_exponent = ''

retail_modulus = 'C12E4877FF0FEDB98AFA6DBE5D7AC86489A4E0901DBA22D2E7BA945C830126CA71976F284A5E57E386B4778FBB78E3F1750ECD8B6ED3984F3CD363F801DD9ED5271EA859F49FF4DA02D88445EFBA774B3A46279167CB905300D0903ECB71772B78CFE918921E2BEF699F87D9CBEF273FBC0E320D177E925021E869BB827F227F17DA206AE9CE5DE262566C4FAB0BE86866F7218FD48C248C0D8CBDEB3ADBC423A6CD1B572F82C0172D23E32B1B4F0C30B04C20749126851DBA39B31A2918C23F91034969E601C0A1099742E8FD31CFE4A02CD9E6BD7B646C52465567C7463B2A9E7E9FA53A0E260D4DAEE68091EE7ADFA1C1491538276161E1C10293F4592A33'
dev_modulus = 'EB686B6005538368435FAE6D2016521E218058B17D9B67EA1B71042835FA39F86BFD0BDD11445846705A49E5D26117528F4DE4A0D5F186D068DE06D8382FBFCAAE403A61540B6E3A24D9A283D340B6BEC6802DAB3441E4260639D66B134E8123FB06748926CC93089E3FE539DA765A7CC0486EBAEC74C26F2B41A8B118E3E9355104BE57DAA370494E09DAC7611D625639EFDE63E02B71EDB00CBA9509E49A3A8C46F19B22107193FF840088BEB3AD510E720B024A36FBBED81B4AC0A0512ED63E1346FDBD8452C0160F8774FFDBE708C761F49DD9CA138B418D947E98774A87A4CDAAED7021C5AA30F88ACC4B07E531BBDFC71A879AE85EDB8214F6D23910DB'


keydata = 'A20F47820E3EEF0336FCFCE37220CFD49BFE5A93BB8755ACD3043FEF5D507ED42459BE6EA64BED00FD28CCFB25A72F9DC74C7359C130E9DE332A6B8C33D18039'
asn1_prefix = '0001FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF003031300D060960864801650304020105000420'

def generate_keys(d, n):
    ''' Prints out the 6.x and 7.x keys generated by using an exponent/modulus.'''

    # Generate ASN.1 data
    asn1 = asn1_prefix + hashlib.sha256(binascii.unhexlify(keydata)).hexdigest()
    asn1 = int(asn1, 16)

    if type(d) == str:
        d = int(d, 16)
    if type(n) == str:
        n = int(n, 16)

    # Exponentiate it.
    keys = hashlib.sha256(binascii.unhexlify('%0512X' % pow(asn1, d, n))).hexdigest().upper()
    print('6.X 0x2F KeyY: %s' % keys[:0x20])
    print('7.x 0x25 KeyX: %s' % keys[0x20:])

def is_key_valid(d, n):
    e = 0x10001
    if type(d) == str:
        d = int(d, 16)
    if type(n) == str:
        n = int(n, 16)

    test_vector = 0xCAFEBABE
    return pow(pow(test_vector, d, n), e, n) == test_vector

if __name__ == '__main__':
    if not retail_exponent and not dev_exponent:
        print('Neither retail nor dev exponents are set up!')
        print('Please insert them before running the program.')
        sys.exit(0)

    if retail_exponent:
        print('Retail:')
        if is_key_valid(retail_exponent, retail_modulus):
            generate_keys(retail_exponent, retail_modulus)
        else:
            print('Retail exponent is incorrect. Double check it?')
    if dev_exponent:
        print('Dev:')
        if is_key_valid(dev_exponent, dev_modulus):
            generate_keys(dev_exponent, dev_modulus)
        else:
            print('Dev exponent is incorrect. Double check it?')