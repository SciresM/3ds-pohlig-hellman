# 3ds-pohlig-hellman
Implements an RSA keyslot recovery attack on the 3DS's hardware RSA engine.

Usage: python pohlig_hellman.py [signature.bin]

Once you have your exponents, insert them into keygen.py, and run to generate the 3ds's 6.x/7.x keys.

# Credits:
-Myria, for discovering (and explaining) [the flaw in the hardware engine](http://3dbrew.org/wiki/3DS_System_Flaws) being exploited.

-[Dr. Andrew Harrington](http://anh.cs.luc.edu/), for the Chinese remainder theorem code.
