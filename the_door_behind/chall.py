#!/usr/bin/env python3
from Crypto.Util.number import bytes_to_long, getPrime, isPrime
from Crypto.Random.random import randint
from secret import FLAG

def weird_prime():
    D = 635
    while 1:
        s = randint(pow(2,1020), pow(2,1021) - 1)
        tmp = D * s ** 2 + 1
        if tmp % 4 == 0 and isPrime((tmp // 4)):
            return tmp // 4

p = weird_prime()
q = getPrime(2048)
n = p * q
e = 65537
m = bytes_to_long(FLAG)
c = pow(m, e, n)

print("n : " , n)
print("c : " , c)