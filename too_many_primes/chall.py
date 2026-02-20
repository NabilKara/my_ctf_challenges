from Crypto.Util.number import bytes_to_long, getPrime
from secret import FLAG
import random

m = bytes_to_long(FLAG)
n = 1
for i in range(40):
    size = random.randint(42, 70) 
    p = getPrime(size)
    n *= p

e = 65537
c = pow(m, e, n)

print("n =", n)
print("c =", c)
