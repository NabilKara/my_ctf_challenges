from secret import FLAG
from Crypto.Util.number import getPrime , bytes_to_long

Kbits1 = 133
Kbits2 = 480


p = Integer(getPrime(512))
q = Integer(getPrime(512))
N = p*q
e = 0x10001

binary_p =  p.digits(2)

leaky_p = 21 * [0] + binary_p[21:Kbits1] +  [0] * 21 + binary_p[Kbits1+21:Kbits2]

c = pow(bytes_to_long(FLAG), e, N)
print(f"(c, e, N) = {(c, e, N)}")
print(f"leaky = {leaky_p}")
