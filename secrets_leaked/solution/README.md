# Secrets Leaked – Crypto Challenge Writeup

## Challenge Description

We are given an RSA modulus `N = p * q` and a ciphertext `c = m^e mod N`.  
However, the prime `p` is **partially leaked**: its binary expansion is revealed with certain windows of bits set to `0`. The goal is to reconstruct `p` from this partial leakage and recover the flag.

The challenge code snippet:

```python
Kbits1 = 133
Kbits2 = 480

p = Integer(getPrime(512))
q = Integer(getPrime(512))
N = p*q
e = 0x10001

binary_p = p.digits(2)

leaky_p = 21 * [0] + binary_p[21:Kbits1] + [0] * 21 + binary_p[Kbits1+21:Kbits2]

c = pow(bytes_to_long(FLAG), e, N)
print(f"(c, e, N) = {(c, e, N)}")
print(f"leaky = {leaky_p}")
```

So, we are given `(c, e, N)` and the list `leaky`, which encodes the partial bits of `p`.

---

## Solution

The idea is to **model the unknown prime `p` as a polynomial equation** over integers, where the missing bit windows are replaced by variables.

We split the unknown parts into three variables:
- `x` for the first 21 missing bits.
- `y` for the next 21 missing bits.
- `z` for the high unknown bits after position `Kbits2`.

This gives a polynomial form:

```python
f = x \
  + (2**21)*sum([2**i * leaky_p[21+i] for i in range(Kbits1-21)]) \
  + (2**Kbits1)*y \
  + (2**(Kbits1+21))*sum([2**i * leaky_p[21+Kbits1+i] for i in range(Kbits2-(Kbits1+21))]) \
  + (2**Kbits2)*z
```

We then set **bounds** on the unknowns:

```python
bounds = {
    x: (0, 2 ** 21),
    y: (0, 2 ** 21),
    z: (0, 2 ** (512 - Kbits2))
}
```

Using a multivariate coppersmith solver [cuso](https://github.com/keeganryan/cuso), we solve for small roots of the polynomial modulo `p`, where `p` divides `N`. Once `p` is found, we compute `q = N // p`, recover φ(N), and decrypt the ciphertext.

---

## Exploit Script

```python
from Crypto.Util.number import *
import cuso
from sage.all import *

Kbits1 = 133
Kbits2 = 480

(c, e, N) = (13634419770374889083256272578186821114534423923123739989359901337685459195401469811781758106566788842899609154818337907713672837955224864665666600168860751083573289898084466328712427340609265852749807484660227911442835876797469759660424439569815500008759188280075511419173567781621133252851220626364585911845, 65537, 113009588028045785374048610180077310313724894042341422294368188308254605442053307016981271704890877088620963398673272368260369447447767333715517412114888012997924126230159326655850821614635593322309392906842170642307794295461458966535454442267741246125047367925648426989735003561004432987080148623952244343043)
leaky = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0]


x,y,z = var('x','y','z')

f = x +\
 (2**21)*sum([2**i * leaky[21+i] for i in range(Kbits1-21)]) +\
 (2**Kbits1)*y +\
 (2**(Kbits1+21))*sum([2**i * leaky[21+Kbits1+i] for i in range(Kbits2-(Kbits1+21))]) + \
 (2**Kbits2)*z

relations = [f]
bounds = {
    x: (0, 2 ** 21),
    y: (0, 2 ** 21),
    z: (0, 2 ** (512 - Kbits2))
}

roots = cuso.find_small_roots(
    relations,
    bounds,
    modulus="p",
    modulus_multiple=N,
    modulus_lower_bound=2**511,
    modulus_upper_bound=2**512-1
)

p = roots[0]["p"]
q = N // p
d = inverse_mod(e, (p-1)*(q-1))
m = pow(c, d, N)
print(long_to_bytes(m))
```
---

## Flag

```
shellmates{w3ll_H3rm4nN_M4y_is_7h3_w4y_t0_s0lve_th3_ch4ll3ng3!_6B27B4CB8A348ED8750D2C85D9E4B367B90C77046E329F813F2A9CCC330D3E40}
```


---

## References : 

https://eprint.iacr.org/2020/1506.pdf
