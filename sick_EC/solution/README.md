# Sick EC

**Author:** 0xn4b1l

**Challenge Overview**

This challenge presents a `ECDH` problem  with scalar multiplication and  addition implemented on a **singular**  elliptic curve over a large prime field `GF(p)`:

- Curve: \(y^2 = x^3 + a x + b\) with `a = -3`, `b = 2`.
- Prime: a prime `p`
- Generator point `G` and public key `Q = d * G` are provided.
- The private key `d` is the **flag**

We must observe that this curve is *singular* i.e ` 4 a^3 + 27 b^2 = 0 `(discriminant is 0 mod `p`), a curve that does satisfy this sequality has a "problematic" point called a **singular point**.


The curve is not a proper elliptic curve group — it has a double root — and it becomes birationally equivalent to a much weaker group (a multiplicative group), allowing an easy reduction to a discrete-log problem in `GF(p)` that we can solve using the  Pohlig-Hellman algorithm (It's already implemented in sage's builtin `discrete_log` function).

---

## Attack strategy (high level)

1. Confirm the curve is singular by checking \((4a^3 + 27b^2) \bmod p = 0\).
2. Find the two roots of \(f(x)=x^3 + a x + b\) in `GF(p)`; identify which root has multiplicity 2 (double root) and which is single.
3. Compute \(t = r_d - r_s\). Compute a square root \(\sqrt{t}\) in `GF(p)`.
4. Apply the transform \(\phi\) to both `G` and `Q` to get `g` and `q`.
5. Solve the discrete log `n = discrete_log(q, g)` in the multiplicative group of `GF(p)` .
6. Convert the resulting integer `n` back to bytes to obtain the FLAG: `long_to_bytes(n)`.


---

## Solver (SageMath)

```python
p = 102360775616927576983385464260307534406913988994641083488371841417601237589487
a = -3
b = 2
assert (4*a^3 + 27*b^2) % p == 0

Gx = 1777671135698746847568710125129424132255529153914112337834835240247819869964
Gy = 6786424314307625790108882554225666781375821855884993473586521771737454762217

Qx =  19432656454182731855889820919208336901319861361037252686193041788066386070606
Qy =  33450649343988672730847959300492782097377595442476118891395518610180059585465

x = GF(p)["x"].gen()
f = x^3 + a*x + b
roots = f.roots()

assert len(roots) == 2
if roots[0][1] == 2:
    double_root = roots[0][0]
    single_root = roots[1][0]
else:
    double_root = roots[1][0]
    single_root = roots[0][0]


# map G and Q to the new "shifted" curve
Gx = (Gx - double_root)
Qx = (Qx - double_root)

# Transform G and Q into numbers g and q, such that q=g^n
t = double_root - single_root
t_sqrt = t.square_root()

def transform(x, y, t_sqrt):
    return (y + t_sqrt * x) / (y - t_sqrt * x)

g = transform(Gx, Gy, t_sqrt)
q = transform(Qx, Qy, t_sqrt)

# Find the private key n
found_key = discrete_log(q, g)
print(found_key)
from Crypto.Util.number import long_to_bytes
print(long_to_bytes(found_key).decode())
```
---
## Flag : 
```
shellmates{yee_s1n6ul4r_curv3!}
```
---

## References : 
- [ECDH](https://en.wikipedia.org/wiki/Elliptic-curve_Diffie%E2%80%93Hellman)
- [Pohling-Hellman algorithm](https://en.wikipedia.org/wiki/Pohlig%E2%80%93Hellman_algorithm)
- [Elleptic curves attacks implementations](https://github.com/elikaski/ECC_Attacks)