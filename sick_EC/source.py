from collections import namedtuple
from Crypto.Util.number import inverse, bytes_to_long
from secret import FLAG

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y
    def __repr__(self):
        return f"Point({self.x}, {self.y})"


O = Point(0, 0)


def check_point(P, p, a=-3, b=2):
    if P == O:
        return True
    else:
        return (P.y**2 - (P.x**3 + a*P.x + b)) % p == 0 and 0 <= P.x < p and 0 <= P.y < p


def point_inverse(P, p):
    if P == O:
        return P
    return Point(P.x, -P.y % p)


def point_addition(P, Q, p, a=-3, b=2):
    if P == O:
        return Q
    elif Q == O:
        return P
    elif Q == point_inverse(P, p):
        return O
    else:
        if P == Q:
            lam = (3*P.x**2 + a)*inverse(2*P.y, p)
            lam %= p
        else:
            lam = (Q.y - P.y) * inverse((Q.x - P.x), p)
            lam %= p
    Rx = (lam**2 - P.x - Q.x) % p
    Ry = (lam*(P.x - Rx) - P.y) % p
    R = Point(Rx, Ry)
    assert check_point(R, p, a, b)
    return R


def double_and_add(P, n, p, a=-3, b=2):
    Q = P
    R = O
    while n > 0:
        if n % 2 == 1:
            R = point_addition(R, Q, p, a, b)
        Q = point_addition(Q, Q, p, a, b)
        n = n // 2
    assert check_point(R, p, a, b)
    return R


def public_key():
    d = bytes_to_long(FLAG)
    return double_and_add(G, d, p, a, b)


p = 102360775616927576983385464260307534406913988994641083488371841417601237589487
a = -3
b = 2

# assert (4*a^3 + 27*b^2) % p == 0

Gx = 1777671135698746847568710125129424132255529153914112337834835240247819869964
Gy = 6786424314307625790108882554225666781375821855884993473586521771737454762217

G = Point(Gx, Gy)
assert check_point(G, p, a, b)
Q = public_key()

print("Q.x: ", Q.x)
print("Q.y: ", Q.y)


