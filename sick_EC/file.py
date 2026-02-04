def add(P, Q, a, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 != y2 or y1 == 0):
        return None
    if x1 == x2:
        m = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        m = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    return (x3, y3)

def mul(k, P, a, p):
    R0 = None
    R1 = P
    for bit in bin(k)[2:]:
        if bit == '0':
            R1 = add(R0, R1, a, p)
            R0 = add(R0, R0, a, p)
        else:
            R0 = add(R0, R1, a, p)
            R1 = add(R1, R1, a, p)
    return R0