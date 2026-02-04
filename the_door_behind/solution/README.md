## Challenge Overview

In this challenge we want to factorize `n` knowing the way one of its primes has been constructed ` p = (635 * s ^ 2 + 1) / 4 `

The Cheng’s 4p − 1 factorization method can be used to solve the challenge as the algorithm described in the paper (Algorithm 1 , section 3 [link to the paper](https://www.scitepress.org/Papers/2019/77866/77866.pdf))

Chengs’s 4p − 1 method works on an elliptic curve (EC) \( E(\mathbb{Z}_n) \) (a set of points \((x,y) \in \mathbb{Z}_n \times \mathbb{Z}_n\)) defined by \( a,b \in \mathbb{Z}_n \) and the Weierstrass equation:

\[
E(\mathbb{Z}_n) : y^2 = x^3 + a x + b,
\]

where \( n \) is the number to be factored and \( 4a^3 + 27b^2 \neq 0 \pmod{n} \).  
It works due to the natural mapping : 

\[
(\mathbb{Z}_n \xrightarrow{\text{mod } p} \mathbb{F}_p)
\]

that induces a homomorphism \( E(\mathbb{Z}_n) \to E(\mathbb{F}_p) \) by reducing the coordinates modulo \( p \).

The methods compute the multiple \( mP \) for a point \( P \) on \( E(\mathbb{Z}_n) \).  

If \( m \in \mathbb{Z} \) is a multiple of \( E(\mathbb{F}_p) \) (order of the EC over \( \mathbb{F}_p \)), then the computation of a certain inversion  
modulo \( n \) fails (inversion of a multiple of \( p \) for \( p \mid n \)) during the scalar multiplication, which reveals \( p \).




## Sage Solver

One proposed solver is the following : 

```python

# adapted from : https://github.com/crocs-muni/cm_factorization/blob/master/cm_factor.sage

import time
import sys
import traceback
from sage.misc.prandom import randrange
from sage.parallel.decorate import fork
from Crypto.Util.number import long_to_bytes , inverse
class AlgException(Exception):
    pass

class NotInvertibleException(AlgException):
    pass

class NotFactoredException(AlgException):
    pass

class FactorRes(object):
    def __init__(self, r=None, c=None, u=None, a=None, th=None, tq=None):
        self.r = r
        self.c = c
        self.u = u
        self.a = a
        self.time_hilbert = th
        self.time_q = tq
        self.time_a = None
        self.time_last_div = None
        self.time_last_gcd = None
        self.time_last_nrm = None
        self.time_total = 0
        self.time_agg_div = 0
        self.time_agg_gcd = 0
        self.time_agg_nrm = 0
        self.time_qinv_char_poly = 0
        self.time_qinv_xgcd = 0
        self.time_qinv_res = 0
        self.rand_elem = None
        self.fact_timeout = None
        self.out_of_time = False
        self.use_quinv2 = False
        self.use_cheng = False

def is_undef(result):
    return result in ['NO DATA (timed out)', 'NO DATA', 'INVALID DATA', 'INVALID DATA ', None]

def class_number(d):
    k = QuadraticField(-d, 'x')
    return k.class_number()

def xgcd(f, g, N=1):
    toswap = False
    if f.degree() < g.degree():
        toswap = True
        f, g = g, f
    r_i = f
    r_i_plus = g
    r_i_plus_plus = f
    s_i, s_i_plus = 1, 0
    t_i, t_i_plus = 0, 1
    while (True):
        lc = r_i.lc().lift()
        lc *= r_i_plus.lc().lift()
        lc *= r_i_plus_plus.lc().lift()
        divisor = gcd(lc, N)
        if divisor > 1:
            print('Divisor of %s is %s'%(N,divisor))
            return divisor, None, None
        q = r_i // r_i_plus
        s_i_plus_plus = s_i - q * s_i_plus
        t_i_plus_plus = t_i - q * t_i_plus
        r_i_plus_plus = r_i - q * r_i_plus
        if r_i_plus.degree() <= r_i_plus_plus.degree() or r_i_plus_plus.degree() == -1:
            if toswap == True:
                assert (r_i_plus == s_i_plus * f + t_i_plus * g)
                return r_i_plus, t_i_plus, s_i_plus
            else:
                assert (r_i_plus == s_i_plus * f + t_i_plus * g)
                return r_i_plus, s_i_plus, t_i_plus
        r_i, r_i_plus = r_i_plus, r_i_plus_plus
        s_i, s_i_plus = s_i_plus, s_i_plus_plus
        t_i, t_i_plus = t_i_plus, t_i_plus_plus
        check_res = r_i == s_i * f + t_i * g
        if not check_res:
            print('Assertion error: %s, %s != %s * %s + %s * %s' % (check_res, r_i, s_i, f, t_i, g))
            raise ValueError('xgcd assertion error')

def Qinverse(Q, a, N):
    j = Q.gens()[0]
    deg = j.charpoly('X').degree()
    A = Q(a).matrix()
    det_a = det(A)
    print('DetA: %s' % det_a)
    factor = gcd(int(det_a), N)
    if factor!=1:
        raise ZeroDivisionError(a)
    else:
        Y = vector([1] + (deg-1)*[0])
        X = A.solve_left(Y)
        jvec = vector([j^i for i in [0..deg-1]])
        Xj = jvec*X
        return Xj

def Qinverse2 (Hx, a, N, time_res):
    ts = time.time()
    r,s,t = xgcd(a.lift(), Hx, N)
    txgcd = time.time()
    if (s,t) == (None, None):
        res = r, 0
    else:
        rinv = r[0]^(-1)
        res = 1, s * rinv
    tres = time.time()
    time_res.time_qinv_char_poly = 0
    time_res.time_qinv_xgcd = txgcd - ts
    time_res.time_qinv_res = tres - txgcd
    return res

def CMfactor(D, N, verb = 1, ctries=10, utries=10, fact_time=None, use_quinv2=False, use_cheng=False):
    ts = time.time()
    Hx = hilbert_class_polynomial(-D)
    tth = time.time()
    if verb == 1:
        print('Hilbert polynomial computed for -%s!' % D)
    res = FactorRes()
    res.use_quinv2 = use_quinv2
    res.use_cheng = use_cheng
    ZN = Integers(N)
    R.<x> = PolynomialRing(ZN)
    ttq = time.time()
    try:
        if use_quinv2:
            Hx = R(Hx)
            Q.<j> = QuotientRing(R, R.ideal(Hx))
            gcd, inverse = Qinverse2(Hx, 1728 - j, N, res)
            if gcd == 1:
                a = Q(j * inverse)
            else:
                print('Early factor found: %s' % gcd)
                res.r = gcd
                return res
        else:
            Q.<j> = ZN.extension(Hx)
            a = j * Qinverse(Q, 1728 - j, N)
    except ZeroDivisionError as noninv:
        print("is not invertible in Q! %s" % noninv)
        raise NotInvertibleException()
    if verb == 1:
        print('Q constructed')
        print('a computed: %s' % a)
    tta = time.time()
    res.time_agg_div = 0
    res.time_agg_gcd = 0
    res.time_hilbert = tth - ts
    res.time_q = ttq - tth
    res.time_a = tta - ttq
    res.a = None
    core_fnc = CMfactor_core
    if fact_time:
        time_left = fact_time - (tta - ts)
        res.fact_timeout = time_left
        core_fnc = fork(CMfactor_core, time_left)
    cres = core_fnc(N, ctries, utries, a, Q, ZN, Hx, res, use_cheng=use_cheng)
    if is_undef(cres):
        res.out_of_time = True
    else:
        res = cres
    tdone = time.time()
    res.time_total = tdone - ts
    return res

def CMfactor_core(N, ctries, utries, a, Q, ZN, Hx, res, use_cheng=False):
    is_done = False
    for c in [1..ctries]:
        E = EllipticCurve(Q, [0, 0, 0, 3 * a * c ^ 2, 2 * a * c ^ 3])
        for u in [1..utries]:
            tcs = time.time()
            rand_elem = ZN.random_element()
            res.rand_elem = int(rand_elem)
            w = E.division_polynomial(N, Q(rand_elem), two_torsion_multiplicity=0)
            ttdiv = time.time()
            print('Division polynomial done')
            if use_cheng:
                poly_gcd = xgcd(w.lift(), Hx, N)[0]
                ttnrm = time.time()
                r = gcd(ZZ(poly_gcd), N)
                ttgcd = time.time()
            else:
                nrm = w.norm()
                ttnrm = time.time()
                r = gcd(nrm, N)
                ttgcd = time.time()
            res.time_agg_div += ttdiv - tcs
            res.time_agg_gcd += ttgcd - ttdiv
            res.time_agg_nrm += ttnrm - ttdiv
            res.c = int(c)
            res.u = int(u)
            res.time_last_div = ttdiv - tcs
            res.time_last_gcd = ttgcd - ttdiv
            res.time_last_nrm = ttnrm - ttdiv
            if r > 1 and r != N:
                res.r = int(r)
                print('A factor of N: %s' % r)
                print('c: %s, u: %s' % (c, u))
                is_done = True
                break
            else:
                print('u failed: %s, next_attempt' % u)
        if is_done:
            break
    return res

def factor_cm(n, D, timeout=4*60, use_quinv2=True, use_cheng=True):
    sys.setrecursionlimit(50000)
    try:
        factor_timeout = timeout if timeout > 0 else None
        res = CMfactor(D, n, 1, fact_time=factor_timeout, use_quinv2=use_quinv2, use_cheng=use_cheng)
        if res.out_of_time or res.r is None or res.r <= 1 or n % res.r != 0:
            return None
        return res.r
    except Exception as e:
        print('Exception: %s' % e)
        if traceback:
            traceback.print_exc()
        return None

n =  Integer(1075046540657649889830822121812142481004900172589716910676941896300770112035833613806943324384856820537897642401504143425686011448660060134792939077259038755691060530460558072805147285506952067010991226349442753805824100134920314311905880232412756454680642140745011552508732431522967978306694355253216358975962699286141678239167514090112723713287088849533847541422311491844910152907522756041573985948763652518867553950597150013925636864933731902866947992254985429069005851872404053051975600226557024503999451256001732483225831966490363737430781147415195494487087477816053880778763615045705354990651038896432101900012602726378987336576840111047782242134571810578373239353239172018305878831018345751398965948495288535278845026057520222933399111518226095839653858538485497378068224665846658804433674863951589307792855529102926501073685589937791948826524914053243686274145601373396621139072089529892020859934892415116649460722031598663894390039790543516075904005918434845288820227205969757081836408090093147146763294040846027065610693225707968917591634725821054425085248698254617316554608952369279836699837869459028513888610378065890545726962865597637524771702347457637080480267895276614358074580938032253251343967631640865137253218743513)
e =  Integer(65537)
c =  Integer(414366691287960539389195595941205127064547999071015644891372937176883390260599519036038912788656021248236298194094763234240142668908315069505706177712822636189511302316617546124677192585145374101625250900530277691686909425848162027274126259750504638354431835111397419325368697062354891110281305616502658399504506546913554165165849035322545127669434079632920391172385602510536512607978541248555069301620300599215311117535027358228712064738411671908451948420348749915435878379489782796515867604138989786148043012735388026639666449851371396271979722344460780691826522682142365408041917380348362833219937795466897468181114132358063598684832633473509767856739005066592389245702046284886062480600569318195985633310095140691051828323920325216724427213360485640192330143468602849206250023512129671709225316028009724222764105512966765320938897343530385757575927860353095360764083077595279891020556404407338845274154723966630228107378806582769755866167544396297049149524907866320007419419715075004844087711631554299806099503017061816890051777073401868038179839535046726975744966378338195969002115284448249854934945172448988566474374143657817376975952293954191139539097551465116832362974748347146098109725634590076403972159688204977620545966698)
D =  Integer(635)

p = factor_cm(n,D)
if p is None :
    print("Error factoring")
    exit()
q = n // p
phi = (p - 1) * (q - 1)
d = inverse(e, phi)
m = pow(c, d, n)
flag = long_to_bytes(m)
print("FLAG : " , flag.decode())
```


---

Otherwise you  can just use the tool in the repo linked [here](https://github.com/crocs-muni/cm_factorization/tree/master) 