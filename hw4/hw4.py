import numpy as np
import matplotlib.pyplot as graph
import scipy.special as ss

ITERS = 200

def p(n, x):
    return ss.legendre(n)(x)

def dp(n, x):
    return n / (1 - x*x) * (p(n-1, x) - x*p(n, x))

def calc_roots(prev_roots):
    n = len(prev_roots) + 1
    pnts = []
    pnts.append(-1)
    pnts.extend(prev_roots)
    pnts.append(1)
    
    res = []
    
    for i in range(1, len(pnts)):
        prev_x = pnts[i-1]
        cur_x = pnts[i]
        x0=(prev_x+cur_x)/2
        for j in range(15):
            x0 = x0 - p(n, x0)/dp(n, x0)
        res.append(x0)
    return res

def simpson(f, a, b, M):
    h = (b - a) / M
    summ = f(a) + f(b) + 4 * f(a + h / 2)
    x = a + h
    for i in range(1, M):
        summ += 2 * f(x) + 4 * f(x + h / 2)
        x += h
    return summ * (h / 6)

def get_roots(n):
    r0 = []
    for i in range(n):
        r0 = calc_roots(r0)
    return r0

def gen_f(i, N, r):
    def f(x):
        ans = 1
        for k in range(N):
            if k != i:
                ans *= ((x - r[k])/(r[i]-r[k]))
        return ans
    return f

def calc_w(roots):
    n = len(roots)
    res = []
    for i in range(n):
        res.append(simpson(gen_f(i, n, roots), -1, 1, ITERS))
    return res

def func(x):
    return 1 / (9 * x ** 2 + 1)

def trapezioid(f, a, b, M):
    summ = (f(a) + f(b)) / 2
    h = (b - a) / M
    x = a + h
    for i in range(1, M):
        summ += f(x)
        x += h
    return summ * h

def calc_integral(f, a, b, n):
    roots = get_roots(n)
    wi = calc_w(roots)
    summ = 0
    for i in range(n):
        summ += wi[i] * f((b - a) / 2 * roots[i] + (a + b) / 2)
    
    return summ * (b - a) / 2

def plot_error():
    real_res = 0.9177579784724424
    y1 = []
    y2 = []
    y3 = []
    x = []
    for n in range(1, 26):
        print(n)
        i1 = calc_integral(func, -1, 5, n)
        i2 = simpson(func, -1, 5, n)
        i3 = trapezioid(func, -1, 5, n)
        err1 = np.log10(abs(i1-real_res))
        err2 = np.log10(abs(i2-real_res))
        err3 = np.log10(abs(i3-real_res))
        x.append(np.log10(n))
        y1.append(err1)
        y2.append(err2)
        y3.append(err3)
    graph.plot(x, y1, label='gauss')
    graph.plot(x, y2, label='simpson')
    graph.plot(x, y3, label='trap')
    graph.xlabel('log10(N)')
    graph.ylabel('log10(eps)')
    graph.title('diff')
    graph.legend()
    graph.savefig("t1_3logx.png")
    graph.close()
    
if __name__ == "__main__":
    #print(calc_integral(func, -1, 5, 20))
    plot_error()
    #print(calc_w(7))
    #print(np.polynomial.legendre.leggauss(7))
    #for i in range(10):
        #r0 = calc_roots(r0)
    #print()
    #calc_roots(3, [0, 0.5])
