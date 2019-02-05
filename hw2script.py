# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 21:40:30 2019

@author: Андрей
"""


import numpy as np
import matplotlib.pyplot as graph
from collections import defaultdict
import math

#task1
def f(x):
    return np.sin(x)

def df(x):
    return np.cos(x)
    
def p(a, h, x0):
    x1 = x0 - h / a
    y1 = f(x1)
    x2 = x0
    y2 = f(x2)
    x3 = x0 + h * a
    y3 = f(x3)
    return 1 / h * ((y2 - y1) * (a - 1 / a) + (y3 - y1) / (a ** 3 + a))

def plot_error(alph, x0):
    x = []
    y = []
    step = 10 ** (-3)
    h = 10 ** (-8)
    b = 100
    ans = df(x0)
    while h <= b:
        x.append(h)
        y.append(abs(p(alph, h, x0) - ans))
        h += step
    graph.title("x0 = " + str(x0))
    graph.plot(x, y)
    graph.ylabel("abs diff")
    graph.xlabel("h")
    
def task1():
    x0 = 2
    arr = [1, 1.2, 1.4, 1.6, 1.8, 2]

    for a in arr:
        plot_error(a, x0)
    graph.legend(list(map(lambda a: "a = " + str(a), arr)))
    graph.savefig("t1e.png")
    graph.close()
    
#task2
    
def func(x):
    return 1.0 / (9 * (x ** 2) + 1)

def func_int_indefinite(x):
    return np.arctan(3 * x) / 3

def func_int_definite(a, b):
    return func_int_indefinite(b) - func_int_indefinite(a)

def trapezioid(f, a, b, M):
    summ = (f(a) + f(b)) / 2
    h = (b - a) / M
    x = a + h
    for i in range(1, M):
        summ += f(x)
        x += h
    return summ * h

def simpson(f, a, b, M):
    h = (b - a) / M
    summ = f(a) + f(b) + 4 * f(a + h / 2)
    x = a + h
    for i in range(1, M):
        summ += 2 * f(x) + 4 * f(x + h / 2)
        x += h
    return summ * (h / 6)

def plot_error_task2(a, b):
    start = 1
    finish = 2000
    x = []
    yt = []
    ys = []
    ans = func_int_definite(a, b)
    while start <= finish:
        x.append(start)
        yt.append(np.log10(abs(ans - trapezioid(func, a, b, start))))
        ys.append(np.log10(abs(ans - simpson(func, a, b, start))))
        start += 1
    graph.plot(x, yt, label='trapezioid')
    graph.plot(x, ys, label='simpson')
    graph.xlabel('M')
    graph.ylabel('log10(eps)')
    graph.title('diff')
    graph.legend()
    graph.savefig("t2a.png")
    graph.close()

#tast2b
def calc_h_runge(eps, a, b):
    h2 = 10 ** (-3)
    h1 = 2 * h2
    sh2 = trapezioid(func, a, b, int((b - a) / h2))
    sh1 = trapezioid(func, a, b, int((b - a) / h1))
    c = abs(1 / (3 * h2 ** 2) * (sh2 - sh1))
    return np.sqrt(eps / c)

def calc_h(eps, a, b):
    start = 1
    finish = 2000
    x = []
    yt = []
    ans = func_int_definite(a, b)
    while start <= finish:
        x.append(start)
        yt.append(abs(ans - trapezioid(func, a, b, start)))
        start += 1 
    for i in range(1, len(yt)):
        if yt[i] < eps:
            return (b - a) / x[i]

#task2c
def gen_f(i, N):
    def f(q):
        ans = 1
        for k in range(1, N + 1):
            if k != i:
                ans *= q - (k - 1)
        return ans
    return f

def fac(n):
    return math.factorial(n)

def calc_weights(a, b, N):
    h = (b - a) / (N - 1)
    arr = []
    for i in range(1, N + 1):
        li = simpson(gen_f(i, N), 0, N - 1, 1000)
        li *= (-1) ** (N - i) * h / fac(i - 1) / fac(N - i) 
        arr.append(li)
    return arr          

def neg_weight_N(a, b):
    N = 2
    w = min(calc_weights(a, b, N))
    while w >= 0:
        N += 1
        w = min(calc_weights(a, b, N))
    return N

def min_weight(a, b, N):
    return min(calc_weights(a, b, N))

def plot_min_weight(a, b):
    x = []
    y = []
    
    for n in range(2, 30):
        x.append(n)
        y.append(min_weight(a, b, n))
    graph.plot(x, y)
    graph.ylabel("min_W")
    graph.xlabel("N")
    graph.savefig("t2c.png")
    graph.close()
    
if __name__ == "__main__":
    #task1()
    """print("Math answer:")
    print(func_int_definite(-1, 5))
    print("Trapezioid answer:")
    print(trapezioid(func, -1, 5, 1000))
    print("Simpson answer:")
    print(simpson(func, -1, 5, 1000))
    print("|math - trapezioid|:")
    print(abs(trapezioid(func, -1, 5, 1000) - func_int_definite(-1, 5)))
    print("|math - simpson|:")
    print(abs(simpson(func, -1, 5, 1000) - func_int_definite(-1, 5)))"""
    #plot_error_task2(-1, 5)
    """print("Runge:")
    print(calc_h_runge(10 ** (-6), -1, 5))
    print("Real:")
    print(calc_h(10 ** (-6), -1, 5))
    print("Diff:")
    print(abs(calc_h_runge(10 ** (-6), -1, 5) - calc_h(10 ** (-6), -1, 5)))"""
    print("First N with negative weight:")
    print(neg_weight_N(-1, 1))
    plot_min_weight(-1, 1)