# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 21:40:30 2019

@author: Андрей
"""


import numpy as np
import matplotlib.pyplot as graph
from collections import defaultdict
import math

#task 1 a

def func(x):
    x = np.longdouble(x)
    return x * (np.sin(2 * x))

def gen_Lk(x, xi):
    xi = np.longdouble(xi)
    x = np.longdouble(x)
    res = np.longdouble([])

    for i, a in enumerate(xi):
        top = np.longdouble(1)
        bottom = np.longdouble(1)
        for j, b in enumerate(xi):
            if i != j:
                top *= (x - b)
                bottom *= (a - b)
        res = np.append(res, top / bottom)

    return res

def calc_lagrange(x, xvals, yvals):
    Lk = gen_Lk(x, xvals)
    result = np.longdouble(0)
    for i, value in enumerate(yvals):
        result += np.longdouble(value) * Lk[i]
    return result

def calc_err(xvals, yvals, a, b, N):
    ans = 0
    left = a
    step = np.longdouble((b - a) / N)
    while left <= b:
        ans = max(ans, np.abs(calc_lagrange(left, xvals, yvals) - func(left)))
        left += step
    return ans

def gen_x_for_lagrange(x0, deg):
    return x0 - 5 + np.arange(0, deg + 1) * 10 / deg

def plot_error(x0, deg, N, fn):
    lagrange_x = gen_x_for_lagrange(x0, deg)
    lagrange_y = [func(pnt) for pnt in lagrange_x]
    xs = np.arange(x0 - 5, x0 + 5, 10.0 / N)
    ys = [np.abs(calc_lagrange(x, lagrange_x, lagrange_y) - func(x)) for x in xs]

    graph.clf()
    graph.plot(xs, ys)
    graph.xlabel('x')
    graph.ylabel('err')
    graph.savefig(fn)
    print('Max error  ' + str(calc_err(lagrange_x, lagrange_y, x0 - 5, x0 + 5, N)))

#task 1b
def task1b():
    xs = range(5, 51)
    x0 = 100
    ys = []
    for deg in range(5, 51):
        lagrange_x = gen_x_for_lagrange(x0, deg)
        lagrange_y = [func(pnt) for pnt in lagrange_x]
        ys.append(calc_err(lagrange_x, lagrange_y, x0 - 5, x0 + 5, 1000))
    graph.clf()
    graph.plot(xs, ys)
    graph.ylabel('max err')
    graph.xlabel('N')
    graph.savefig('1b.png')    

def task1bscale():
    xs = range(20, 51)
    x0 = 100
    ys = []
    for deg in range(20, 51):
        lagrange_x = gen_x_for_lagrange(x0, deg)
        lagrange_y = [func(pnt) for pnt in lagrange_x]
        ys.append(calc_err(lagrange_x, lagrange_y, x0 - 5, x0 + 5, 1000))
    graph.clf()
    graph.plot(xs, ys)
    graph.ylabel('max err')
    graph.xlabel('N')
    graph.savefig('1bscale.png')
    
def task1blog():
    xs = range(5, 51)
    x0 = 100
    ys = []
    for deg in range(5, 51):
        lagrange_x = gen_x_for_lagrange(x0, deg)
        lagrange_y = [func(pnt) for pnt in lagrange_x]
        ys.append(np.log10(calc_err(lagrange_x, lagrange_y, x0 - 5, x0 + 5, 1000)))
    graph.clf()
    graph.plot(xs, ys)
    graph.ylabel('max err')
    graph.xlabel('N')
    graph.savefig('1blog.png')  
    
#task 1с

def cheb(k, deg):
    return np.longdouble(math.cos(np.longdouble(math.pi) / 2 * (2 * k - 1) / deg))

def move_segment(a, b, t):
    return np.longdouble(0.5) * (a + b) + np.longdouble(0.5) * (b - a) * t

def gen_arr(N):
   return np.longdouble(np.array([np.cos((np.pi * (2 * k - 1)) / (2 * N)) for k in range(1, N + 1)]))

def task1c():
    xs = range(5, 51)
    x0 = 100
    ys = []
    for deg in range(5, 51):
        cheb_arr = gen_arr(deg+1)
        lagrange_x = [move_segment(95, 105, t) for t in cheb_arr]
        lagrange_y = [func(pnt) for pnt in lagrange_x]
        ys.append(calc_err(lagrange_x, lagrange_y, x0 - 5, x0 + 5, 1000))
    graph.clf()
    graph.plot(xs, ys)
    graph.ylabel('max err')
    graph.xlabel('N')
    graph.savefig('1c.png')

def task1c_both():
    xs = range(5, 51)
    x0 = 100
    ys = []
    ys1 = []
    for deg in range(5, 51):
        cheb_arr = gen_arr(deg+1)
        lagrange_x = [move_segment(95, 105, t) for t in cheb_arr]
        lagrange_y = [func(pnt) for pnt in lagrange_x]
        ys.append(calc_err(lagrange_x, lagrange_y, x0 - 5, x0 + 5, 1000))
        lagrange_x1 = gen_x_for_lagrange(x0, deg)
        lagrange_y1 = [func(pnt) for pnt in lagrange_x1]
        ys1.append(calc_err(lagrange_x1, lagrange_y1, x0 - 5, x0 + 5, 1000))
    graph.clf()
    graph.plot(xs, ys)
    graph.plot(xs, ys1)
    graph.ylabel('max err')
    graph.xlabel('N')
    graph.savefig('1c_both.png')

def task1c_both_vals():
    xs = range(5, 51)
    x0 = 100
    ys = []
    ys1 = []
    for deg in range(5, 51):
        cheb_arr = gen_arr(deg+1)
        lagrange_x = [move_segment(95, 105, t) for t in cheb_arr]
        lagrange_y = [func(pnt) for pnt in lagrange_x]
        ys.append(calc_lagrange(x0, lagrange_x, lagrange_y))
        lagrange_x1 = gen_x_for_lagrange(x0, deg)
        lagrange_y1 = [func(pnt) for pnt in lagrange_x1]
        ys1.append(calc_lagrange(x0, lagrange_x1, lagrange_y1))
    graph.clf()
    graph.plot(xs, ys)
    graph.plot(xs, ys1)
    graph.ylabel('value in x0')
    graph.xlabel('N')
    graph.savefig('1c_vals.png')

#task 1 d
def fm(x):
    return abs(x - 1)

def gen_x_for_lagrange_fm(x0, deg):
    return x0 - 1 + np.arange(0, deg + 1) * 2 / deg

def calc_err_fm(xvals, yvals, a, b, N):
    ans = 0
    left = a
    step = np.longdouble((b - a) / N)
    while left <= b:
        ans = max(ans, np.abs(calc_lagrange(left, xvals, yvals) - fm(left)))
        left += step
    return ans

def task1d_std():
    xs = range(5, 51)
    x0 = 1
    ys = []
    for deg in range(5, 51):
        lagrange_x = gen_x_for_lagrange_fm(x0, deg)
        lagrange_y = [fm(pnt) for pnt in lagrange_x]
        ys.append(np.log10(calc_err_fm(lagrange_x, lagrange_y, x0 - 1, x0 + 1, 1000)))
    graph.clf()
    graph.plot(xs, ys)
    graph.ylabel('max err')
    graph.xlabel('N')
    graph.savefig('1d_stdlog.png') 

def task1d_cheb():
    xs = range(5, 51)
    ys = []
    for deg in range(5, 51):
        cheb_arr = gen_arr(deg+1)
        lagrange_x = [move_segment(0, 2, t) for t in cheb_arr]
        lagrange_y = [fm(pnt) for pnt in lagrange_x]
        ys.append(calc_err_fm(lagrange_x, lagrange_y, 0, 2, 1000))
    graph.clf()
    graph.plot(xs, ys)
    graph.ylabel('max err')
    graph.xlabel('N')
    graph.savefig('1d_cheb.png')
 
def task1d_all():
    xs = range(5, 100)
    ys1 = []
    ys2 = []
    ys3 = []
    ys4 = []
    for deg in range(5, 100):
        """cheb_arr = gen_arr(deg+1)
        lagrange_x1 = [move_segment(0, 2, t) for t in cheb_arr]
        lagrange_y1 = [fm(pnt) for pnt in lagrange_x1]
        ys1.append(np.log10(calc_err_fm(lagrange_x1, lagrange_y1, 0, 2, 1000)))
        
        lagrange_x2 = gen_x_for_lagrange_fm(1, deg)
        lagrange_y2 = [fm(pnt) for pnt in lagrange_x2]
        ys2.append(np.log10(calc_err_fm(lagrange_x2, lagrange_y2, 0, 2, 1000)))
        
        lagrange_x3 = [move_segment(95, 105, t) for t in cheb_arr]
        lagrange_y3 = [func(pnt) for pnt in lagrange_x3]
        ys3.append(np.log10(calc_err(lagrange_x3, lagrange_y3, 95, 105, 1000)))"""
        
        lagrange_x4 = gen_x_for_lagrange(100, deg)
        lagrange_y4 = [func(pnt) for pnt in lagrange_x4]
        ys4.append((calc_err(lagrange_x4, lagrange_y4, 95, 105, 1000)))
    graph.clf()
    #graph.plot(xs, ys1)
    #graph.plot(xs, ys2)
    #graph.plot(xs, ys3)
    graph.plot(xs, ys4)
    #graph.legend(['fM cheb', 'fM uni', 'fS cheb', 'fS uni'], loc='upper left')
    graph.ylabel('max err')
    graph.xlabel('N')
    graph.savefig('1d_4.png')

if __name__ == "__main__":
    task1d_all()
    #task1d_std()
    #task1d_cheb()
    #task1с()
    #task1с_both()
    #task1c_both_vals()
    #task1b()
    #task1blog()
    #task1bscale()
    #plot_error(100, 5, 1000, 'err5.png')
    #plot_error(100, 10, 1000, 'err10.png')
    #plot_error(100, 15, 1000, 'err15.png')
