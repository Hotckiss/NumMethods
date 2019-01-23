# -*- coding: utf-8 -*-
"""
Редактор Spyder

Это временный скриптовый файл.
"""

import numpy as np
import matplotlib.pyplot as graph
from collections import defaultdict
import math

# task 2
def Wa(z):
    sum = 0
    for i in range(1, 1100000):
        sum += 1.0 / (i * i - i - z)
    return sum
        

def plotWa(fn):
    a1 = []
    a2 = []
    
    for i in range(3, 197): # start point 3/100 = 0.03, end point: 197/100=1.97
        arg = i / 100.0
        a1.append(arg)
        a2.append(Wa(arg))
        
    graph.plot(a1, a2, label=fn)
    graph.savefig(fn)
    graph.close() 

def Wb(z):
    sum = 1-1.0/z
    for i in range(2, 120):
        sum += z / ((i * i - i) * (i * i - i - z))
    return sum    

def plotWb(fn):
    a1 = []
    a2 = []
    
    for i in range(3, 197): # start point 3/100 = 0.03, end point: 197/100=1.97
        arg = i / 100.0
        a1.append(arg)
        a2.append(Wb(arg))
        
    graph.plot(a1, a2, label=fn)
    graph.savefig(fn)
    graph.close() 

def diffW(z):
    return Wa(z) - Wb(z)

def plotDiffW(fn):
    a1 = []
    a2 = []
    
    for i in range(1, 199): # start point 1/100 = 0.01, end point: 199/100=1.99
        arg = i / 100.0
        a1.append(arg)
        a2.append(diffW(arg))
    
    max_diff_abs = max(abs(max(a2)), abs(min(a2)))
    print(max_diff_abs) # difference less than 2*eps
    print(max_diff_abs < 2.0 / 1000000)
    graph.xlabel('x')
    graph.ylabel('y')
    graph.plot(a1, a2, label=fn)
    graph.savefig(fn)
    graph.close() 

# task 3
def step(z, i):
    return z ** i / (2 * i - 1)

def eitk(z, steps):
    s1 = z 
    s2 = s1 + step(z, 2) 
    s3 = s2 + step(z, 3)
    diff2 = s3 - 2 * s2 + s1
    diff_s = s3 - (s3 - s2) ** 2 / diff2 
    for i in range(4, steps): 
        s1 = s2 
        s2 = s3 
        s3 += step(z, i)
        diff2 = s3 - 2 * s2 + s1
        if diff2 != 0: 
            diff_s = s3 - (s3 - s2) ** 2 / diff2 
    return diff_s

def generate_first_steps(z): # gen sj
    res = []
    s = z
    for i in range(2, 100): # 100 steps
        res.append(s)
        s += step(z, i)
    return res

def generate_diff_steps(prev): # [sj]->[s'j]->[s''j]..., длина списка сокращается всегда
    res = []
    for i in range(2, len(prev)):
        diff1 = prev[i] - prev[i-1]
        diff2 = prev[i] - 2 * prev[i-1] + prev[i-2]
        if diff2 != 0: 
            diff_s = prev[i] - diff1 ** 2 / diff2
            res.append(diff_s)
    return res
    
def eitk_iterated(z):
    start = generate_first_steps(z)
    while len(start) > 4: # can generate non-empty sequence
        start = generate_diff_steps(start)
    return start[-1]

def P2R(radii, angles):
    return radii * np.exp(1j*angles)

def convergecce_speed():
    x = []
    y = []
    
    for i in range(1, 100): # radius=1, angle from 0 to pi with step pi/100
        angle = math.pi * i / 100.0
        rad = 1
        z = P2R(rad, angle)
        ans = eitk(z, 10000)
        res = eitk(z, 300)
        x.append(i)
        y.append(abs(ans - res))
    graph.plot(x, y, label="speed.png")
    graph.savefig("speed.png")
    graph.close()

def convergecce_speed_iter():
    x = []
    y = []
    
    for i in range(1, 100): # radius=1, angle fron 0 to pi with step pi/100
        angle = math.pi * i / 100.0
        rad = 1
        z = P2R(rad, angle)
        ans = eitk(z, 10000)
        res = eitk_iterated(z)
        x.append(i)
        y.append(abs(ans - res))
    graph.plot(x, y, label="speed_iter.png")
    graph.savefig("speed_iter.png")
    graph.close()

def convergecce_speed_in_Radius():
    x = []
    y = []
    
    for i in range(1, 190): 
        arg = (i - 100) / 100.0
        z = arg
        if z == 0.0:
            continue
        ans = eitk(z, 10000)
        res = eitk_iterated(z)
        x.append(i)
        y.append(abs(ans - res))
    graph.plot(x, y, label="speed_iter_R.png")
    graph.savefig("speed_iter_R.png")
    graph.close()

if __name__ == "__main__":
    plotWa("Wa.png")
    plotWb("Wb.png")
    plotDiffW("DiffW.png")
    print('Iterated  ' + str(eitk_iterated(-0.9)))
    print('Simple    ' + str(eitk(-0.9, 300)))
    print('\n')
    print('Iterated  ' + str(eitk_iterated(1j)))
    print('Simple    ' + str(eitk(1j, 300)))
    print('\n')
    print('Iterated  ' + str(eitk_iterated(-1)))
    print('Simple    ' + str(eitk(-1, 300)))
    print('\n')
    print('Iterated  ' + str(eitk_iterated((1j - 1) / math.sqrt(2))))
    print('Simple    ' + str(eitk((1j - 1) / math.sqrt(2), 300)))
    convergecce_speed()
    convergecce_speed_iter()
    convergecce_speed_in_Radius()
