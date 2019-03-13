# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 21:40:30 2019

@author: Андрей
"""


import numpy as np
import matplotlib.pyplot as graph
import math

def p(r, l):
    return l * (l + 1) * r ** (-2) - 2 / r
     
def rho(r):
    return 1

def matr(l, R, N):
    h = R / N
    res = np.zeros((N - 1, N - 1))
    #rho=1 always
    for i in range(0, N - 1):
        if i != 0:
            res[i][i - 1] = h ** (-2)
        res[i][i] = - (2 / h ** 2 + p((i+1) * h, l))
        if i != N - 2:
            res[i][i + 1] = h ** (-2)
    return -res / 2

def task1():
    R = 100
    M = matr(0, R, 5000)
    eig = np.linalg.eig(M)
    spec, vec = eig
    idx = np.argsort(spec)[:5]
    print("5 eig vals")
    print(spec[idx])
    print("5 eig functions")
    print(vec[idx])

def drawPlot(eigs, R):
    N = (len(eigs) + 1)
    h = R / N
    x = [i * h for i in range(N + 1)]
    y = [0]
    y.extend(eigs)
    y.append(0)
    
    graph.plot(x, y)
    graph.ylabel("eig_f(N)")
    graph.xlabel("N")
    graph.show()
    
def task1graph():
    R = 100
    M = matr(0, R, 5000)
    eig = np.linalg.eig(M)
    spec, vec = eig
    idx = np.argsort(spec)[:5]
    print(spec[idx][0])
    drawPlot(vec[idx][0], R)
    print(spec[idx][1])
    drawPlot(vec[idx][1], R)
    print(spec[idx][2])
    drawPlot(vec[idx][2], R)
    print(spec[idx][3])
    drawPlot(vec[idx][3], R)
    print(spec[idx][4])
    drawPlot(vec[idx][4], R)
    
    #graph.savefig("t1g.png")
    #graph.close()

def num_a(l, R, N):
    h = R / N
    A = np.zeros((N - 1, N - 1))
    for i in range(0, N - 1):
        A[i][i] = - (2 * h ** (-4) +  p((i + 1) * h, l) - 1 / 6 * p((i + 1) * h, l) * h ** (-2))
        if i != 0:
            A[i][i - 1] = 1 * h ** (-4) - 1 / 12 * p(i * h, l) * h ** (-2)
        if i != N - 2:
            A[i][i + 1] = 1 * h ** (-4) - 1 / 12 *  p((i + 2) * h, l) * h ** (-2)
    return A

def num_b(l, R, N):
    h = R / N
    B = np.zeros((N - 1, N - 1))
    for i in range(0, N - 1):
        B[i][i] = 1 - 1 / 6 * h ** (-2)
        if i != 0:
            B[i][i - 1] = h ** (-2) / 12
        if i != N - 2:
            B[i][i + 1] = h ** (-2) / 12
              
    return B

def numerov(l, R, N):
    return -(np.linalg.inv(num_b(l, R, N)) @ num_a(l, R, N)) / 2

def task2plot():
    R = 1000
    count = 5
    result = [-1 / (2 * (n + 1) ** 2) for n in range(count)]
    x = []
    y1 = []
    y2 = []
    for i in range(6, 300):
        x.append(i)
        m1 = matr(0, R, i)
        m2 = numerov(0, R, i)
        y1.append(np.log10(max(abs(np.sort(np.linalg.eig(m1)[0])[:5] - result))))
        y2.append(np.log10(max(abs(np.sort(np.linalg.eig(m2)[0])[:5] - result))))
    graph.plot(x, y1, y2)
    graph.ylabel("log10(err)")
    graph.xlabel("N")
    graph.legend(("standard", "numerov"))
    graph.show()

if __name__ == "__main__":
    task1graph()
    #task2plot()
    """R = 10
    M = matr(0, R, 1000)
    eig = np.linalg.eig(M)
    spec = eig[0]
    vec = eig[1]
    indexes = np.argsort(spec)[:5]
    
    #print(matr(0, 10, 5))
    print(spec)
    print(spec[indexes][0])
    task1graph(vec[indexes][0], R)"""
