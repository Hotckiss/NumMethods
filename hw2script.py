# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 21:40:30 2019

@author: Андрей
"""


import numpy as np
import matplotlib.pyplot as graph
import math
import time
    
def getA(n): 
    return np.fromfunction(lambda i, j: 1 / (i + j + 1), (n, n))

def forwardIterations(A, n, iters=1000):
    prev = np.ones(n)
    cur = A @ prev
    res = []
    for i in range(iters):
        prev = cur
        cur = A @ cur
        norm = np.linalg.norm(cur)
        if norm > 1e9:
            cur /= norm
        res.append(np.dot(prev, cur) / np.dot(prev, prev))
    return res

def convSpeed():
    itersNum = 300
    iters = forwardIterations(4, itersNum)
    realResult = max(np.linalg.eig(getA(4))[0])
    
    x = []
    y = []
    for i in range(itersNum):
        x.append(np.log10(i))
        y.append((abs(realResult - iters[i])))
        
    graph.plot(x, y)
    graph.xlabel('log10(N)')
    graph.ylabel('log10(eps)')
    graph.title('max eig')
    #graph.savefig('maxconv.png')
    graph.show()

def plot_max_eigenvalue_draft():
    x = []
    y2 = []
    for i in range(1, 501):
        if i == 1:
            print(i ** (max(np.linalg.eigh(getA(i))[0]) / np.pi))
        if i == 500:
            print(i ** (max(np.linalg.eigh(getA(i))[0]) / np.pi))
        x.append(i)
        y2.append(i ** (max(np.linalg.eigh(getA(i))[0]) / np.pi))
       
    graph.plot(x, y2, label='numpy')
    graph.xlabel('N')
    graph.ylabel('max eig')
    graph.title('Max eig')
    graph.legend()
    graph.savefig('draft.png')
    graph.show()

def func(n):
    return np.pi * np.log(109.15 * n / 499 + 389.85 / 499) / np.log(n)

def plotDiff():
    x = []
    y2 = []
    for i in range(1, 1001):
        x.append(i)
        y2.append(abs(max(np.linalg.eigh(getA(i))[0]) - func(i) ))
       
    graph.plot(x, y2, label='numpy')
    graph.xlabel('N')
    graph.ylabel('|ans - func|')
    graph.title('diff')
    graph.legend()
    graph.savefig('draft.png')
    graph.show()

def plot_max_eigenvalue():
    x = []
    y1 = []
    y2 = []
    for i in range(1, 1001):
       x.append((i))
       y1.append(i ** (forwardIterations(i)[-1] / 3))
       y2.append(i ** (max(np.linalg.eigh(getA(i))[0]) / 3))
       
    graph.plot(x, y1, label='iterations')
    graph.plot(x, y2, label='numpy')
    graph.xlabel('N')
    graph.ylabel('N ** (eigenvalue / 3)')
    graph.title('Max eigenvalue')
    graph.legend()
    graph.show()

def getMinEig(n, iters=1000, alpha=None):
    alpha = forwardIterations(getA(n), n, iters)[-1] + 1e-6 if alpha is None else alpha
    A = getA(n) - alpha * np.identity(n)
    
    return forwardIterations(A, n, iters)[-1] + alpha

def plotMinEig():
    x = []
    y1 = []
    y2 = []
    for i in range(1, 501):
       x.append(i)
       y1.append(np.log10(getMinEig(i)))
       y2.append(np.log10(min(np.linalg.eigh(getA(i))[0])))
       
    graph.plot(x, y1, label='iterations')
    graph.plot(x, y2, label='numpy')
    graph.xlabel('N')
    graph.ylabel('log10(min eig)')
    graph.title('min eig')
    graph.legend()
    graph.savefig('minn.png')
    graph.show()

def getK(n):
    return forwardIterations(getA(n), n)[-1] / getMinEig(n)

def plotK():
    x = []
    y1 = []
    y2 = []
    for i in range(1, 11):
       x.append(i)
       y1.append(np.log10(getK(i)))
       y2.append(np.log10(np.linalg.cond(getA(i))))
       
    graph.plot(x, y1, label='iterations')
    graph.plot(x, y2, label='numpy')
    graph.xlabel('N')
    graph.ylabel('log10(kappa)')
    graph.title('kappa')
    graph.legend()
    graph.savefig('cap.png')
    graph.show()

def calcNext(cur, prev, alpha):
    return np.dot(prev, cur) / np.dot(prev, prev) + alpha

def eitk(n, iters = 1000):
    alpha = forwardIterations(getA(n), n, iters)[-1] + 1e-6
    A = getA(n) - alpha * np.identity(n)
    
    prev = np.ones(n)
    cur = A @ prev
    alpha = forwardIterations(getA(n), n)[-1] + 1e-6
    
    s1 = calcNext(cur, prev, alpha)
    prev = cur
    cur = A @ prev
    
    s2 = calcNext(cur, prev, alpha)
    prev = cur
    cur = A @ prev
    
    s3 = calcNext(cur, prev, alpha)
    diff = s3 - (s3 - s2) ** 2 / (s3 - 2 * s2 + s1)
    
    for i in range(iters):
        prev = cur
        cur = A @ cur
        norm = np.linalg.norm(cur)
        if norm > 1e9:
            cur /= norm
        s1 = s2
        s2 = s3
        s3 = calcNext(cur, prev, alpha)
        if s3 - 2 * s2 + s1 != 0:
            diff = s3 - (s3 - s2) ** 2 / (s3 - 2 * s2 + s1)
    return diff

def plotMinEig2():
    x = []
    y1 = []
    y2 = []
    y3 = []
    for i in range(1, 6):
       x.append(i)
       y1.append(np.log10(getMinEig(i)))
       y2.append(np.log10(min(np.linalg.eigh(getA(i))[0])))
       y3.append(np.log10(eitk(i)))
       
    graph.plot(x, y1, label='iterations')
    graph.plot(x, y2, label='numpy')
    graph.plot(x, y3, label='eitken')
    graph.xlabel('N')
    graph.ylabel('log10(min eig)')
    graph.title('min eig')
    graph.legend()
    graph.savefig('minne.png')
    graph.show()

if __name__ == "__main__":
    #plotK()
    plotMinEig2()
    #convSpeed()
    #plotDiff()
    #print(max(np.linalg.eigh(getA(10000))[0]))