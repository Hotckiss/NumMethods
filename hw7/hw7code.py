# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 21:40:30 2019

@author: Андрей
"""


import numpy as np
import matplotlib.pyplot as graph
import math

T=10
#a=[1, np.sqrt(2), 2, 2 * np.sqrt(2), 3, 4, math.sqrt(20), 5, 6, 7, 8, 16]
a=[1, math.sqrt(20)]
#a=[8]
N=[10, 50, 100]

def euler(f0, f1, y0, y1, N):
    h = T / N
    x = 0
    y = [y0]
    z = [y1]
    
    for i in range(N):
        yc = y[-1] + h * f0(x, y[-1], z[-1])
        zc = z[-1] + h * f1(x, y[-1], z[-1])
        y.append(yc)
        z.append(zc)
        x += h
    
    return y

def task2():
    for cur_a in a:
        for cur_n in N:
            xs = [i * (T / cur_n) for i in range(cur_n + 1)]
            ys = euler(lambda x, y, z: z, lambda x, y, z: cur_a ** 2 * y, 1, -cur_a, cur_n)
            graph.plot(xs, ys, label = "N = " + str(cur_n))
            graph.ylabel("y")
            graph.xlabel("x")
            
        graph.title("T = 10; a = " + ("1" if cur_a == 1 else "sqrt(20)"))
        graph.legend()
        graph.savefig("hw7_t1_a" + ("1" if cur_a == 1 else "sqrt(20)") + ".png")
        graph.close()

def sol(a):
    return lambda x: np.e ** (-a * x)

def task3():
    for cur_a in a:
        xs = []
        ys = []
        for n in range(1, 201):
            xs.append(n)
            h = T / n
            math_sol = [sol(cur_a)(i * h) for i in range(n + 1)] 
            solution = euler(lambda x, y, z: z, lambda x, y, z: cur_a ** 2 * y, 1, -cur_a, n)
            ys.append(np.log10(max([abs(solution[i] - math_sol[i]) for i in range(len(math_sol))])))
        graph.plot(xs, ys)
        graph.ylabel("log10(err)")
        graph.xlabel("N")
        graph.plot([1, 200], [-3, -3])
        graph.title("T = 10; a = " + ("1" if cur_a == 1 else "sqrt(20)"))
        graph.show()
        #graph.savefig("hw7_t3_a" + ("1" if cur_a == 1 else "sqrt(20)") + ".png")
        #graph.close()

def runge(f0, f1, y0, y1, N):
    beta = 0.5
    h = T / N
    x = 0
    y = [y0]
    z = [y1]
    for i in range(N):
        yc = y[-1] + h * ( (1 - beta) * f0(x, y[-1], z[-1]) +
                                beta * f0(x + h / (2 * beta), 
                                         y[-1] + h / (2 * beta) * f0(x, y[-1], z[-1]),
                                         z[-1] + h / (2 * beta) * f1(x, y[-1], z[-1])) )
        zc = z[-1] + h * ( (1 - beta) * f1(x, y[-1], z[-1]) +
                                beta * f1(x + h / (2 * beta), 
                                         y[-1] + h / (2 * beta) * f0(x, y[-1], z[-1]),
                                         z[-1] + h / (2 * beta) * f1(x, y[-1], z[-1])) )
        y.append(yc)
        z.append(zc)
        x += h
    return y

def runge4(f0, f1, y0, y1, N):
    h = T / N
    x = 0
    y = [y0]
    z = [y1]
    for i in range(N):
        yk1 = f0(x, y[-1], z[-1])
        zk1 = f1(x, y[-1], z[-1])
        yk2 = f0(x + h / 2, y[-1] + h * yk1 / 2, z[-1] + h * zk1 / 2)
        zk2 = f1(x + h / 2, y[-1] + h * yk1 / 2, z[-1] + h * zk1 / 2)
        yk3 = f0(x + h / 2, y[-1] + h * yk2 / 2, z[-1] + h * zk2 / 2)
        zk3 = f1(x + h / 2, y[-1] + h * yk2 / 2, z[-1] + h * zk2 / 2)
        yk4 = f0(x + h, y[-1] + h * yk2, z[-1] + h * zk2)
        zk4 = f1(x + h, y[-1] + h * yk2, z[-1] + h * zk2)
        yc = y[-1] + (h / 6.0) * (yk1 + 2 * yk2 + 2 * yk3 + yk4)
        
        zc = z[-1] + (h / 6.0) * (zk1 + 2 * zk2 + 2 * zk3 + zk4)
        y.append(yc)
        z.append(zc)
        x += h
    return y

def task4():
    for cur_a in a:
        xs = []
        ys1 = []
        ys2 = []
        for n in range(1, 501):
            xs.append(n)
            h = T / n
            math_sol = [sol(cur_a)(i * h) for i in range(n + 1)] 
            solution1 = euler(lambda x, y, z: z, lambda x, y, z: cur_a ** 2 * y, 1, -cur_a, n)
            solution2 = runge(lambda x, y, z: z, lambda x, y, z: cur_a ** 2 * y, 1, -cur_a, n)
            ys1.append(np.log10(max([abs(solution1[i] - math_sol[i]) for i in range(len(math_sol))])))
            ys2.append(np.log10(max([abs(solution2[i] - math_sol[i]) for i in range(len(math_sol))])))
        graph.plot(xs, ys1, ys2)
        graph.ylabel("log10(err)")
        graph.xlabel("N")
        graph.legend(("euler", "runge"))
        graph.plot([1, 500], [-3, -3])
        graph.title("T = 10; a = " + str(cur_a))
        graph.show()
        #graph.savefig("hw7_t4_a" + ("1" if cur_a == 1 else "sqrt(20)") + ".png")
        #graph.close()

def task4b():
    for cur_a in a:
        xs = []
        ys1 = []
        ys2 = []
        ys3 = []
        for n in range(1, 501):
            xs.append(n)
            h = T / n
            math_sol = [sol(cur_a)(i * h) for i in range(n + 1)] 
            solution1 = euler(lambda x, y, z: z, lambda x, y, z: cur_a ** 2 * y, 1, -cur_a, n)
            solution2 = runge(lambda x, y, z: z, lambda x, y, z: cur_a ** 2 * y, 1, -cur_a, n)
            solution3 = runge4(lambda x, y, z: z, lambda x, y, z: cur_a ** 2 * y, 1, -cur_a, n)
            ys1.append(np.log10(max([abs(solution1[i] - math_sol[i]) for i in range(len(math_sol))])))
            ys2.append(np.log10(max([abs(solution2[i] - math_sol[i]) for i in range(len(math_sol))])))
            ys3.append(np.log10(max([abs(solution3[i] - math_sol[i]) for i in range(len(math_sol))])))
        graph.plot(xs, ys1)
        graph.plot(xs, ys2)
        graph.plot(xs, ys3)
        graph.ylabel("log10(err)")
        graph.xlabel("N")
        graph.legend(("euler", "runge", "runge4"))
        graph.plot([1, 500], [-3, -3])
        graph.title("T = 10; a = " + str(cur_a))
        graph.show()
        #graph.savefig("hw7_t4b_a" + str(int(cur_a * 1000)) + ".png")
        #graph.close()
        
if __name__ == "__main__":
    task4b()
