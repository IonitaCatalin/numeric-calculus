import math
import numpy as np


def horner(a, x):
    d = a[len(a)-1]

    i = len(a)-2
    while i >= 0:
        d = a[i] + d*x
        i = i-1

    return d


def polynom(x0, xn, n, x, function):
    values = []
    values.append((x0, function(x0)))
    diff = (xn - x0) / n

    for i in range(1, n):
        values.append((x0 + i*diff, function(x0 + i*diff)))

    values.append((xn, function(xn)))

    B = np.array([[values[i][0]**j for j in range(n+1)] for i in range(n+1)])
    f = np.array([values[i][1] for i in range(n+1)])

    a = np.dot(np.linalg.inv(B), f)


    print("Solutie folosind polinoame: ")
    sol = horner(a, x)
    print(sol)

    print("")

    print("Solutie reala: ")
    print(function(x))

    print("")
    print("Eroare: ")
    print(abs(sol - function(x)))



def compute_A(i, values, da):
    if(i == 0):
        return da
    else:
        return -compute_A(i-1, values, da) + (2*(values[i][1] - values[i-1][1])/(values[i][0] - values[i-1][0]))


def spline(x0, xn, n, x, function, derivate):
    da = derivate(x0)

    values = []
    values.append((x0, function(x0)))
    diff = (xn - x0) / n

    for i in range(1, n):
        values.append((x0 + i * diff, function(x0 + i * diff)))

    values.append((xn, function(xn)))

    i = 0
    while x > values[i][0]:
        i = i+1

    sol = ((compute_A(i, values, da) - compute_A(i-1, values, da))/(2*(values[i][0] - values[i-1][0])))*((x - values[i-1][0])**2) + compute_A(i-1, values, da) * (x - values[i-1][0]) + values[i-1][1]

    print("Solutie folosind spline: ")
    print(sol)

    print("")

    print("Solutie reala: ")
    print(function(x))

    print("")
    print("Eroare: ")
    print(abs(sol - function(x)))



def f1(x):
    return x**2 - 12*x + 30


def f1_drv(x):
    return x - 12


def f2(x):
    return math.sin(x) - math.cos(x)


def f2_drv(x):
    return math.cos(x) + math.sin(x)


def f3(x):
    return 2*(x**2) - 3*x + 15


def f3_drv(x):
    return 4*x - 3



if __name__ == '__main__':
    x0 = float(input("x_0: "))
    xn = float(input("x_n: "))
    n = int(input("n: "))

    x = float(input("x: "))

    if not x0 < x < xn:
        exit(0)

    print("\n\n Choose a function: ")
    print("1: x^2 -12x + 30")
    print("2: sin(x) - cos(x)")
    print("3: 2x^3 - 3x + 15")

    fun_case = int(input())

    function = None
    derivate = None

    if fun_case == 1:
        function = f1
        derivate = f1_drv
    elif fun_case == 2:
        function = f2
        derivate = f2_drv
    elif fun_case == 3:
        function = f3
        derivate = f3_drv
    else:
        exit(0)

    print("\n\n\n\n")
    polynom(x0, xn, n, x, function)

    print("\n\n\n\n")
    spline(x0, xn, n, x, function, derivate)






