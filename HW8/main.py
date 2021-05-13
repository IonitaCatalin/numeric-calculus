import math


def isclose(a, b, rel_tol=1e-10, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def f1(x):
    return 1.0/3.0*(x**3) - 2*(x**2) + 2*x + 3


def g1(f, x):
    h = 1e-05
    return (3*f(x) - 4*f(x-h) + f(x-2*h))/(2*h)


def g2(f, x):
    h = 1e-05
    return (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h))/(12*h)


def f1_d2(f, x):
    h = 1e-05
    return (-f(x+2*h) + 16*f(x+h) - 30*f(x) + 16*f(x-h) - f(x-2*h))/(12*(h**2))


def f2(x):
    return x**2 + math.sin(x)


def dh(f, g, eps=1e-10):
    x = 0
    k = 0
    k_max = 1000
    delta = eps
    max_value = pow(10, 8)


    while eps <= delta <= max_value and k < k_max:
        if abs(g(f, x + g(f, x)) - g(f, x)) < eps:
            return (x,k)
        z = x + (g(f, x) ** 2) / (g(f, x + g(f, x)) - g(f, x))
        temp = x - (g(f, x) * (g(f, z) - g(f, x)))/(g(f, x + g(f, x)) - g(f, x))
        delta = abs(temp - x)

        x = temp
        k += 1

    if delta < eps:
        return (x,k)
    else:
        return (None,k)


if __name__ == '__main__':
    value1,k1 = (dh(f1, g1))
    value2,k2 = (dh(f1, g2))

    print("Rezultat cu F'(x) = G1(x): ")
    print(value1)
    print("G1(rez):")
    print(g1(f1, value1))
    print("Iteratii necesare:")
    print(k1)

    print("")

    print("Rezultat cu F'(x) = G2(x): ")
    print(value2)
    print("G2(rez):")
    print(g2(f1, value2))
    print("Iteratii necesare:")
    print(k2)

    print("")

    print('Verificare F"(x) > 0')
    print(f1_d2(f1, value1))
