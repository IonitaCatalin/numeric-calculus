import numpy as np


def isclose(a, b, rel_tol=1e-07, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def contains_isClose(roots, value):
    for r in roots:
        if isclose(r, value):
            return True
    return False


def horner(a, x):
    d = a[0]
    for a_i in a[1:]:
        d = a_i + d*x

    return d


def compute_c(a, x):
    a_d1 = np.polyder(a)
    a_d2 = np.polyder(a, 2)

    return ((horner(a, x)**2) * horner(a_d2, x))/(horner(a_d1, x)**3)


def olver(a, value, eps=1e-10):
    x = value
    k = 0
    delta = eps

    a_d1 = np.polyder(a)

    max_val = pow(10, 8)
    k_max = 1000

    while eps <= delta <= max_val and k < k_max:
        temp = x - (horner(a, x) / horner(a_d1, x)) - 0.5*compute_c(a, x)
        delta = abs(temp - x)
        x = temp
        k += 1

    return x

if __name__ == '__main__':
    a = [1.0, -6.0, 11.0, -6.0]
    R = (abs(a[0]) + max([abs(a[i]) for i in range(1, len(a))]))/abs(a[0])

    step = 0.01
    poly_roots = []

    value = -R
    while value < R:
        temp = olver(a, value)
        if not contains_isClose(poly_roots, temp):
            poly_roots.append(temp)
        value += step

    with open("results.txt", "w+") as file:
        for value in poly_roots:
            file.write(str(value) + "\n")

    print(poly_roots)




