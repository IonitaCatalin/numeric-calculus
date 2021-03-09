import numpy as np
from scipy.linalg import lu
import math
import sys
import random
from sklearn import datasets
from cmath import sqrt

A = list()
b = list()


def convert_to_L():
    global A

    for p in range(len(A)):
        if A[p][p] - sum([A[p][j] ** 2 for j in range(p)]) < 0:
            print(A[p][p] - sum([A[p][j] ** 2 for j in range(p)]))
            return False

        A[p][p] = math.sqrt(A[p][p] - sum([A[p][j] ** 2 for j in range(p)]))

        if p < (len(A) - 1):
            for i in range(p + 1, len(A)):
                A[i][p] = (A[i][p] - sum([A[i][j] * A[p][j] for j in range(p)])) / A[p][p]

    return True


def getDet(A):
    return np.prod([A[i][i] for i in range(len(A))]) ** 2


def getY(A, b):
    y = list()
    for i in range(len(A)):
        y.append((b[i] - sum([A[i][j] * y[j] for j in range(0, i)])) / A[i][i])

    return y


def getX(A, y):
    x = [0.0 for i in range(len(A))]

    for i in range(len(A) - 1, -1, -1):
        x[i] = ((y[i] - sum([A[j][i] * x[j] for j in range(len(A) - 1, i, -1)])) / A[i][i])

    return x


def validate_solution(A, x, b, diag):
    y = [0.0 for i in range(len(A))]

    for i in range(len(A)):
        for j in range(len(A[i])):
            if j > i:
                y[i] += A[i][j] * x[i]
            elif i > j:
                y[i] += A[j][i] * x[i]
            else:
                y[i] += diag[i] * x[i]

    return sqrt(sum([y[i] - b[i] for i in range(len(y))])).real


def validate_inv(A_chol, A_bibl):
    return np.linalg.norm(A_chol - A_bibl)


def lu_cpy(A):
    n = A.shape[0]
    U = A.copy()
    L = np.eye(n, dtype=np.double)

    for i in range(n):
        factor = U[i + 1:, i] / U[i, i]
        L[i + 1:, i] = factor
        U[i + 1:] -= factor[:, np.newaxis] * U[i]

    return L, U


def compute_inverse(A):
    A_inverse = list()
    for i in range(len(A)):
        e = [0 for i in range(len(A))]
        e[i] = 1
        x = getX(A, getY(A, e))
        A_inverse.append(x)

    return np.array(A_inverse).T


def isclose(a, b, rel_tol=1e-05, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def check_simetry(A):
    for i in range(len(A)):
        for j in range(i, len(A[i])):
            if not isclose(A[i][j], A[j][i]):
                print(A[i][j], A[j][i])
                return False

    return True


def get_A_init(A, diag):
    A_init = [[0.0 for i in range(len(A))] for j in range(len(A))]

    for i in range(len(A)):
        for j in range(len(A[i])):
            if j > i:
                A_init[i][j] = A[i][j]
            elif i > j:
                A_init[i][j] = A[j][i]
            else:
                A_init[i][j] = diag[i]

    return A_init


if __name__ == '__main__':

    '''
        python main.py [type] [aux]
            [type] 1 => file reading [aux] = file name / default: input.txt
            [type] 2 => random generating [aux] = matrix size / default: 100
            [type] 3 => keyboard input [aux] = matrix size / default: 3
    '''

    if sys.argv[1] == '1':

        file_name = 'input.txt'
        if len(sys.argv) > 1:
            file_name = sys.argv[2]

        with open(file_name) as f:
            b = [float(n) for n in f.readline().strip().split(' ')]
            A = [[float(n) for n in line.strip().split(' ')] for line in f]
    elif sys.argv[1] == '2':

        count = 100
        if len(sys.argv) > 1:
            count = int(sys.argv[2])

        A = [[0.0 for j in range(0, int(count))] for i in range(0, int(count))]
        b = [0.0 for i in range(0, int(count))]

        for i in range(count):
            x = random.uniform(0.0, 100.0)
            b[i] = x
            for j in range(count):
                x = random.uniform(0.0, 100.0)
                A[i][j] = x

            A = np.dot(np.array(A), np.array(A).T)

        A = datasets.make_spd_matrix(count)

    elif sys.argv[1] == '3':
        count = 3
        if len(sys.argv) > 1:
            count = int(sys.argv[2])

        print(f"Introduceti matricea A pe linii, o linie = {count} inputs")
        A = [[int(input()) for j in range(0, int(count))] for i in range(0, int(count))]

        print(f"Introduceti vectorul B: {count} inputs")
        b = [int(input()) for i in range(0, int(count))]

    if check_simetry(A):
        diag = [A[i][i] for i in range(len(A))]

        if convert_to_L():

            print("Determinant:")
            print(getDet(A))

            print(" ")

            print("Solutia X pentru sistemul A * X = B, unde A = L * L.T este:")
            print(getX(A, getY(A, b)))

            print(" ")

            print("Eroare la calcul A_chol * X = B:")
            print(validate_solution(A, getX(A, getY(A, b)), b, diag))

            print(" ")

            (L, U) = lu_cpy(np.array(get_A_init(A, diag)))
            print("Descompunerea LU:")
            print("L:")
            print(L)
            print("U:")
            print(U)

            print(" ")

            print("Solutia X pentru sistemul A*X = b, unde A = L * U este:")
            print(getX(U.T, getY(L, b)))

            print(" ")

            print("Eroare la calcul A^(-1)")
            print(validate_inv(compute_inverse(A), np.linalg.inv(np.array(get_A_init(A, diag)))))
        else:
            print("Matricea nu e pozitiv definita")
    else:
        print("Matricea nu e simetrica")
