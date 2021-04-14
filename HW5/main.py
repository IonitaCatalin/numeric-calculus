from sklearn import datasets
import numpy as np
import math


def isclose(a, b, rel_tol=1e-05, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def sign(x):
    if x < 0:
        return -1
    else:
        return 1


def compute_theta(A, p, q):
    alpha = (A[p][p] - A[q][q])/(2*A[p][q])
    t = - alpha + sign(alpha)*math.sqrt(alpha**2 + 1)
    c = 1/math.sqrt(1 + t**2)
    s = t/math.sqrt(1 + t**2)

    return t, s, c



def compute_p_q(A):
    p, q = 0, 0
    max = - math.inf
    for i in range(1, len(A)):
        for j in range(i):
            if abs(A[i][j]) > max:
                max = abs(A[i][j])
                p = i
                q = j

    return p, q


def compute_rotation(A, p, q, t, s, c):
    R = [[0.0 for i in range(len(A))] for j in range(len(A))]

    for i in range(len(A)):
        for j in range(len(A)):
            if i == j and i != p and i != q:
                R[i][j] = 1
            elif i == j and (i == p or i == q):
                R[i][j] = c
            elif i == p and j == q:
                R[i][j] = s
            elif i == q and j == p:
                R[i][j] = -s
            else:
                R[i][j] = 0

    return np.array(R)


def compare(A_init_U, U_A, norm):
    return np.linalg.norm(A_init_U - U_A, norm)


def Jacobi(A):
    k = 0
    U = np.identity(len(A))
    p, q = compute_p_q(A)
    t, s, c = compute_theta(A, p, q)
    A_np = np.array(A)

    while not isclose(A_np[p][q], 0) and k < 1000:
        R = compute_rotation(A_np, p, q, t, s, c)
        A_np = np.dot(np.dot(R, A_np), R.T)
        U = np.dot(U, R.T)
        p, q = compute_p_q(A_np)
        t, s, c = compute_theta(A_np, p, q)
        k = k+1

    return compare(np.dot(A, U), np.dot(U, A_np), 2)


def diagonal(A):
    return [A[i][i] for i in range(len(A))]



def lu_cpy(A):
    n = A.shape[0]
    U = A.copy()
    L = np.eye(n, dtype=np.double)

    for i in range(n):
        factor = U[i+1:, i] / U[i, i]
        L[i+1:, i] = factor
        U[i+1:] -= factor[:, np.newaxis] * U[i]

    return L, U


def Cholenski(A):
    A_k = np.array([row[:] for row in A])
    L, U = lu_cpy(A)
    A_j = np.dot(U, L)

    k = 0

    while not isclose(compare(A_k, A_j, 2), 0) and k < 1000:
        A_k = np.array([row[:] for row in A_j])
        L, U = lu_cpy(A_k)
        A_j = np.dot(U, L)

        k = k + 1

    return A_j


if __name__ == '__main__':
    A = datasets.make_spd_matrix(4)

    print("Eroare Jacobi: ")
    print(Jacobi(A), end="\n \n \n")

    print("Descompunerea Cholenski: ")
    print(Cholenski(A), end = "\n \n \n")

    A_2 = np.random.rand(8, 4)
    svd = np.linalg.svd(A_2, full_matrices=False)
    print("Matricea folosita in partea a doua")
    print(A_2, end="\n \n \n")

    print("Valorile singulare: ")
    print(svd[1], end="\n \n \n")

    print("Rank-ul matricei: ")
    print(len(svd[1]), end="\n \n \n")

    print("Numarul de conditionare: ")
    print(max(svd[1]) / min([i for i in svd[1] if i > 0]), end="\n \n \n")

    print("pseudo inversa Moore-Penrose")
    A_I = np.dot(np.dot(svd[2].T, np.linalg.inv(np.diag(svd[1]))), svd[0].T)
    print(A_I, end="\n \n \n")

    print("Pseudo inversa in sensul celor mai mici patrate")
    A_J = np.dot(np.linalg.inv(np.dot(A_2.T, A_2)), A_2.T)
    print(A_J, end="\n \n \n")

    print("Norma")
    print(compare(A_I, A_J, 1))
