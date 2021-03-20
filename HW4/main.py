from math import sqrt

def read_f(path):
    x = list()

    with open(path) as f:
        n = int(f.readline())

        for i in range(n):
            x.append(float(f.readline()))

    return x


def read_tridiagonal_matrix(path):
    a = list()
    b = list()
    c = list()
    p = 0
    q = 0

    with open(path) as f:
        n = int(f.readline())
        p = int(f.readline())
        q = int(f.readline())

        for i in range(n):
            a.append(float(f.readline()))
        for i in range(n - q):
            b.append(float(f.readline()))
        for i in range(n - p):
            c.append(float(f.readline()))

    return a, b, c, p, q


def check_if_not_zero(a):
    for el in a[0]:
        if el == 0:
            return False
    for el in a[1]:
        if el == 0:
            return False
    for el in a[2]:
        if el == 0:
            return False

    return True




def gauss_seidel(path_a, path_f):
    a = read_tridiagonal_matrix(path_a)
    f = read_f(path_f)

    x_gs = [0.0 for i in range(len(f))]
    k = 0
    k_max = 1000
    eps = pow(10, -a[3])
    delta_x = eps

    while eps <= delta_x <= pow(10, 8) and k <= k_max:
        delta_x = 0
        for i in range(len(f)):
            sum_1 = 0
            sum_2 = 0
            temp = x_gs[i]
            if i > a[4] - 1:
                sum_1 = a[2][i-a[4]]*x_gs[i]
            if i < len(f) - a[3]:
                sum_2 = a[1][i + a[3] - 1]*x_gs[i]

            x_gs[i] = (f[i] - sum_1 - sum_2)/a[0][i]

            delta_x += pow(x_gs[i] - temp, 2)

        k += 1
        delta_x = sqrt(delta_x)

    return precision(a, x_gs, f)


def precision(a, x, f):
    sol = [0.0 for i in range(len(f))]

    sol[0] = a[0][0]*x[0] + a[1][0] * x[1]
    sol[len(f) - 1] = a[2][len(f) - 2]*x[len(f) - 2] + a[0][len(f) - 1] * x[len(f) - 1]
    for i in range(1, len(f) - 1):
            sol[i] = a[2][i-1] * x[i-1] + a[0][i] * x[i] + a[1][i] * x[i+1]

    print_solution_into_file("sol.txt", sol)
    inf_norm = max([abs(sol[i] - f[i]) for i in range(len(f))])

    return inf_norm





def print_solution_into_file(path, solution):
    with open(path, 'w+') as f:
        for i in range(len(solution)):
            f.write(str(solution[i]) + "\n")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print(gauss_seidel("a_1.txt", "f_1.txt"))