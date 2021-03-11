from copy import deepcopy


def is_close(a, b, rel_tol=1e-15, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def check_colision(line, j):
    for i, element in enumerate(line):
        if element[1] == j:
            return i
    return -1


def read_rare_matrix(path):
    m = list()

    with open(path) as f:
        n = int(f.readline())
        m = [list() for i in range(n)]

        for line in f.readlines():
            line = line.split(", ")
            value = float(line[0])
            i = int(line[1])
            j = int(line[2])

            collision_index = check_colision(m[i], j)
            if collision_index >= 0:
                m[i][collision_index][0] = m[i][collision_index][0] + value
            else:
                m[i].append([value, j])

    return m


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


def check(a, b):
    for i in range(len(a)):
        if len(a[i]) != len(b[i]):
            print(f"Linia nu are acelasi numar de elemente nenule linia {i}")
            return False
        for element in a[i]:
            collision_index = check_colision(b[i], element[1])
            if collision_index >= 0:
                if not is_close(element[0], b[i][collision_index][0]):
                    print("Elementele nu sunt egale")
                    return False
    return True




def addition(m, d):
    add = deepcopy(m)

    for i, element in enumerate(d[0]):
        collision_index = check_colision(add[i], i)
        if collision_index >= 0:
            add[i][collision_index][0] = add[i][collision_index][0] + element
        else:
            add[i].append([element, i])

    for i, element in enumerate(d[1]):
        collision_index = check_colision(add[i], i + d[4])
        if collision_index >= 0:
            add[i][collision_index][0] = add[i][collision_index][0] + element
        else:
            add[i].append([element, i + d[4]])

    for i, element in enumerate(d[2]):
        index = i + 1
        collision_index = check_colision(add[index], index - d[3])
        if collision_index >= 0:
            add[index][collision_index][0] = add[index][collision_index][0] + element
        else:
            add[index].append([element, index - d[3]])

    return add



def multiplication(m, d):
    mult = [list() for i in range(len(m))]

    for i, line in enumerate(m):
        for j in range(len(m)):
            sum = 0
            for el in line:
                if j == el[1]:
                    sum += el[0]*d[0][el[1]]
                elif j == el[1]+d[4]:
                    sum += el[0]*d[1][el[1]]
                elif j == el[1]-d[3]:
                    sum += el[0]*d[2][el[1]-1]

            if sum > 0:
                mult[i].append([sum, j])

    return mult






if __name__ == '__main__':

    m = read_rare_matrix("a.txt")
    d = read_tridiagonal_matrix("b.txt")

    print(check(addition(m, d), read_rare_matrix("aplusb.txt")))
    print(check(multiplication(m, d), read_rare_matrix("aorib.txt")))

