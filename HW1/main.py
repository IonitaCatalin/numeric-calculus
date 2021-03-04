import math
import random as rand
import time


def precision():
    m = 1
    while 1.0 + pow(10, -m) != 1.0:
        m += 1
    return pow(10, -(m - 1))


def check_addition(u):
    return (1.0 + u / 10) + u / 10 == 1.0 + (u / 10 + u / 10)


def check_multiplication(u):
    return ((1.0 / u) * u ) * u  == (1.0 / u) * ((u / 10) * (u / 10))


def convert_with_period(x):
    while x > math.pi / 2:
        x = x - math.pi

    while x < -math.pi / 2:
        x = x + math.pi

    return x


def tan_poly(x):
    negativity = False
    external = False

    x = convert_with_period(x)

    if x > math.pi / 4:
        x = math.pi / 2 - x
        external = True
    elif x < math.pi / 4:
        x = math.pi / 2 + x
        external = True
        negativity = True

    c1 = 0.33333333333333333
    c2 = 0.133333333333333333
    c3 = 0.053968253968254
    c4 = 0.0218694885361552

    result = x + c1 * pow(x, 3) + c2 * pow(x, 5) + c3 * pow(x, 7) + c4 * pow(x, 9)

    if external:
        if negativity:
            return -1 / result
        else:
            return 1 / result
    else:
        if negativity:
            return -result
        else:
            return result


def tan_lentz(x, eps):
    negativity = False

    x = convert_with_period(x)

    if x < 0:
        x = -x
        negativity = True

    f = 0
    mic = pow(10, -12)
    if f == 0:
        f = mic
    C = f
    D = 0
    delta = 0
    b = 1
    a = x

    while abs(delta - 1) > eps:
        D = b + a * D
        if D == 0:
            D = mic

        C = b + a / C
        if C == 0:
            C = mic

        D = 1 / D
        delta = C * D
        f = delta * f

        a = -(pow(x, 2))
        b = b + 2

    if negativity:
        return -f
    else:
        return f


if __name__ == '__main__':
    u = precision()
    print(u)
    print(check_addition(u))
    print(check_multiplication(u))

    checksum_lentz = 0
    checksum_poly = 0

    total_time_lentz = 0
    total_time_poly = 0

    for i in range(1, 10_000):
        x = rand.uniform(-math.pi / 2, math.pi / 2)
        start = time.perf_counter()
        my_tan = tan_lentz(x, u)
        finish = time.perf_counter()

        total_time_lentz = finish - start
        checksum_lentz = checksum_lentz + abs(my_tan - math.tan(x))

    for i in range(1, 10_000):
        x = rand.uniform(-math.pi / 2, math.pi / 2)
        start = time.perf_counter()
        my_tan = tan_poly(x)
        finish = time.perf_counter()

        total_time_poly = finish - start
        checksum_poly = checksum_poly + abs(my_tan - math.tan(x))

    print(f'Lentz: \n Error: {checksum_lentz / 10_000} \n Time: {total_time_lentz}')
    print(f'Poly: \n Error: {checksum_poly / 10_000} \n Time: {total_time_poly}')
