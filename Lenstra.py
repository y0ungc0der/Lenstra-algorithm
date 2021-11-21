import time
from Crypto.Util.number import *
from Crypto.Random import *
from math import log, floor
from sympy import nextprime

# ------- Lenstra's method ------- #

def main(n, m):
    it = 0
    while True:
        # select a random elliptic curve E(Z⁄nZ) and a point Q on it
        temp = random_curve(n)
        if len(temp) == 1:
            print("Iterations: {}".format(it))
            return temp[0]

        a = temp[0]
        b = temp[1]
        x = temp[2]
        y = temp[3]
        i = 0
        Q = [x, y]
        print(f'Elliptic curve defined by y^2 = x^3 + {a} * x + {b} over ring of integers modulo {n}')
        print(f'Point Q: ({x} : {y})')

        while i <= m:
            if not (isPrime(i)):
                i = nextprime(i)
                continue

            # α = largest integer not greater than 0,5 * (log(n)/log(p_i))
            alpha = floor((log(n, 2) / log(i, 2)) * (1 / 2))

            if alpha == 0:
                alpha = 1

            # multiply the point Q_i sequentially by the current prime number in degree α
            p = i ** alpha
            # calculate Q_i = p_i * Q_i
            for j in range(1, p):
                it += 1
                # Р = (x, y) и Q = (x_1, y_1) - calculate Q + P = (x_2, y_2)

                if Q[0] == x and Q[1] == y: # P = Q
                    # λ = (3x^(2) + a)/2y
                    _lambda = ((3 * Q[0] ** 2 + a) % n) * inverse((2 * Q[1]), n) % n

                    d = GCD(_lambda, n)
                    # for 1 < d < n result: d
                    if  d > 1 and d < n:
                        print("Iterations: {}".format(it))
                        return d

                    # x_2 = λ^(2) - 2* x_1 (mod n)
                    x_2 = (_lambda ** 2 - 2 * Q[0]) % n
                    # y_2 = λ*(x - x_2) - y_1 (mod n)
                    y_2 = (_lambda * (Q[0] - x_2) - Q[1]) % n
                    Q[1] = y_2
                    Q[0] = x_2

                else:
                    # λ = (y - y_1)/(x - x_1)
                    _lambda = ((y - Q[1]) % n) * inverse((x - Q[0]), n) % n

                    d = GCD(_lambda, n)
                    # for 1 < d < n result: d
                    if  d > 1 and d < n:
                        print("Iterations: {}".format(it))
                        return d

                    # x_2 = λ^(2) - 2 * x_1 (mod n)
                    x_2 = (_lambda ** 2 - Q[0] - x) % n
                    # y_2 = λ*(x_1 - x_2) - y (mod n)
                    y_2 = (_lambda * (Q[0] - x_2) - Q[1]) % n
                    Q[1] = y_2
                    Q[0] = x_2

def random_curve(n):
    while True:
        # select a random value of x, y, a from Z_n.  Weierstrass elliptic curve form is y^(2) = x^(3) + ax + b
        a = getRandomRange(0, n - 1, randfunc = get_random_bytes)
        x = getRandomRange(0, n - 1, randfunc = get_random_bytes)
        y = getRandomRange(0, n - 1, randfunc = get_random_bytes)
        # calculate b 
        b = (y ** 2 - x ** 3 - a * x) % n
 
        c = 4 * a ** 3 + 27 * b ** 2
        res = []
        g = GCD(n, c)

        if g == n:
            continue

        if g > 1:
            res.append(g)
            return (res)

        res.append(a)
        res.append(b)
        res.append(x)
        res.append(y)
        return (res)

if __name__ == '__main__':
    sbase = int(input('Select size of base: '))
    bits = int(input('Select length of n (in bits): '))
    
    start_time = time.time()
    
    # form a composite number n 
    p_get = getPrime(bits//2, randfunc = get_random_bytes)
    q_get = getPrime(bits//2, randfunc = get_random_bytes)
    n = p_get * q_get
    
    # calculate divisor of number n
    p = main(n, sbase)
    
    print("{} * {} = {}".format(p, n // p, n))
    print(f'\nRuntime = {time.time() - start_time} seconds\n')