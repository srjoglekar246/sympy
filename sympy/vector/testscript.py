if __name__ == '__main__':
    from sympy import *
    from vector import i, j, k
    from scalar import x, y, z
    a = Symbol('a')
    v = a * i + 4 * j - k
    q = Symbol('q')
    e = i*(sin(q) + cos(q))**2 + j
    u1, u2, u3 = symbols('u1 u2 u3')
    v1, v2, v3 = symbols('v1 v2 v3')
    vect1 = u1*i + u2*j + u3*k
    vect2 = v1*i + v2*j + v3*k
