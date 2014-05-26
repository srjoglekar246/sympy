if __name__ == '__main__':
    from sympy import *
    from vector import i, j, k
    from scalar import x, y, z
    from deloperator import delop
    a = Symbol('a')
    v = a * i + 4 * j - k
    q = Symbol('q')
    e = i*(sin(q) + cos(q))**2 + j
    u1, u2, u3 = symbols('u1 u2 u3')
    v1, v2, v3 = symbols('v1 v2 v3')
    vect1 = x*i + y**2*j + 3*x*z*k
    vect2 = 4*y*x*z*i + z**2*y*3*j + x*z*k
