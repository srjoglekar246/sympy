from sympy.simplify import simplify, trigsimp
from sympy import pi, sqrt, symbols, ImmutableMatrix as Matrix, \
     sin, cos
from sympy.vector.vector import i, j, k, Vector, BaseVector, VectorAdd, \
     VectorMul, VectorZero

a, b, c = symbols('a b c')


def test_vector():
    assert isinstance(i, BaseVector)
    assert i != j
    assert j != k
    assert k != i
    assert i - i == Vector.Zero
    assert i + Vector.Zero == i
    assert Vector.Zero != 0

    v1 = a*i + b*j + c*k
    v2 = a**2*i + b**2*j + c**2*k
    v3 = v1 + v2
    v4 = 2 * v1
    v5 = a * i

    assert isinstance(v1, VectorAdd)
    assert v1 - v1 == Vector.Zero
    assert v1 + Vector.Zero == v1
    assert v1.dot(i) == a
    assert v1.dot(j) == b
    assert v1.dot(k) == c
    assert i.dot(v2) == a**2
    assert j.dot(v2) == b**2
    assert k.dot(v2) == c**2
    assert v3.dot(i) == a**2 + a
    assert v3.dot(j) == b**2 + b
    assert v3.dot(k) == c**2 + c

    assert v1 + v2 == v2 + v1
    assert v1 - v2 == -1 * (v2 - v1)
    assert a * v1 == v1 * a

    assert isinstance(v5, VectorMul)
    assert isinstance(v4, Vector)
    assert isinstance(v4, VectorAdd)
    assert isinstance(v4, Vector)
    assert isinstance(Vector.Zero, VectorZero)
    assert isinstance(Vector.Zero, Vector)
    assert isinstance(v1 * 0, VectorZero)

    assert v1.to_matrix() == Matrix([[a], [b], [c]])

    assert i.components == {i: 1}
    assert v5.components == {i: a}
    assert v1.components == {i: a, j: b, k: c}

    assert VectorAdd(v1, Vector.Zero) == v1
    assert VectorMul(a, v1) == v1*a
    assert VectorMul(1, i) == i
    assert VectorAdd(v1, Vector.Zero) == v1
    assert VectorMul(0, Vector.Zero) == Vector.Zero


def test_vector_magnitude_normalize():
    assert Vector.Zero.magnitude() == 0
    assert Vector.Zero.normalize() == Vector.Zero

    assert i.magnitude() == 1
    assert j.magnitude() == 1
    assert k.magnitude() == 1
    assert i.normalize() == i
    assert j.normalize() == j
    assert k.normalize() == k

    v1 = a * i
    assert v1.normalize() == (a/sqrt(a**2))*i
    assert v1.magnitude() == sqrt(a**2)

    v2 = a*i + b*j + c*k
    assert v2.magnitude() == sqrt(a**2 + b**2 + c**2)
    assert v2.normalize() == v2 / v2.magnitude()


def test_vector_simplify():
    A, s, k, m = symbols('A, s, k, m')

    test1 = (1 / a + 1 / b) * i
    assert (test1 & i) != (a + b) / (a * b)
    test1 = simplify(test1)
    assert (test1 & i) == (a + b) / (a * b)
    assert test1.simplify() == simplify(test1)

    test2 = (A**2 * s**4 / (4 * pi * k * m**3)) * i
    test2 = simplify(test2)
    assert (test2 & i) == (A**2 * s**4 / (4 * pi * k * m**3))

    test3 = ((4 + 4 * a - 2 * (2 + 2 * a)) / (2 + 2 * a)) * i
    test3 = simplify(test3)
    assert (test3 & i) == 0

    test4 = ((-4 * a * b**2 - 2 * b**3 - 2 * a**2 * b) / (a + b)**2) * i
    test4 = simplify(test4)
    assert (test4 & i) == -2 * b

    v = (sin(a)+cos(a))**2*i - j
    assert trigsimp(v) == (2*sin(a + pi/4)**2)*i + (-1)*j
    assert trigsimp(v) == v.trigsimp()

    assert simplify(Vector.Zero) == Vector.Zero


def test_vector_dot():
    assert i.dot(Vector.Zero) == 0
    assert Vector.Zero.dot(i) == 0
    assert i & Vector.Zero == 0

    assert i.dot(i) == 1
    assert i.dot(j) == 0
    assert i.dot(k) == 0
    assert i & i == 1
    assert i & j == 0
    assert i & k == 0

    assert j.dot(i) == 0
    assert j.dot(j) == 1
    assert j.dot(k) == 0
    assert j & i == 0
    assert j & j == 1
    assert j & k == 0

    assert k.dot(i) == 0
    assert k.dot(j) == 0
    assert k.dot(k) == 1
    assert k & i == 0
    assert k & j == 0
    assert k & k == 1


def test_vector_cross():
    assert i.cross(Vector.Zero) == Vector.Zero
    assert Vector.Zero.cross(i) == Vector.Zero

    assert i.cross(i) == Vector.Zero
    assert i.cross(j) == k
    assert i.cross(k) == -j
    assert i ^ i == Vector.Zero
    assert i ^ j == k
    assert i ^ k == -j
    
    assert j.cross(i) == -k
    assert j.cross(j) == Vector.Zero
    assert j.cross(k) == i
    assert j ^ i == -k
    assert j ^ j == Vector.Zero
    assert j ^ k == i

    assert k.cross(i) == j
    assert k.cross(j) == -i
    assert k.cross(k) == Vector.Zero
    assert k ^ i == j
    assert k ^ j == -i
    assert k ^ k == Vector.Zero