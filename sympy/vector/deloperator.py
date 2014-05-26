#from sympy.vector import express
from sympy.core import Basic, Symbol
from sympy import diff
from vector import Vector
from vector import i, j, k
from scalar import x, y, z


class Del(Basic):
    """
    Represents the vector differential operator, usually represented in
    mathematical expressions as the 'nabla' symbol.
    """

    def __init__(self, system):
        self._system = system
        #self._x, self._y, self._z = self.system.base_scalars()
        #self._i, self._j, self._k = self.system.base_vectors()
        self._x, self._y, self._z = x, y, z
        self._i, self._j, self._k = i, j, k

    @property
    def system(self):
        return self._system

    def __call__(self, scalar_field):
        """
        Represents the gradient of the given scalar field with respect
        to the coordinate system this operator instance belongs to.

        Parameters
        ==========

        scalar_field : SymPy expression
            The scalar field to calculate the gradient of.

        Examples
        ========

        """

        #scalar_field = express(scalar_field, self.system, \
        #                       variables = True)
        vx = diff(scalar_field, self._x)
        vy = diff(scalar_field, self._y)
        vz = diff(scalar_field, self._z)

        return vx*self._i + vy*self._j + vz*self._k

    def dot(self, vect):
        """
        Represents the dot product between this operator and a given
        vector - equal to the divergence of the vector field in the
        coordinate system corresponding to this operator.

        Parameters
        ==========

        vect : Vector
            The vector whose divergence is to be calculated.

        Examples
        ========

        """

        #vect = express(vect, self._system, variables = True)
        vx = diff(vect.dot(self._i), self._x)
        vy = diff(vect.dot(self._j), self._y)
        vz = diff(vect.dot(self._k), self._z)

        return vx + vy + vz

    __and__ = dot

    def cross(self, vect):
        """
        Represents the cross product between this operator and a given
        vector - equal to the curl of the vector field in the
        coordinate system corresponding to this operator.

        Parameters
        ==========

        vect : Vector
            The vector whose curl is to be calculated.

        Examples
        ========

        """

        #vect = express(vect, self._system, variables = True)
        vectx = vect.dot(self._i)
        vecty = vect.dot(self._j)
        vectz = vect.dot(self._k)
        outvec = Vector.Zero
        outvec += (diff(vectz, self._y) - diff(vecty, self._z)) * self._i
        outvec += (diff(vectx, self._z) - diff(vectz, self._x)) * self._j
        outvec += (diff(vecty, self._x) - diff(vectx, self._y)) * self._k

        return outvec

    __xor__ = cross

delop = Del(Symbol('R'))
