from sympy.core.symbol import Dummy, Symbol


class BaseScalar(Dummy):
    """
    A coordinate symbol/base scalar.

    Ideally, users should not instantiate this class.

    """

    def __new__(cls, name, index):
        obj = super(BaseScalar, cls).__new__(cls, name)
        if index not in range(0, 3):
            raise ValueError("Invalid index specified.")
        obj._id = (index, Symbol('R'))
        obj._name = name

        return obj

    def __eq__(self, other):
        #Check if the other object is a BaseScalar of same index
        if isinstance(other, BaseScalar):
            if other._id == self._id:
                return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return tuple(self._id).__hash__()

    def __str__(self, printer=None):
        return self._name

    __repr__ = __str__
    _sympystr = __str__


#Just some hacks for now
x = BaseScalar('x', 0)
y = BaseScalar('y', 1)
z = BaseScalar('z', 2)