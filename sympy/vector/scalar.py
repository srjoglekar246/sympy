from sympy.core.symbol import Dummy


class BaseScalar(Dummy):
    """
    A coordinate symbol/base scalar associated with a coordinate
    system.

    Ideally, users should not instantiate this class. Instances of
    this class must only be accessed through the corresponding system,
    for e.g. R.x for the X-coordinate symbol of system R.

    """

    def __new__(cls, system, index, name):
        obj = super(BaseScalar, cls).__new__(cls, name)
        if not isinstance(system, Symbol):
            raise TypeError(str(system) + " is not a valid system.")
        if index not in range(0, 3):
            raise ValueError("Invalid index specified.")
        obj._id = (system, index)
        obj._name = name

        return obj

    @property
    def system(self):
        return self._id[0]

    def __eq__(self, other):
        #Check if the other object is a BaseScalar of the same frame
        #and same index
        if isinstance(other, BaseScalar):
            if other._id == self._id:
                return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return tuple((self._id[0].__hash__(), self._id[1])).__hash__()

    def __str__(self, printer=None):
        return self._name

    __repr__ = __str__
    _sympystr = __str__


#Just some hacks for now
from sympy.core.symbol import Symbol
x = BaseScalar(Symbol('R'), 0, 'R.x')
y = BaseScalar(Symbol('R'), 1, 'R.y')
z = BaseScalar(Symbol('R'), 2, 'R.z')
