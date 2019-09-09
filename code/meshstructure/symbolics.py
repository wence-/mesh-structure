import itertools
import numbers

from pymbolic import primitives as pym

from .utils import lazyattr


class Index(pym.Variable):
    count = itertools.count()

    def __init__(self, lo, hi):
        """An index representing the range [lo, hi)"""
        super().__init__("i{}".format(next(self.count)))
        assert isinstance(lo, numbers.Integral) and isinstance(hi, numbers.Integral)
        # [lo, hi)
        self.lo = lo
        self.hi = hi

    def __str__(self):
        return self.name

    def __repr__(self):
        return "Index({}, lo={}, hi={})".format(self.name, self.lo, self.hi)

    @lazyattr
    def extent(self):
        return self.hi - self.lo


class Point(tuple):
    """A point in an entity set.

    :arg multiindex: a tuple of index expressions to index the entity set.
    :arg entity_set: The set being indexed."""
    def __new__(cls, multiindex, entity_set):
        multiindex = tuple(multiindex)
        assert len(multiindex) == len(entity_set.indices)
        return super().__new__(cls, (multiindex, entity_set))

    @lazyattr
    def multiindex(self):
        return self[0]

    @lazyattr
    def entity_set(self):
        return self[1]

    def __str__(self):
        return "Point(%s, %s)" % self

    def __repr__(self):
        return "Point(%r, %r)" % self
