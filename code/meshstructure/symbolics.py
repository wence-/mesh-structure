import numbers
import itertools
from pymbolic import primitives as pym
from .utils import lazyprop


class Index(pym.Variable):
    count = itertools.count()

    def __init__(self, lo, hi):
        super().__init__("i{}".format(next(self.count)))
        assert isinstance(lo, numbers.Integral) and isinstance(hi, numbers.Integral)
        # [lo, hi)
        self.lo = lo
        self.hi = hi

    def __str__(self):
        return self.name

    def __repr__(self):
        return "Index({}, lo={}, hi={})".format(self.name, self.lo, self.hi)

    @lazyprop
    def extent(self):
        return self.hi - self.lo
