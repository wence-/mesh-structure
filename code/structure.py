import islpy as isl
import itertools
from pymbolic import primitives as pym
from pymbolic.mapper.evaluator import evaluate
import abc
import numpy
import numbers
from functools import reduce, singledispatch
import operator
import ufl


class lazyprop(object):

    '''A read-only @property that is only evaluated once. The value is cached
    on the object itself rather than the function or class; this should prevent
    memory leakage.'''

    def __init__(self, fget, doc=None):
        self.fget = fget
        self.__doc__ = doc or fget.__doc__
        self.__name__ = fget.__name__
        self.__module__ = fget.__module__

    def __get__(self, obj, cls):
        if obj is None:
            return self
        obj.__dict__[self.__name__] = result = self.fget(obj)
        return result


class Index(pym.Variable):
    def __init__(self, name, lo, hi):
        super().__init__(name)
        assert isinstance(lo, numbers.Integral) and isinstance(hi, numbers.Integral)
        # [lo, hi)
        self.lo = lo
        self.hi = hi

    @lazyprop
    def extent(self):
        return self.hi - self.lo


def triangular_linear_index_map(i, j, n):
    """Given (i, j) : 0 <= i, j < n and i + j < n, produce a
    linear index."""
    return (n*(n-1)) // 2 - ((n - i)*(n - i - 1))//2 + i + j


def tetrahedral_linear_index_map(i, j, k, n):
    """Given (i, j, k): 0 <= i, j, k < n and i + j + k < n, produce a
    linear index."""
    ioff = (n*(n+1)*(n+2)) // 6 - ((n - i)*(n - i + 1)*(n - i + 2)) // 6
    return ioff + triangular_linear_index_map(j, k, n - i)


class EntitySet(metaclass=abc.ABCMeta):
    def __init__(self, indices, constraints, variant_tag=None):
        self.indices = tuple(indices)
        self.constraints = tuple(constraints)
        self.variant_tag = variant_tag

    @lazyprop
    def index_extents(self):
        """The extent of each index in the set."""
        return tuple(i.extent for i in self.indices)

    @abc.abstractmethod
    def linear_index_map(self, index_exprs, index_order):
        """Produce a map from indices in the set into a linear index.

        :arg index_exprs: Expressions for each index the map should apply to.
        :arg index_order: The order in which to apply the map to the indices.
        """

    @lazyprop
    def size(self):
        """The total number of points in the set."""
        n = 0
        for point in numpy.ndindex(*(i.extent for i in self.indices)):
            bindings = dict((i.name, i.lo + p) for (i, p) in zip(self.indices, point))
            if all(evaluate(constraint, bindings) for constraint in self.constraints):
                n += 1
        assert n == self.isl_set.count_val().get_num_si()
        return n

    @lazyprop
    def isl_set(self):
        v = isl.make_zero_and_vars(tuple(i.name for i in self.indices))
        exprs = []
        for index in self.indices:
            expr = ((index.lo + v[0]).le_set(v[index.name]) &
                    v[index.name].lt_set(index.hi + v[0]))
            exprs.append(expr)

        @singledispatch
        def translate(expr, v):
            raise AssertionError("Unhandled type %r" % type(expr))

        @translate.register(pym.Sum)
        def translate_sum(expr, v):
            return operator.add(*(translate(c, v) for c in expr.children))

        @translate.register(pym.Variable)
        def translate_variable(expr, v):
            return v[expr.name]

        @translate.register(numbers.Integral)
        def translate_number(expr, v):
            return v[0] + expr

        @translate.register(pym.Comparison)
        def translate_comparison(expr, v):
            left = translate(expr.left, v)
            right = translate(expr.right, v)
            fn = {">": "gt_set",
                  ">=": "ge_set",
                  "==": "eq_set",
                  "!=": "ne_set",
                  "<": "lt_set",
                  "<=": "le_set"}[expr.operator]
            return getattr(left, fn)(right)

        for constraint in self.constraints:
            expr = translate(constraint, v)
            exprs.append(expr)

        return reduce(operator.and_, exprs)


class SimplexEntitySet(EntitySet):
    def __init__(self, index_names, extent, variant_tag=None):
        assert isinstance(extent, numbers.Integral)
        self.extent = extent
        indices = tuple(Index(name, 0, extent) for name in index_names)
        constraints = (pym.Comparison(reduce(operator.add, indices), "<", extent), )
        super().__init__(indices, constraints, variant_tag=variant_tag)


class PointEntitySet(SimplexEntitySet):
    def linear_index_map(self, index_exprs, index_order):
        assert len(index_exprs) == 0
        return 0


class IntervalEntitySet(SimplexEntitySet):
    def linear_index_map(self, index_exprs, index_order):
        i, = index_exprs
        return i


class TriangleEntitySet(SimplexEntitySet):
    def linear_index_map(self, index_exprs, index_order):
        """index_order is ignored (just swap the order in index_exprs)"""
        i, j = index_exprs
        return triangular_linear_index_map(i, j, self.extent)


class TetrahedronEntitySet(SimplexEntitySet):
    def linear_index_map(self, index_exprs, index_order):
        """index_order is ignored (just swap the order in index_exprs)"""
        i, j, k = index_exprs
        return tetrahedral_linear_index_map(i, j, k, self.extent)


class TensorProductEntitySet(EntitySet):
    def __init__(self, *factors, variant_tag=None):
        self.factors = tuple(factors)
        indices = itertools.chain(*(f.indices for f in self.factors))
        assert len(set(i.name for i in indices)) == len(indices), "Must provide unique index names"
        constraints = itertools.chain(*(f.constraints for f in self.factors))
        super().__init__(indices, constraints, variant_tag=variant_tag)

    def linear_index_map(self, index_exprs, index_order):
        assert len(index_exprs) == sum(len(f.indices) for f in self.factors)
        assert len(index_order) == len(self.factors)
        factors = tuple(self.factors[i] for i in index_order)
        strides = tuple(numpy.cumprod((1, ) + tuple(f.size for f in reversed(factors)))[:-1][::-1])
        expr = 0
        for factor, stride in zip(factors, strides):
            nindex = len(factor.indices)
            index_expr = index_exprs[:nindex]
            index_exprs = index_exprs[nindex:]
            expr = expr + factor.linear_index_map(index_expr, None)*stride
        return expr


class StructureBase(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def entity_variants(self, codim=None):
        pass

    @abc.abstractmethod
    def subentity_map(self, entity_set, subentity_set, multiindex, subentity_multiindex):
        pass

    @abc.abstractmethod
    def dual_subentity_map(self, subentity_set, entity_set, multiindex):
        pass


class Extrusion(StructureBase):
    def __init__(self, nlevel):
        super().__init__()
        self.nlevel = nlevel
        # Extruded triangular prism
        cells = TensorProductEntitySet(TriangleEntitySet(("i", "j"), 1,
                                                         variant_tag=ufl.triangle),
                                       IntervalEntitySet(("k", ), nlevel,
                                                         variant_tag=ufl.interval),
                                       variant_tag=ufl.TensorProductCell(ufl.triangle, ufl.interval))
        hfaces = TensorProductEntitySet(TriangleEntitySet(("i", "j"), 1,
                                                         variant_tag=ufl.triangle),
                                       IntervalEntitySet(("k", ), nlevel+1,
                                                         variant_tag=ufl.interval),
                                       variant_tag=ufl.TensorProductCell(ufl.triangle, ufl.point))
        vfaces = TensorProductEntitySet(IntervalEntitySet
        self.entities = {(0, 0): cells,
                         (0, 1): hfaces,
                         (1, 0): vfaces,
                         (1, 1): hedges,
                         (2, 0): vedges,
                         (2, 1): vertices}

        def cellhface(mi, s):
            assert s in {0, 1}
            mi, = mi
            return (mi + s, )

        def cellvface(mi, s):
            assert s in {0, 1, 2}
            mi, = mi
            return (mi, s)

        def hfacehedge(mi, s):
            assert s in {0, 1, 2}
            mi, = mi
            return (mi, s)

        def vfacehedge(mi, s):
            assert s in {0, 1}
            i, j = mi
            return (i + s, j)

        def vfacevedge(mi, s):
            assert s in {0, 1}
            i, j = mi
            return (i, (j + s) % 3)

        def hedgevertex(mi, s):
            assert s in {0, 1}
            i, j = mi
            return (i, (j + s) % 3)

        def vedgevertex(mi, s):
            assert s in {0, 1}
            i, j = mi
            return (i + s, j)

        self.maps_ = {(cells, hfaces): cellhface,
                      (cells, vfaces): cellvface,
                      (hfaces, hedges): hfacehedge,
                      (vfaces, hedges): vfacehedge,
                      (vfaces, vedges): vfacevedge,
                      (hedges, vertices): hedgevertex,
                      (vedges, vertices): vedgevertex}

        def hfacecell(mi):
            i, = mi
            return ((i - 1, ), (i, ))

        def vfacecell(mi):
            i, _ = mi
            return ((i, ), )

        def hedgehface(mi):
            i, j = mi
            return ((i, ), )

        def hedgevface(mi):
            i, j = mi
            return ((i - 1, j), (i, j))

        def vedgevface(mi):
            i, j = mi
            return ((i, (j + 2) % 3), (i, j))

        def vertexhedge(mi):
            i, j = mi
            return ((i, (j + 2) % 3), (i, j))

        def vertexvedge(mi):
            i, j = mi
            return ((i - 1, j), (i, j))

        self.dual_maps_ = {(hfaces, cells): hfacecell,
                           (vfaces, cells): vfacecell,
                           (hedges, hfaces): hedgehface,
                           (hedges, vfaces): hedgevface,
                           (vedges, vfaces): vedgevface,
                           (vertices, hedges): vertexhedge,
                           (vertices, vedges): vertexvedge}

    def entity_variants(self, codim=None):
        if codim is None:
            return tuple(self.entities.values())
        try:
            b, e = codim
            return (self.entities[codim], )
        except TypeError:
            return tuple(v for k, v in self.entities.items() if sum(k) == codim)

    def subentity_map(self, eset, sset, multiindex, subentity):
        try:
            return self.maps_[(eset, sset)](multiindex, subentity)
        except KeyError:
            raise ValueError("No map between {} and {}".format(eset, sset))

    def dual_subentity_map(self, sset, eset, multiindex):
        try:
            return self.dual_maps_[(sset, eset)](multiindex)
        except KeyError:
            raise ValueError("No dual map between {} and {}".format(sset, eset))
