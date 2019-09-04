import abc
import itertools
import numbers
import operator
from functools import reduce, singledispatch

import islpy as isl
import numpy
from pymbolic import primitives as pym
from pymbolic.mapper.evaluator import evaluate

from .symbolics import Index
from .utils import lazyattr


__all__ = ("PointEntitySet", "IntervalEntitySet", "TriangleEntitySet", "TetrahedronEntitySet",
           "TensorProductEntitySet")


class EntitySet(metaclass=abc.ABCMeta):
    """A representation of some set of entities.

    :arg indices: :class:`Index` objects encoding points in the set.
    :arg constraints: constraints on the indices (expressions), must
        be expressible in Presburger arithmetic.
    :arg The codimension of the entities this entityset represents
    :arg variant_tag: Arbitrary data used to distinguish this set."""
    def __init__(self, indices, constraints, codimension, variant_tag=None):
        self.indices = tuple(indices)
        self.constraints = tuple(constraints)
        self.variant_tag = variant_tag
        self.codimension = codimension

    @lazyattr
    def index_extents(self):
        """The extent of each index in the set."""
        return tuple(i.extent for i in self.indices)

    @abc.abstractmethod
    def linear_index_map(self, index_exprs, index_order):
        """A map from indices in the set into a linear index.

        :arg index_exprs: Expressions for each index the map should apply to.
        :arg index_order: The order in which to apply the map to the indices.
        """

    @lazyattr
    def size(self):
        """The total number of points in the set."""
        n = 0
        for point in numpy.ndindex(*(i.extent for i in self.indices)):
            bindings = dict((i.name, i.lo + p) for (i, p) in zip(self.indices, point))
            if all(evaluate(constraint, bindings) for constraint in self.constraints):
                n += 1
        assert n == self.isl_set.count_val().get_num_si()
        return n

    @lazyattr
    def isl_set(self):
        """An ISLPY representation of the index set for this entity set."""
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

        if len(exprs) == 0:
            # Hack!
            return v[0].eq_set(v[0])
        return reduce(operator.and_, exprs)

    def __str__(self):
        return "{}({})".format(type(self).__name__, self.isl_set)

    __repr__ = __str__


def triangular_linear_index_map(i, j, n):
    """Given (i, j) : 0 <= i, j < n and i + j < n, produce a
    linear index."""
    return (n*(n-1)) // 2 - ((n - i)*(n - i - 1))//2 + i + j


def tetrahedral_linear_index_map(i, j, k, n):
    """Given (i, j, k): 0 <= i, j, k < n and i + j + k < n, produce a
    linear index."""
    ioff = (n*(n+1)*(n+2)) // 6 - ((n - i)*(n - i + 1)*(n - i + 2)) // 6
    return ioff + triangular_linear_index_map(j, k, n - i)


class SimplexEntitySet(EntitySet):
    """An entity set describing some number of simplices.

    :arg extent: Max extent in each direction.

    Produces an index set with

    0 <= indices < extent, sum(indices) < extent

    Where len(indices) == dimension
    """
    def __init__(self, extent, codimension, variant_tag=None):
        assert isinstance(extent, numbers.Integral)
        self.extent = extent
        indices = tuple(Index(0, extent) for _ in range(self.dimension))
        if indices:
            constraints = (pym.Comparison(reduce(operator.add, indices), "<", extent), )
        else:
            constraints = ()
        super().__init__(indices, constraints, codimension, variant_tag=variant_tag)

    @property
    @abc.abstractmethod
    def dimension(self):
        pass


class PointEntitySet(SimplexEntitySet):

    """A representation of a single (zero-dimensional) point."""

    dimension = 0

    def linear_index_map(self, index_exprs, index_order):
        assert len(index_exprs) == 0
        return 0


class IntervalEntitySet(SimplexEntitySet):

    """A representation of some number of intervals."""
    dimension = 1

    def linear_index_map(self, index_exprs, index_order):
        i, = index_exprs
        return i


class PeriodicIntervalEntitySet(SimplexEntitySet):

    dimension = 1

    def linear_index_map(self, index_exprs, index_order):
        i, = index_exprs
        return i % self.extent


class TriangleEntitySet(SimplexEntitySet):
    """A representation of entities with a "triangle" constraint."""
    dimension = 2

    def linear_index_map(self, index_exprs, index_order):
        """index_order is ignored (just swap the order in index_exprs)"""
        i, j = index_exprs
        return triangular_linear_index_map(i, j, self.extent)


class TetrahedronEntitySet(SimplexEntitySet):
    """A representation of entities with a "tetrahedron" constraint."""

    dimension = 3

    def linear_index_map(self, index_exprs, index_order):
        """index_order is ignored (just swap the order in index_exprs)"""
        i, j, k = index_exprs
        return tetrahedral_linear_index_map(i, j, k, self.extent)


class TensorProductEntitySet(EntitySet):
    """A representation of a set of indices with tensor product
    structure.

    :arg factors: index sets for the directions in the tensor product.
    :arg variant_tag: Arbitrary distinguishing tag."""
    def __init__(self, *factors, variant_tag=None):
        self.factors = tuple(factors)
        if any(isinstance(f, TensorProductEntitySet) for f in self.factors):
            raise ValueError("Can't deal with nested tensor products sorry")
        indices = tuple(itertools.chain(*(f.indices for f in self.factors)))
        assert len(set(i.name for i in indices)) == len(indices), "Must provide unique index names"
        constraints = tuple(itertools.chain(*(f.constraints for f in self.factors)))
        codimension = sum(f.codimension for f in self.factors)
        super().__init__(indices, constraints, codimension, variant_tag=variant_tag)

    def linear_index_map(self, index_exprs, index_order):
        assert len(index_exprs) == sum(len(f.indices) for f in self.factors)
        assert len(set(index_order)) == len(self.factors)

        # FIXME: Debug this more thoroughly
        indices = index_exprs
        index_exprs = []
        for factor in self.factors:
            nindex = len(factor.indices)
            index_exprs.append(indices[:nindex])
            indices = indices[nindex:]
        index_exprs = tuple(index_exprs[i] for i in index_order)
        factors = tuple(self.factors[i] for i in index_order)
        strides = numpy.cumprod((1, ) + tuple(f.size for f in reversed(factors)))
        strides = tuple(reversed(strides[:-1]))

        expr = 0
        for stride, index_expr, factor in zip(strides, index_exprs, self.factors):
            expr = expr + factor.linear_index_map(index_expr, None)*stride
        return expr

    def __str__(self):
        factors = ", ".join(str(f) for f in self.factors)
        return "TensorProductEntitySet({}: {})".format(factors, self.isl_set)


class MeshTopology(metaclass=abc.ABCMeta):
    """Object representing some structured mesh pattern."""

    @property
    @abc.abstractmethod
    def embedding_dimension(self):
        pass

    @abc.abstractmethod
    def entity_variants(self, codim=None):
        pass

    @abc.abstractmethod
    def subentity_map(self, entity_set, subentity_set, multiindex, subentity_multiindex):
        pass

    @abc.abstractmethod
    def dual_subentity_map(self, subentity_set, entity_set, multiindex):
        pass
