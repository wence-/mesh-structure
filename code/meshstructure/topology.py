import abc
import itertools
import numbers
import operator
from functools import reduce, singledispatch

import islpy as isl
import numpy
import ufl
from pymbolic import primitives as pym
from pymbolic.mapper.evaluator import evaluate

from .symbolics import Index
from .utils import lazyattr

__all__ = ("UnstructuredEntitySet", "IntervalEntitySet",
           "TriangleEntitySet", "TetrahedronEntitySet",
           "TensorProductEntitySet")


class EntitySet(metaclass=abc.ABCMeta):
    """A representation of some set of entities.

    :arg indices: :class:`Index` objects encoding points in the set.
    :arg constraints: constraints on the indices (expressions), must
        be expressible in Presburger arithmetic.
    :arg cell: The type of cells this set contains, a UFL cell.
    :arg codimension: The codimension of the entities this entity set represents.
    :arg variant_tag: Arbitrary data used to distinguish this set."""
    def __init__(self, indices, constraints, *, cell=None, codimension=None, variant_tag=None):
        self.indices = tuple(indices)
        self.constraints = tuple(constraints)
        self.variant_tag = variant_tag
        self.cell = cell
        self.codimension = codimension

    @lazyattr
    def dimension(self):
        return self.cell.topological_dimension()

    @lazyattr
    def index_extents(self):
        """The extent of each index in the set."""
        return tuple(i.extent for i in self.indices)

    @abc.abstractmethod
    def linear_index_map(self, index_exprs):
        """A map from indices in the set into a linear index.

        :arg index_exprs: Expressions for each index the map should apply to.
        """

    def boundaries(self):
        """Return subsets of this entity set on the boundary of the polyhedral domain."""
        raise NotImplementedError

    def interior(self):
        """Return a subset of this entity set on the interior of the polyhedral domain."""
        raise NotImplementedError

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
            return reduce(operator.add, (translate(c, v) for c in expr.children))

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
        return "{}({}, {})".format(type(self).__name__, self.isl_set, self.cell)

    __repr__ = __str__


class UnstructuredEntitySet(EntitySet):
    def __init__(self, size, *, cell, codimension, variant_tag=None):
        indices = (Index(0, size), )
        super().__init__(indices, (), cell=cell, codimension=codimension, variant_tag=variant_tag)

    def linear_index_map(self, index_exprs):
        index, = index_exprs
        return index


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
    """An entity set describing some number of entities with a simplex
    constraint on the indexing.

    :arg extent: Max extent in each direction.

    Produces an index set with

    0 <= indices < extent, sum(indices) < extent

    Where len(indices) == dimension
    """
    def __init__(self, extent, *, cell, codimension, variant_tag=None):
        assert isinstance(extent, numbers.Integral)
        self.extent = extent
        dimension = cell.topological_dimension()
        indices = tuple(Index(0, extent) for _ in range(dimension))
        if indices:
            constraints = (pym.Comparison(reduce(operator.add, indices), "<", extent), )
        else:
            constraints = ()
        super().__init__(indices, constraints, cell=cell, codimension=codimension, variant_tag=variant_tag)


class IntervalEntitySet(EntitySet):
    def __init__(self, extent_or_index, *, cell, codimension, variant_tag=None):
        assert isinstance(extent_or_index, (numbers.Integral, Index))
        if isinstance(extent_or_index, numbers.Integral):
            self.extent = extent_or_index
            indices = Index(0, extent_or_index)
        else:
            self.extent = extent_or_index.extent
            indices = extent_or_index
        if indices:
            constraints = (pym.Comparison(indices, "<", indices.hi),
                           pym.Comparison(indices, ">=", indices.lo))
        else:
            constraints = ()
        super().__init__((indices,), constraints, cell=cell, codimension=codimension, variant_tag=variant_tag)

    """A representation of some number of intervals."""

    def linear_index_map(self, index_exprs):
        i, = index_exprs
        return i

    def boundaries(self):
        lo = self.indices[0].lo
        hi = self.indices[0].hi

        return (IntervalEntitySet(Index(lo, lo + 1), cell=self.cell, codimension=self.codimension,
                                  variant_tag=self.variant_tag),
                IntervalEntitySet(Index(hi - 1, hi), cell=self.cell, codimension=self.codimension,
                                  variant_tag=self.variant_tag),)


class TriangleEntitySet(SimplexEntitySet):
    """A representation of entities with a "triangle" constraint."""

    def linear_index_map(self, index_exprs):
        i, j = index_exprs
        return triangular_linear_index_map(i, j, self.extent)


class TetrahedronEntitySet(SimplexEntitySet):
    """A representation of entities with a "tetrahedron" constraint."""

    def linear_index_map(self, index_exprs):
        i, j, k = index_exprs
        return tetrahedral_linear_index_map(i, j, k, self.extent)


class TensorProductEntitySet(EntitySet):
    """A representation of a set of indices with tensor product
    structure.

    :arg factors: index sets for the directions in the tensor product.
    :arg variant_tag: Arbitrary distinguishing tag."""
    def __init__(self, *factors, variant_tag=None):
        self.factors = tuple()

        def flatten_factor(factor):
            if isinstance(factor, TensorProductEntitySet):
                for f in factor.factors:
                    yield from flatten_factor(f)
            else:
                yield factor

        for factor in factors:
            for f_ in flatten_factor(factor):
                self.factors += (f_,)

        if any(isinstance(f, TensorProductEntitySet) for f in self.factors):
            raise ValueError("Can't deal with nested tensor products sorry")
        indices = tuple(itertools.chain(*(f.indices for f in self.factors)))
        assert len(set(i.name for i in indices)) == len(indices), "Must provide unique index names"
        constraints = tuple(itertools.chain(*(f.constraints for f in self.factors)))
        codimension = sum(f.codimension for f in self.factors)
        cell = ufl.TensorProductCell(*(f.cell for f in self.factors))
        super().__init__(indices, constraints, cell=cell, codimension=codimension, variant_tag=variant_tag)

    def linear_index_map(self, index_exprs):
        assert len(index_exprs) == sum(len(f.indices) for f in self.factors)

        strides = numpy.cumprod((1, ) + tuple(f.size for f in reversed(self.factors)))
        strides = tuple(reversed(strides[:-1]))

        expr = 0
        for stride, factor in zip(strides, self.factors):
            nindex = len(factor.indices)
            index_expr = index_exprs[:nindex]
            index_exprs = index_exprs[nindex:]
            expr = expr + factor.linear_index_map(index_expr)*stride
        return expr

    def __str__(self):
        factors = ", ".join(str(f) for f in self.factors)
        return "TensorProductEntitySet({}: {})".format(factors, self.isl_set)

    def boundaries(self):
        factors_set = set(self.factors)

        boundaries = tuple()
        for factor in self.factors:
            factor_boundaries = factor.boundaries()
            for boundary in factor_boundaries:
                boundaries = boundaries + (TensorProductEntitySet(boundary, *(tuple(factors_set - {factor})),
                                                                  variant_tag=self.variant_tag),)
        return boundaries


class MeshTopology(metaclass=abc.ABCMeta):
    def __init__(self, cell):
        """Object representing some mesh topology.

        :arg cell: The cell type for the mesh, a UFL cell."""
        assert cell.topological_dimension() == cell.geometric_dimension()
        self.cell = cell

    @lazyattr
    def dimension(self):
        return self.cell.topological_dimension()

    @property
    @abc.abstractmethod
    def entities(self):
        """A dict mapping codimension to tuples of entity sets of given codimension"""

    @lazyattr
    def valid_entities(self):
        return frozenset(self.entity_variants())

    def entity_variants(self, *, codimension=None):
        """Return a tuple of entity sets of given codimension, separated by variant.

        :arg codimension: The codimension to select, or None for all entity sets."""
        if codimension is None:
            return tuple(itertools.chain.from_iterable(self.entities.values()))
        else:
            return self.entities.get(codimension, ())

    @abc.abstractmethod
    def cone(self, point):
        """Return the codimension + 1 neighbours of a point.

        :arg point: a :class:`Point` object.
        :returns: A (possibly empty) tuple of points
        """

    @abc.abstractmethod
    def support(self, point):
        """Return the codimension - 1 neighbours of a point.

        :arg point: a :class:`Point` object.
        :returns: A (possibly empty) tuple of (support_point,
            local_subentity_index) pairs. Where local_subentity_index
            is the local index in the support entity of the point.
        """

    def closure(self, point):
        """All points in the closure of the given point.

        :arg point: The point.
        :returns: A tuple of points in the topological closure of the
            given point."""
        seen = set([])
        fifo = [point]
        closure = []
        while fifo:
            p = fifo.pop()
            if p in seen:
                continue
            closure.append(p)
            seen.add(p)
            cone = self.cone(p)
            for c in cone:
                fifo.insert(0, c)
        return tuple(closure)

    def star(self, point):
        """All points in the star (dual of closure) of the given
        point.

        :arg point: The point:
        :returns: A tuple of points in the topological star of the
            given point."""
        seen = set([])
        fifo = [point]
        star = []
        while fifo:
            p = fifo.pop()
            if p in seen:
                continue
            star.append(p)
            seen.add(p)
            try:
                support, _ = zip(*self.support(p))
            except ValueError:
                support = ()
            for s in support:
                fifo.insert(0, s)
        return tuple(star)

    def index_relation(self, point, target):
        """Compute points for all entities in the target which
        are reachable from the source point.

        :arg multiindex: A :class:`Point` object.
        :arg target: An entity set.
        :returns: A (possibly empty) tuple of multiindices describing
            points in the target set.
        """
        _, source = point
        if source.codimension == target.codimension:
            return (point, )
        elif source.codimension > target.codimension:
            try:
                points, _ = zip(*self.support(point))
            except ValueError:
                points = ()
        else:
            points = self.cone(point)
        seen = set()
        points = tuple(itertools.chain(*(self.index_relation(p, target)
                                         for p in points)))
        filtered_points = []
        for p in points:
            if p not in seen:
                filtered_points.append(p)
                seen.add(p)
        return tuple(filtered_points)


class UnstructuredMeshTopology(MeshTopology):
    """An unstructured topology"""
    pass


class StructuredMeshTopology(MeshTopology):
    """A structured topology

    :arg base: The base topology being "refined".
    :arg cell: The cell type of the topology."""
    def __init__(self, base, cell):
        super().__init__(cell)
        self.base = base
