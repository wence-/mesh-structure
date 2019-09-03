import islpy as isl
import itertools
from pymbolic import primitives as pym
from pymbolic.mapper.evaluator import evaluate
import abc
import numpy
import numbers
from functools import reduce, singledispatch, partial
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
    """A representation of some set of entities.

    :arg indices: :class:`Index` objects encoding points in the set.
    :arg constraints: constraints on the indices (expressions), must
        be expressible in Presburger arithmetic.
    :arg variant_tag: Arbitrary data used to distinguish this set."""
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
        """A map from indices in the set into a linear index.

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


class SimplexEntitySet(EntitySet):
    """An entity set describing some number of simplices.

    :arg extent: Max extent in each direction.

    Produces an index set with

    0 <= indices < extent, sum(indices) < extent

    Where len(indices) == dimension
    """
    def __init__(self, extent, variant_tag=None):
        assert isinstance(extent, numbers.Integral)
        self.extent = extent
        indices = tuple(Index(0, extent) for _ in range(self.dimension))
        if indices:
            constraints = (pym.Comparison(reduce(operator.add, indices), "<", extent), )
        else:
            constraints = ()
        super().__init__(indices, constraints, variant_tag=variant_tag)

    @abc.abstractproperty
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


# This would allow us to represent entities in a block-structured grid
# in a staggered fashion, I think. Consider a 2D quad, we need
#
# cells = TensorProductEntitySet(IntervalEntitySet(i, N), IntervalEntitySet(j, N))
# vertices = TensorProductEntitySet(IntervalEntitySet(k, N+1), IntervalEntitySet(l, N+1))
# | faces = TensorProductEntitySet(IntervalEntitySet(k, N+1), IntervalEntitySet(j, N))
# - faces = TensorProductEntitySet(IntervalEntitySet(i, N), IntervalEntitySet(l+1, N+1))
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
        super().__init__(indices, constraints, variant_tag=variant_tag)

    def linear_index_map(self, index_exprs, index_order):
        assert len(index_exprs) == sum(len(f.indices) for f in self.factors)
        assert len(index_order) == len(self.factors)
        factors = tuple(self.factors[i] for i in index_order)
        strides = tuple(numpy.cumprod((1, ) + tuple(f.size for f in reversed(factors)))[:-1][::-1][index_order])
        expr = 0
        for factor, stride in zip(factors, strides):
            nindex = len(factor.indices)
            index_expr = index_exprs[:nindex]
            index_exprs = index_exprs[nindex:]
            expr = expr + factor.linear_index_map(index_expr, None)*stride
        return expr

    def __str__(self):
        factors = ", ".join(str(f) for f in self.factors)
        return "TensorProductEntitySet({}: {})".format(factors, self.isl_set)


class StructureBase(metaclass=abc.ABCMeta):
    """Object representing some structured mesh pattern."""
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
    def __init__(self, base, nlevel):
        super().__init__()
        self.nlevel = nlevel
        if base in {ufl.triangle}:
            # OK
            base_cell = PointEntitySet(1)
            base_edge = IntervalEntitySet(base.num_edges())
            base_vertex = IntervalEntitySet(base.num_vertices())
        else:
            # But then what about if the base thing is structured?
            raise NotImplementedError()

        cells = TensorProductEntitySet(base_cell,
                                       IntervalEntitySet(nlevel),
                                       variant_tag=ufl.TensorProductCell(base,
                                                                         ufl.interval))
        hfaces = TensorProductEntitySet(base_cell,
                                        IntervalEntitySet(nlevel+1),
                                        variant_tag=ufl.TensorProductCell(base,
                                                                          ufl.vertex))
        vfaces = TensorProductEntitySet(base_edge,
                                        IntervalEntitySet(nlevel),
                                        variant_tag=ufl.TensorProductCell(ufl.interval,
                                                                          ufl.interval))
        vedges = TensorProductEntitySet(base_vertex,
                                        IntervalEntitySet(nlevel),
                                        variant_tag=ufl.TensorProductCell(ufl.vertex,
                                                                          ufl.interval))
        hedges = TensorProductEntitySet(base_edge,
                                        IntervalEntitySet(nlevel+1),
                                        variant_tag=ufl.TensorProductCell(ufl.interval,
                                                                          ufl.vertex))
        vertices = TensorProductEntitySet(base_vertex,
                                          IntervalEntitySet(nlevel+1),
                                          variant_tag=ufl.TensorProductCell(ufl.vertex,
                                                                            ufl.vertex))
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


# TODO: use an enum
CELL = object()
VERTEX = object()


class HypercubeRefinement(StructureBase):

    """Representation of structured refinement of a hypercube.

    :arg cells_per_dimension: Number of cells in each direction."""
    def __init__(self, *cells_per_dimension):
        self.entities = {}
        cells = tuple(IntervalEntitySet(n, variant_tag=CELL) for n in cells_per_dimension)
        vertices = tuple(IntervalEntitySet(n+1, variant_tag=VERTEX) for n in cells_per_dimension)
        dimension = len(cells_per_dimension)
        for codim in range(dimension+1):
            ents = []
            # Sets of entities of given codim are created by selecting
            # cell and vertex intervals such that n_vertex_intervals  == codim
            #
            # This is a multiset permutation of [0] * (dimension-codim) + [1] * codim
            for vtx in itertools.combinations(range(dimension), codim):
                idx = list(vtx)
                factors = numpy.asarray(cells)
                factors[idx] = numpy.asarray(vertices)[idx]
                ents.append(TensorProductEntitySet(*factors))
            self.entities[codim] = tuple(ents)

    def cone(self, indices, eset):
        """Given indices into an entity set, produce the index
        expressions for the cone of the entity, that is, the entities
        with codimension 1 greater.

        :arg indices: The indices describing a point in the set.
        :arg eset: The entity set
        :returns: A tuple of two-tuples, each of the form
            (([index_expr1, ...]), codim+1-set)
        """
        assert len(indices) == len(eset.factors)
        codim = sum(f.variant_tag == VERTEX for f in eset.factors) + 1
        targets = self.entity_variants(codim=codim)
        exprs = []
        for target in targets:
            expr = zip(*((index, index if parent.variant_tag == child.variant_tag else index + 1)
                         for index, parent, child in zip(indices, eset.factors, target.factors)))
            exprs.append((tuple(expr), target))
        return tuple(exprs)

    def support(self, indices, eset):
        """Given indices into an entity set, produce index expression
        for the support of the entity, that is, the entities with
        codimension 1 less.

        :arg indices: The indices describing a point in the set.
        :arg eset: The entity set
        :returns: A tuple of two-tuples, each of the form
            (([index_expr1, ...]), codim-1-set)
        """
        assert len(indices) == len(eset.factors)
        codim = sum(f.variant_tag == VERTEX for f in eset.factors) - 1
        targets = self.entity_variants(codim=codim)
        exprs = []
        for target in targets:
            expr = zip(*((index, index if parent.variant_tag == child.variant_tag else index - 1)
                         for index, parent, child in zip(indices, eset.factors, target.factors)))
            exprs.append((tuple(expr), target))
        return tuple(exprs)

    def subentity_map(self, eset, sset, indices, subentity):
        """Map some indices on a given entity set into indices on an
        immediately neighbouring subentity set, picking out a given
        subentity.

        :arg eset: The source entity set
        :arg sset: The target entity set, with codimension 1 greater.
        :arg indices: Indices into eset
        :arg subentity: Which subentity in sset from the map to pick
            out.
        :returns: Indices into sset.
        """
        candidates = self.cone(indices, eset)
        cone, = (cone for cone, candidate in candidates if candidate == sset)
        return cone[subentity]

    def dual_subentity_map(self, sset, eset, indices):
        """Map some indices on a given entity set into indices on an
        immediately neighbouring super-entity set.

        :arg sset: The source entity set
        :arg eset: The target entity set, with codimension 1 less.
        :arg indices: Indices into eset
        :returns: A tuple of indices into eset which index all the
            neighbouring points of the point provided in sset.
        """
        candidates = self.support(indices, sset)
        support, = (support for support, candidate in candidates if candidate == eset)
        return support

    def entity_variants(self, codim=None):
        if codim is None:
            return tuple(self.entities.values())
        else:
            return self.entities.get(codim, ())
