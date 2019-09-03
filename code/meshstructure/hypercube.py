import itertools
from enum import Enum

import numpy

from .topology import IntervalEntitySet, StructureBase, TensorProductEntitySet
from .utils import lazyprop


__all__ = ("HypercubeRefinement", )


class Tag(Enum):
    CELL = 0
    VERTEX = 1


class HypercubeRefinement(StructureBase):

    """Representation of structured refinement of a hypercube.

    :arg cells_per_dimension: Number of cells in each direction."""
    def __init__(self, *cells_per_dimension):
        self.entities = {}
        cells = tuple(IntervalEntitySet(n, codimension=0, variant_tag=Tag.CELL)
                      for n in cells_per_dimension)
        vertices = tuple(IntervalEntitySet(n+1, codimension=1, variant_tag=Tag.VERTEX)
                         for n in cells_per_dimension)
        dimension = len(cells_per_dimension)
        self.embedding_dimension = dimension
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

    @lazyprop
    def embedding_dimension(self):
        pass

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
        codim = sum(f.variant_tag == Tag.VERTEX for f in eset.factors) + 1
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
        codim = sum(f.variant_tag == Tag.VERTEX for f in eset.factors) - 1
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
