import itertools
from enum import Enum

import numpy

from .topology import IntervalEntitySet, MeshTopology, TensorProductEntitySet
from .utils import lazyattr


__all__ = ("HypercubeRefinement", )


class Tag(Enum):
    CELL = 0
    VERTEX = 1


class HypercubeRefinement(MeshTopology):

    """Representation of structured refinement of a hypercube.

    :arg cells_per_dimension: Number of cells in each direction."""
    def __init__(self, base, *cells_per_dimension):
        super().__init__(base, base.dimension)
        assert len(cells_per_dimension) == self.dimension
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

    @lazyattr
    def embedding_dimension(self):
        pass

    def cone(self, multiindex):
        """Given indices into an entity set, produce the index
        expressions for the cone of the entity, that is, the entities
        with codimension 1 greater.

        :arg multiindex: The index describing a point in the set.
        :returns: A tuple of multiindices.
        """
        indices, eset = multiindex
        assert len(indices) == len(eset.factors)
        codim = eset.codimension + 1
        targets = self.entity_variants(codim=codim)
        exprs = []
        for target in targets:
            expr = zip(*((index, index if parent.variant_tag == child.variant_tag else index + 1)
                         for index, parent, child in zip(indices, eset.factors, target.factors)))
            exprs.extend(MultiIndex(e, target) for e in expr)
        return tuple(exprs)

    def support(self, multiindex):
        """Given indices into an entity set, produce index expression
        for the support of the entity, that is, the entities with
        codimension 1 less.

        :arg multiindex: The index describing a point in the set.
        :returns: A tuple of multiindices. FIXME: not correct.
        """
        assert len(indices) == len(eset.factors)
        codim = eset.codimension - 1
        targets = self.entity_variants(codim=codim)
        exprs = []
        for target in targets:
            expr = zip(*((index, index if parent.variant_tag == child.variant_tag else index - 1)
                         for index, parent, child in zip(indices, eset.factors, target.factors)))
            exprs.extend(MultiIndex(e, target) for e in expr)
        raise NotImplementedError("Need to also determine local subentity")
        return tuple(exprs)

    def entity_variants(self, codim=None):
        if codim is None:
            return tuple(self.entities.values())
        else:
            return self.entities.get(codim, ())
