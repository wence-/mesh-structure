import itertools
from enum import Enum

import numpy
import ufl

from .symbolics import Point
from .topology import (IntervalEntitySet, StructuredMeshTopology,
                       TensorProductEntitySet)
from .unstructured import UnstructuredHyperCube
from .utils import lazyattr

__all__ = ("HyperCubeRefinement", )


class Tag(Enum):
    CELL = 0
    VERTEX = 1


class HyperCubeRefinement(StructuredMeshTopology):

    def __init__(self, base, *cells_per_dimension):
        """Representation of structured refinement of a hypercube.

        :arg base: The base topology to refine.
        :arg cells_per_dimension: Number of cells in each direction."""
        # FIXME: Base is ignored everywhere else.
        assert isinstance(base, UnstructuredHyperCube)
        super().__init__(base, ufl.cell.hypercube(base.dimension))
        assert len(cells_per_dimension) == self.dimension
        self.cells_per_dimension = cells_per_dimension

    @lazyattr
    def entities(self):
        entities = {}
        for codim in range(self.dimension+1):
            ents = []
            # Sets of entities of given codim are created by selecting
            # cell and vertex intervals such that n_vertex_intervals  == codim
            #
            # This is a multiset permutation of [0] * (dimension-codim) + [1] * codim
            for vtx in itertools.combinations(range(self.dimension), codim):
                cells = tuple(IntervalEntitySet(n, cell=ufl.interval, codimension=0, variant_tag=Tag.CELL)
                              for n in self.cells_per_dimension)
                vertices = tuple(IntervalEntitySet(n + 1, cell=ufl.interval, codimension=1, variant_tag=Tag.VERTEX)
                                 for n in self.cells_per_dimension)
                idx = list(vtx)
                factors = numpy.asarray(cells)
                factors[idx] = numpy.asarray(vertices)[idx]
                ents.append(TensorProductEntitySet(*factors))
            entities[codim] = tuple(ents)
        return entities

    def cone(self, point):
        """Given indices into an entity set, produce the index
        expressions for the cone of the entity, that is, the entities
        with codimension 1 greater.

        :arg point: The point in the set.
        :returns: A tuple of points.
        """
        indices, eset = point
        assert len(indices) == len(eset.factors)
        codim = eset.codimension + 1
        targets = self.entity_variants(codimension=codim)
        exprs = []
        for target in targets:
            expr = zip(*((index, index if parent.variant_tag == child.variant_tag else index + 1)
                         for index, parent, child in zip(indices, eset.factors, target.factors)))
            exprs.extend(Point(e, target) for e in expr)
        return tuple(exprs)

    def support(self, point):
        """Given indices into an entity set, produce index expression
        for the support of the entity, that is, the entities with
        codimension 1 less.

        :arg point: The point in the set.
        :returns: A tuple of points. FIXME: not correct.
        """
        indices, eset = point
        assert len(indices) == len(eset.factors)
        codim = eset.codimension - 1
        targets = self.entity_variants(codimension=codim)
        exprs = []
        for target in targets:
            expr = zip(*((index, index if parent.variant_tag == child.variant_tag else index - 1)
                         for index, parent, child in zip(indices, eset.factors, target.factors)))
            exprs.extend(Point(e, target) for e in expr)
        raise NotImplementedError("Need to also determine local subentity")
        return tuple(exprs)
