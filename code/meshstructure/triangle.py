import itertools
from enum import Enum

import numpy
import ufl

from .symbolics import Point
from .topology import (IntervalEntitySet, StructuredMeshTopology,
                       TensorProductEntitySet)
from .unstructured import UnstructuredSimplex
from .utils import lazyattr

__all__ = ("TriangleRefinement", )


class Tag(Enum):
    CELL_LL = 0 # lower left
    CELL_UR = 1 # upper right
    EDGE_X  = 2 # edge in x direction
    EDGE_Y  = 3 # edge in y direction
    EDGE_XY = 4 # diagonal edge in xy direction
    VERTEX  = 5 #  


class TriangleRefinement(StructuredMeshTopology):

    def __init__(self, base, cells_per_dimension):
        """Representation of structured refinement of a traingle.

        :arg base: The base topology to refine.
        :arg cells_per_dimension: Number of cells in each direction (scalar)."""
        # FIXME: Base is ignored everywhere else.
        assert isinstance(base, UnstructuredSimplex)
        assert base.dimension == 2
        super().__init__(base, ufl.cell.hypercube(base.dimension))
        assert type(cells_per_dimension) is int
        self.cells_per_dimension = cells_per_dimension

    @lazyattr
    def entities(self):
        entities = {}
        # cells
        entities[0] = tuple(
            TriangularEntitySet(self.cells_per_dimension,   ufl.triangle, codimension=0, variant_tag=Tag.CELL_LL),
            TriangularEntitySet(self.cells_per_dimension-1, ufl.triangle, codimension=0, variant_tag=Tag.CELL_UR))
        # faces
        entities[1] = tuple(
            TriangularEntitySet(self.cells_per_dimension,   ufl.triangle, codimension=1, variant_tag=Tag.EDGE_X),
            TriangularEntitySet(self.cells_per_dimension,   ufl.triangle, codimension=1, variant_tag=Tag.EDGE_Y),
            TriangularEntitySet(self.cells_per_dimension,   ufl.triangle, codimension=1, variant_tag=Tag.EDGE_XY))
        # vertices
        entities[2] = tuple(
            TriangularEntitySet(self.cells_per_dimension+1, ufl.triangle, codimension=2, variant_tag=Tag.Tag.VERTEX))
        return entities

    def cone(self, point):
        """Given indices into an entity set, produce the index
        expressions for the cone of the entity, that is, the entities
        with codimension 1 greater.

        :arg point: The point in the set.
        :returns: A tuple of points.
        """
        raise NotImplementedError
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
        raise NotImplementedError
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
