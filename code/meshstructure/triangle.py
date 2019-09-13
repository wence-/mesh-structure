import itertools
from enum import Enum

import numpy
import ufl

from .symbolics import Point
from .topology import (TriangleEntitySet, StructuredMeshTopology)
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
        super().__init__(base, ufl.cell.triangle())
        assert isinstance(cells_per_dimension, numbers.Integral)
        self.cells_per_dimension = cells_per_dimension
        ########################
        ### generate entity sets
        self.entities = {}
        # cells
        self.entities[0] = tuple(
            TriangleEntitySet(self.cells_per_dimension,   ufl.triangle, codimension=0, variant_tag=Tag.CELL_LL),
            TriangleEntitySet(self.cells_per_dimension-1, ufl.triangle, codimension=0, variant_tag=Tag.CELL_UR))
        # faces
        self.entities[1] = tuple(
            TriangleEntitySet(self.cells_per_dimension,   ufl.interval, codimension=1, variant_tag=Tag.EDGE_X),
            TriangleEntitySet(self.cells_per_dimension,   ufl.interval, codimension=1, variant_tag=Tag.EDGE_Y),
            TriangleEntitySet(self.cells_per_dimension,   ufl.interval, codimension=1, variant_tag=Tag.EDGE_XY))
        # vertices
        self.entities[2] = tuple(
            TriangleEntitySet(self.cells_per_dimension+1, ufl.vertex,   codimension=2, variant_tag=Tag.VERTEX))
        #######################
        ### generate cone rules
        self.cones = {}
        # cell -> edge
        self.cones[self.entities[0][0]] = [([0,0],self.entities[1][2]), ([0,0],self.entities[1][0]), ([0,0],self.entities[1][1])] # ll
        self.cones[self.entities[0][1]] = [([0,1],self.entities[1][0]), ([1,0],self.entities[1][1]), ([0,0],self.entities[1][2])] # ur
        # edge -> vertex
        self.cones[self.entities[1][0]] = [([0,0],self.entities[2][0]), ([1,0],self.entities[2][0])] # x
        self.cones[self.entities[1][1]] = [([0,0],self.entities[2][0]), ([0,1],self.entities[2][0])] # y
        self.cones[self.entities[1][2]] = [([1,0],self.entities[2][0]), ([0,1],self.entities[2][0])] # xy
        ##########################
        ### generate support rules
        # edge -> cell
        self.supps[self.entities[1][0]] = [([-1,0],self.entities[0][0]), ([0,0],self.entities[0][1])]
        self.supps[self.entities[1][1]] = [([0,-1],self.entities[0][0]), ([0,0],self.entities[0][1])]
        self.supps[self.entities[1][2]] = [([0,0],self.entities[0][0]), ([0,0],self.entities[0][1])]
        # vertex -> edge
        self.supps[self.entities[2][0]] = [([-1,0],self.entities[1][0]), ([0,0],self.entities[1][0]),  # y
                                           ([0,-1],self.entities[1][1]), ([0,0],self.entities[1][1]),  # x
                                           ([-1,0],self.entities[1][2]), ([0,-1],self.entities[1][2])] # xy
        
    def cone(self, point):
        """Given indices into an entity set, produce the index
        expressions for the cone of the entity, that is, the entities
        with codimension 1 greater.

        :arg point: The point in the set.
        :returns: A tuple of points.
        """
        indices, eset = point
        return tuple(indices + shift for shift in cones[eset])

    def support(self, point):
        """Given indices into an entity set, produce index expression
        for the support of the entity, that is, the entities with
        codimension 1 less.

        :arg point: The point in the set.
        :returns: A tuple of points. FIXME: not correct.
        """
        raise NotImplementedError
        indices, eset = point
        return tuple(indices + shift for shift in supps[eset])
