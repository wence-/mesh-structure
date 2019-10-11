from collections import defaultdict
from enum import Enum

import ufl

from .topology import (IntervalEntitySet, StructuredMeshTopology,
                       TensorProductEntitySet)
from .utils import lazyattr

__all__ = ("MeshExtrusion", )


class Tag(Enum):
    CELL = 0
    VERTEX = 1
    HORIZONTAL = 2
    VERTICAL = 3


class MeshExtrusion(StructuredMeshTopology):
    def __init__(self, base, nlevel):
        """An extruded topology.

        :arg base: The base topology to be extruded.
        :arg nlevel: The number of levels (cells) in the extruded topology."""
        super().__init__(base, ufl.TensorProductCell(base.cell, ufl.interval))
        self.nlevel = nlevel

    @lazyattr
    def entities(self):
        cells = IntervalEntitySet(self.nlevel, cell=ufl.interval, codimension=0, variant_tag=Tag.CELL)
        vertices = IntervalEntitySet(self.nlevel+1, cell=ufl.interval, codimension=1, variant_tag=Tag.VERTEX)
        entities = defaultdict(list)
        for (codim, esets) in self.base.entities.items():
            entities[codim].extend(TensorProductEntitySet(eset, cells,
                                                          variant_tag=Tag.VERTICAL)
                                   for eset in esets)
            entities[codim+1].extend(TensorProductEntitySet(eset, vertices,
                                                            variant_tag=Tag.HORIZONTAL)
                                     for eset in esets)
        return dict((k, tuple(v)) for k, v in entities.items())

    def cone(self, multiindex):
        raise NotImplementedError

    def support(self, multiindex):
        raise NotImplementedError
