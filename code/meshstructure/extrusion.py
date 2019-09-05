import itertools
from enum import Enum

from .topology import IntervalEntitySet, MeshTopology, TensorProductEntitySet


__all__ = ("MeshExtrusion", )


class Tag(Enum):
    CELL = 0
    VERTEX = 1


class MeshExtrusion(MeshTopology):
    def __init__(self, base, nlevel):
        super().__init__(base, base.dimension + 1)
        self.nlevel = nlevel
        cells = IntervalEntitySet(nlevel, codimension=0, variant_tag=Tag.CELL)
        vertices = IntervalEntitySet(nlevel+1, codimension=1, variant_tag=Tag.VERTEX)

        entities = {}
        for (codim, esets) in topology.entities.items():
            # TODO: Flatten codim tuple?
            entities[(codim, 0)] = tuple(TensorProductEntitySet(eset, cells)
                                         for eset in esets)
            entities[(codim, 1)] = tuple(TensorProductEntitySet(eset, vertices)
                                         for eset in esets)

    def entity_variants(self, codim=None):
        if codim is None:
            return tuple(self.entities.values())
        try:
            b, e = codim
            return (self.entities[codim], )
        except TypeError:
            return tuple(v for k, v in self.entities.items() if sum(k) == codim)

    def cone(self, multiindex):
        raise NotImplementedError

    def support(self, multiindex):
        raise NotImplementedError
