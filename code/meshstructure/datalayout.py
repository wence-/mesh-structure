import itertools
from enum import Enum

import numpy

from .topology import IntervalEntitySet, TensorProductEntitySet
from .utils import lazyattr


__all__ = ("DataLayout", )


class Tag(Enum):
    DOF = 0
    DATA_LAYOUT = 1


class DataSet(TensorProductEntitySet):
    def __init__(self, entity_set, dofs_per_entity):
        self.entity_set = entity_set
        factors = entity_set.factors + (IntervalEntitySet(dofs_per_entity, codimension=0,
                                                          variant_tag=Tag.DOF), )
        super().__init__(*factors, variant_tag=Tag.DATA_LAYOUT)


class DataLayout(object):

    def __init__(self, topology, dofs_per_entity):
        self.topology = topology
        self.datasets = dict((codim, tuple(DataSet(eset, dofs_per_entity[codim])
                                           for eset in esets))
                             for codim, esets in topology.entities.items())

    @lazyattr
    def size(self):
        return sum(d.size for d in itertools.chain(*self.datasets.values()))

    @lazyattr
    def ranges(self):
        sizes = []
        for codim in range(self.topology.dimension + 1):
            sizes.append(sum(d.size for d in itertools.chain(*self.datasets.values())
                             if d.codimension == codim))
        return (0, ) + tuple(numpy.cumsum(sizes))

    def range(self, codim):
        """Range of dofs [start, end) in local vector for data
        attached to entities of given codimension."""
        # FIXME: Doesn't work for offsetting in variants of same codim
        # FIXME: What if we want to interleave entities?
        return self.ranges[codim:codim+1]
