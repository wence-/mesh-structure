from enum import Enum

import ufl

from .symbolics import Point
from .topology import TensorProductEntitySet, UnstructuredEntitySet
from .utils import as_tuple, lazyattr

__all__ = ("DataLayout", "FIATSimplexDataLayout", "FIATHyperCubeDataLayout")


class Tag(Enum):
    DOF = 0
    DATA_LAYOUT = 1


# FIXME: Is this the right way to think about data, just as another direction for a set?
# It might get in the way later?
class DataSet(TensorProductEntitySet):
    def __init__(self, entity_set, dofs_per_point):
        # FIXME: Sort of doesn't align with the rest.
        factors = (entity_set, UnstructuredEntitySet(dofs_per_point,
                                                     cell=ufl.vertex,
                                                     codimension=0,
                                                     variant_tag=Tag.DOF), )
        super().__init__(*factors, variant_tag=Tag.DATA_LAYOUT)


class DataLayout(object):

    def __init__(self, topology, dofs_per_codim_entity):
        """Build a data layout object (like a function space).

        :arg topology: The mesh topology.
        :arg dofs_per_codim_entity: Number of degrees of freedom on
            each entity of a given codimension."""
        self.topology = topology
        self.datasets = dict((eset, DataSet(eset, dofs_per_codim_entity.get(eset.codimension, 0)))
                             for eset in topology.entity_variants())

    @lazyattr
    def size(self):
        return sum(e.size*d.size for e, d in self.datasets.items())

    def closure(self, point):
        """Return expressions for all the dofs supported on the
        topological closure of a given point.

        :arg point: The point (should be a point in the topology).
        :returns: A tuple of expressions indexing all the dofs.
        """
        tclosure = self.topology.closure(point)
        exprs = []
        for p in tclosure:
            indices, tset = p
            dset = self.datasets.get(tset, None)
            if dset is None:
                continue
            dofs = dset.factors[-1]
            if dofs.size > 0:
                exprs.append(Point(indices + dofs.indices, dset))
        return tuple(exprs)


class FIATDataLayout(DataLayout):
    def __init__(self, topology, element):
        from FIAT.reference_element import flatten_entities, TensorProductCell
        assert element.cell == topology.fiat_cell
        entity_dofs = element.entity_dofs()
        # UGH
        if isinstance(element.cell, TensorProductCell):
            entity_dofs = flatten_entities(entity_dofs)
        dofs_per_codim_entity = {}
        for e, v in entity_dofs.items():
            ndof, = set(map(len, v.values()))
            codimension = topology.dimension - sum(as_tuple(e))
            dofs_per_codim_entity[codimension] = ndof
        super().__init__(topology, dofs_per_codim_entity)
        self.element = element


class FIATSimplexDataLayout(FIATDataLayout):
    pass


class FIATHyperCubeDataLayout(FIATDataLayout):
    # FIAT has a really weird numbering for hypercube elements
    def closure(self, point):
        from FIAT.reference_element import flatten_entities, TensorProductCell
        tclosure = self.topology.closure(point)
        exprs = [None for _ in range(self.element.space_dimension())]
        entity_dofs = self.element.entity_dofs()
        if isinstance(self.element.cell, TensorProductCell):  # UGH!
            entity_dofs = flatten_entities(entity_dofs)
        for p in tclosure:
            indices, tset = p
            dset = self.datasets.get(tset, None)
            if dset is None:
                continue
            dofs = dset.factors[-1]
            if dofs.size > 0:
                i, = indices
                idxs = entity_dofs[tset.dimension][i]
                for k, j in enumerate(idxs):
                    exprs[j] = Point(indices + (k, ), dset)
        assert not any(e is None for e in exprs)
        return tuple(exprs)
