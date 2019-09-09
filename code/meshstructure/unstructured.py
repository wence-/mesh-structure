import abc
import numbers
import operator
from functools import reduce

import FIAT
import ufl

from .symbolics import Point
from .topology import UnstructuredEntitySet, UnstructuredMeshTopology
from .utils import lazyattr


def binomial(n, k):
    return int(reduce(operator.mul, ((n + 1 - i)/i for i in range(1, k+1)), 1))


class UnstructuredSimplex(UnstructuredMeshTopology):
    def __init__(self, dimension):
        super().__init__(dimension)
        self.cell = ufl.cell.simplex(dimension)

    @lazyattr
    def entities(self):
        dim = self.dimension
        return dict((codim, UnstructuredEntitySet(binomial(dim+1, codim), dim, codim))
                    for codim in range(dim+1))

    def cone(self, point):
        raise NotImplementedError()

    def support(self, point):
        raise NotImplementedError()


class UnstructuredHyperCube(UnstructuredMeshTopology):
    def __init__(self, dimension):
        super().__init__(dimension)
        self.cell = ufl.cell.hypercube(dimension)

    @lazyattr
    def entities(self):
        dim = self.dimension
        nent = lambda d: 2**(dim - d)*binomial(dim, d)
        return dict((codim, UnstructuredEntitySet(nent(dim - codim), dim, codim))
                    for codim in range(dim+1))

    def cone(self, point):
        raise NotImplementedError()

    def support(self, point):
        raise NotImplementedError()


class FIATCell(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def fiat_cell(self):
        """The FIAT reference cell."""

    @lazyattr
    def cones(self):
        topology = self.fiat_cell.get_topology()
        cones = {}
        dimension = self.dimension
        for dim, entities in topology.items():
            cones[dimension - dim] = {}
            for ent, vertices in entities.items():
                vertices = frozenset(vertices)
                cone = []
                for e, verts in topology.get(dim - 1, {}).items():
                    if vertices.issuperset(verts):
                        cone.append(e)
                cones[dimension - dim][ent] = tuple(sorted(cone))
        return cones

    @lazyattr
    def supports(self):
        topology = self.fiat_cell.get_topology()
        supports = {}
        dimension = self.dimension
        for dim, entities in topology.items():
            supports[dimension - dim] = {}
            for ent, vertices in entities.items():
                vertices = frozenset(vertices)
                support = []
                for e, verts in topology.get(dim + 1, {}).items():
                    if vertices.issubset(verts):
                        support.append(e)
                supports[dimension - dim][ent] = tuple(sorted(support))
        return supports

    def cone(self, point):
        indices, eset = point
        index, = indices
        assert isinstance(index, numbers.Integral)
        assert index < eset.size
        codim = eset.codimension
        target = self.entities[codim + 1]
        return tuple(Point((c,), target) for c in self.cones[codim][index])

    def support(self, point):
        indices, eset = point
        index, = indices
        assert isinstance(index, numbers.Integral)
        assert index < eset.size
        codim = eset.codimension
        target = self.entities[codim - 1]
        # FIXME: Also return local entity
        return tuple(Point((s,), target) for s in self.supports[codim][index])


class FIATSimplex(FIATCell, UnstructuredSimplex):

    @lazyattr
    def fiat_cell(self):
        return FIAT.ufc_simplex(self.dimension)


class FIATHyperCube(FIATCell, UnstructuredHyperCube):

    @lazyattr
    def fiat_cell(self):
        return {0: FIAT.reference_element.Point(),
                1: FIAT.reference_element.UFCInterval(),
                2: FIAT.reference_element.UFCQuadrilateral(),
                3: FIAT.reference_element.UFCHexahedron()}[self.dimension]
