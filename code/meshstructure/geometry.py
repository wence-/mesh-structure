import abc
import operator
from functools import reduce

import ufl

from .hypercube import Tag, HypercubeRefinement

__all__ = ("AffineMappedHyperCube", )


class MeshGeometry(metaclass=abc.ABCMeta):

    def __init__(self, topology):
        self.topology = topology

    @abc.abstractmethod
    def geometric_quantity(self, entity_set, multiindex, geometric_quantity):
        """Return geometric information, to be interpreted by code generation.

        :arg entity_set: The set of entities to get the geometry for.
        :arg multiindex: An index describing a point in the entity set.
        :arg geometric_quantity: A UFL GeometricQuantity object indicating what
           geometric information is required. """
        pass


class AffineMappedHyperCube(MeshGeometry, ufl.Mesh):

    def __init__(self, topology):
        assert isinstance(topology, HypercubeRefinement)
        super().__init__(topology)
        # FIXME: This is really topological dimension
        gdim = topology.embedding_dimension
        cell = ufl.TensorProductCell(*(ufl.interval for _ in range(gdim)))
        element = ufl.TensorProductElement(*(ufl.FiniteElement("P", ufl.interval, 1)
                                             for _ in range(gdim)),
                                           cell=cell)
        ufl.Mesh.__init__(self, ufl.VectorElement(element, dim=gdim))

    def geometric_quantity(self, entity_set, multiindex, geometric_quantity):
        """Returns an expression, in terms of a "global" reference
        entity geometric quantity, for the same geometric quantity on
        the refined entity."""
        if isinstance(geometric_quantity, ufl.JacobianDeterminant):
            # FIXME: error checking that geometric_quantity matches codim of entity_set
            # Number of entities in the projection of the hypercube onto the current submanifold.
            scaling = reduce(operator.mul,
                             (f.size for f in entity_set.factors if f.variant_tag == Tag.CELL))
            return (1/scaling) * geometric_quantity
        else:
            raise NotImplementedError("Only detJ provided so far")
