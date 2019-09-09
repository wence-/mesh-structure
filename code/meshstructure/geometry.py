import abc

import ufl

from .utils import lazyattr

__all__ = ("MeshGeometry", )


class MeshGeometry(ufl.Mesh, metaclass=abc.ABCMeta):

    def __init__(self, topology, element):
        """A representation of """
        assert topology.dimension == element.topological_dimension()
        assert isinstance(element, ufl.VectorElement)
        self.topology = topology
        self.element = element
        ufl.Mesh.__init__(self, element)

    @lazyattr
    def geometric_dimension(self):
        """The embedding dimension."""
        return self.element.geometric_dimension()

    @lazyattr
    def topological_dimension(self):
        return self.topology.dimension

    @abc.abstractmethod
    def geometry_dofs(self, point):
        """An expression for the coordinate dofs in the closure
    of the provided point.

        :arg point: The entity to provide dofs for.
        :returns: An expression of shape (num_basis_functions, geometric_dimension)
        """
        pass
