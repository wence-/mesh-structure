import islpy as isl
from pymbolic import primitives as pym
import abc
import ufl


class EntitySet(object):
    def __init__(self, polyhedral_set, variant_tag):
        self.polyhedral_set = polyhedral_set
        self.variant_tag = variant_tag


class StructureBase(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def entity_variants(self, codim=None):
        pass

    @abc.abstractmethod
    def subentity_map(self, entity_set, subentity_set, multiindex, subentity_multiindex):
        pass

    @abc.abstractmethod
    def dual_subentity_map(self, subentity_set, entity_set, multiindex):
        pass


class Extrusion(StructureBase):
    def __init__(self, nlevel):
        super().__init__()
        self.nlevel = nlevel
        vars = isl.make_zero_and_vars("i,j")
        # Extruded triangular prism
        cells = EntitySet(vars[0].le_set(vars["i"]) & vars["i"].lt_set(0 + nlevel),
                          ufl.TensorProductCell(ufl.triangle, ufl.interval))
        hfaces = EntitySet(vars[0].le_set(vars["i"]) & vars["i"].le_set(0 + nlevel),
                           ufl.triangle)
        vfaces = EntitySet(vars[0].le_set(vars["i"]) & vars["i"].lt_set(0 + nlevel)
                           & vars[0].le_set(vars["j"]) & vars["j"].lt_set(vars[0] + 3),
                           ufl.quadrilateral)
        hedges = EntitySet(vars[0].le_set(vars["i"]) & vars["i"].le_set(0 + nlevel)
                           & vars[0].le_set(vars["j"]) & vars["j"].lt_set(vars[0] + 3),
                           ufl.interval)
        vedges = EntitySet(vars[0].le_set(vars["i"]) & vars["i"].lt_set(0 + nlevel)
                           & vars[0].le_set(vars["j"]) & vars["j"].lt_set(vars[0] + 3),
                           ufl.interval)
        vertices = EntitySet(vars[0].le_set(vars["i"]) & vars["i"].le_set(0 + nlevel)
                             & vars[0].le_set(vars["j"]) & vars["j"].lt_set(vars[0] + 3),
                             ufl.vertex)

        self.entities = {(0, 0): cells,
                         (0, 1): hfaces,
                         (1, 0): vfaces,
                         (1, 1): hedges,
                         (2, 0): vedges,
                         (2, 1): vertices}

        def cellhface(mi, s):
            assert s in {0, 1}
            mi, = mi
            return (mi + s, )

        def cellvface(mi, s):
            assert s in {0, 1, 2}
            mi, = mi
            return (mi, s)

        def hfacehedge(mi, s):
            assert s in {0, 1, 2}
            mi, = mi
            return (mi, s)

        def vfacehedge(mi, s):
            assert s in {0, 1}
            i, j = mi
            return (i + s, j)

        def vfacevedge(mi, s):
            assert s in {0, 1}
            i, j = mi
            return (i, (j + s) % 3)

        def hedgevertex(mi, s):
            assert s in {0, 1}
            i, j = mi
            return (i, (j + s) % 3)

        def vedgevertex(mi, s):
            assert s in {0, 1}
            i, j = mi
            return (i + s, j)

        self.maps_ = {(cells, hfaces): cellhface,
                      (cells, vfaces): cellvface,
                      (hfaces, hedges): hfacehedge,
                      (vfaces, hedges): vfacehedge,
                      (vfaces, vedges): vfacevedge,
                      (hedges, vertices): hedgevertex,
                      (vedges, vertices): vedgevertex}

        def hfacecell(mi):
            i, = mi
            return ((i - 1, ), (i, ))

        def vfacecell(mi):
            i, _ = mi
            return ((i, ), )

        def hedgehface(mi):
            i, j = mi
            return ((i, ), )

        def hedgevface(mi):
            i, j = mi
            return ((i - 1, j), (i, j))

        def vedgevface(mi):
            i, j = mi
            return ((i, (j + 2) % 3), (i, j))

        def vertexhedge(mi):
            i, j = mi
            return ((i, (j + 2) % 3), (i, j))

        def vertexvedge(mi):
            i, j = mi
            return ((i - 1, j), (i, j))

        self.dual_maps_ = {(hfaces, cells): hfacecell,
                           (vfaces, cells): vfacecell,
                           (hedges, hfaces): hedgehface,
                           (hedges, vfaces): hedgevface,
                           (vedges, vfaces): vedgevface,
                           (vertices, hedges): vertexhedge,
                           (vertices, vedges): vertexvedge}

    def entity_variants(self, codim=None):
        if codim is None:
            return tuple(self.entities.values())
        try:
            b, e = codim
            return (self.entities[codim], )
        except TypeError:
            return tuple(v for k, v in self.entities.items() if sum(k) == codim)

    def subentity_map(self, eset, sset, multiindex, subentity):
        try:
            return self.maps_[(eset, sset)](multiindex, subentity)
        except KeyError:
            raise ValueError("No map between {} and {}".format(eset, sset))

    def dual_subentity_map(self, sset, eset, multiindex):
        try:
            return self.dual_maps_[(sset, eset)](multiindex)
        except KeyError:
            raise ValueError("No dual map between {} and {}".format(sset, eset))
