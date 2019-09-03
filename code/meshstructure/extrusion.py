import ufl

from .topology import (IntervalEntitySet, PointEntitySet, StructureBase,
                       TensorProductEntitySet)


__all__ = ("Extrusion", )


class Extrusion(StructureBase):
    def __init__(self, base, nlevel):
        super().__init__()
        self.nlevel = nlevel
        if base in {ufl.triangle}:
            # OK
            base_cell = PointEntitySet(1)
            base_edge = IntervalEntitySet(base.num_edges())
            base_vertex = IntervalEntitySet(base.num_vertices())
        else:
            # But then what about if the base thing is structured?
            raise NotImplementedError()

        cells = TensorProductEntitySet(base_cell,
                                       IntervalEntitySet(nlevel),
                                       variant_tag=ufl.TensorProductCell(base,
                                                                         ufl.interval))
        hfaces = TensorProductEntitySet(base_cell,
                                        IntervalEntitySet(nlevel+1),
                                        variant_tag=ufl.TensorProductCell(base,
                                                                          ufl.vertex))
        vfaces = TensorProductEntitySet(base_edge,
                                        IntervalEntitySet(nlevel),
                                        variant_tag=ufl.TensorProductCell(ufl.interval,
                                                                          ufl.interval))
        vedges = TensorProductEntitySet(base_vertex,
                                        IntervalEntitySet(nlevel),
                                        variant_tag=ufl.TensorProductCell(ufl.vertex,
                                                                          ufl.interval))
        hedges = TensorProductEntitySet(base_edge,
                                        IntervalEntitySet(nlevel+1),
                                        variant_tag=ufl.TensorProductCell(ufl.interval,
                                                                          ufl.vertex))
        vertices = TensorProductEntitySet(base_vertex,
                                          IntervalEntitySet(nlevel+1),
                                          variant_tag=ufl.TensorProductCell(ufl.vertex,
                                                                            ufl.vertex))
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
