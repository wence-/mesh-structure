import pytest

from meshstructure import (HyperCubeRefinement, MeshExtrusion)
from meshstructure.unstructured import (UnstructuredHyperCube, UnstructuredSimplex)
from meshstructure.extrusion import Tag


@pytest.mark.parametrize('dim', range(1, 4))
def test_hypercube_size(dim):
    h = HyperCubeRefinement(UnstructuredHyperCube(dim), *([5]*dim))
    cell_set, = h.entity_variants(codimension=0)
    assert cell_set.size == 5**dim

    vertex_set, = h.entity_variants(codimension=dim)
    assert vertex_set.size == 6**dim

    for codim in range(1, dim):
        sets = h.entity_variants(codimension=codim)
        set_sizes = [set.size for set in sets]
        assert all(set_sizes[0] == s for s in set_sizes)
        assert set_sizes[0] == 6**codim * 5**(dim-codim)


@pytest.mark.parametrize('base', [UnstructuredSimplex(2), UnstructuredHyperCube(2)])
def test_extrusion_size(base):
    extrusion = MeshExtrusion(base, 10)

    cell_set, = extrusion.entity_variants(codimension=0)
    assert cell_set.size == 10

    base_vertex_set, = base.entity_variants(codimension=2)
    vertex_set, = extrusion.entity_variants(codimension=3)
    assert vertex_set.size == base_vertex_set.size * 11

    face_sets = extrusion.entity_variants(codimension=1)
    h_face_set, = [s for s in face_sets if s.variant_tag == Tag.HORIZONTAL]
    assert h_face_set.size == 11

    base_face_set, = base.entity_variants(codimension=1)
    v_face_set, = [s for s in face_sets if s.variant_tag == Tag.VERTICAL]
    assert v_face_set.size == base_face_set.size * 10

    edge_sets = extrusion.entity_variants(codimension=2)
    h_edge_set, = [s for s in edge_sets if s.variant_tag == Tag.HORIZONTAL]
    assert h_edge_set.size == base_face_set.size * 11

    v_edge_set, = [s for s in edge_sets if s.variant_tag == Tag.VERTICAL]
    assert v_edge_set.size == base_vertex_set.size * 10
