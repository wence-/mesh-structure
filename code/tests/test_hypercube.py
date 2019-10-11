import pytest

from meshstructure import (HyperCubeRefinement)
from meshstructure.unstructured import UnstructuredHyperCube


@pytest.mark.parametrize('dim', range(1, 4))
def test_hypercube_size(dim):
    h = HyperCubeRefinement(UnstructuredHyperCube(dim), *([5]*dim))
    cell_set, = h.entity_variants(codimension=0)
    assert cell_set.size == 5**dim

    vertex_set, = h.entity_variants(codimension=dim)
    assert vertex_set.size == 6**dim

    for codim in range(1,dim):
        sets = h.entity_variants(codimension=codim)
        set_sizes = [set.size for set in sets]
        assert all(set_sizes[0] == s for s in set_sizes)
        assert set_sizes[0] == 6**codim * 5**(dim-codim)
