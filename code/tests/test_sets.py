import pytest
import ufl

from meshstructure import (Index, IntervalEntitySet, TensorProductEntitySet,
                           TriangleEntitySet)


@pytest.mark.parametrize("lo", range(4))
@pytest.mark.parametrize("hi", range(4, 8))
def test_index_extent(lo, hi):
    i = Index(lo, hi)
    assert i.extent == hi - lo


def test_triangle_size():
    set_ = TriangleEntitySet(4, cell=ufl.triangle, codimension=0)
    assert len(set_.indices) == 2
    assert set_.size == 10


def test_tensor_product_size():
    a = TriangleEntitySet(4, cell=ufl.triangle, codimension=0)
    b = IntervalEntitySet(3, cell=ufl.interval, codimension=0)
    set_ = TensorProductEntitySet(a, b)
    assert len(set_.indices) == 3
    assert set_.size == 10*3
