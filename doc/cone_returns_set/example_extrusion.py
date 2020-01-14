import islpy as isl


def align(minor, major, verbose=False):
    major_dims = len(major.get_var_names(isl.dim_type.set))
    aligned1, aligned2 = isl.align_two(minor, major)
    if verbose:
        print(aligned1)
        print(aligned2)
    aligned = aligned1 & aligned2
    if verbose:
        print(aligned)
    return aligned.remove_dims(isl.dim_type.set, 0, major_dims)


def realize_params(set, params):
    n_params = len(set.get_var_names(isl.dim_type.param))
    assert len(params) <= n_params
    realized_set = set
    for i, p in enumerate(params):
        realized_set = realized_set.fix_val(isl.dim_type.param, i, p)
    return realized_set.remove_dims(isl.dim_type.param, 0, len(params))


# cells to faces
cell_to_h_face = isl.BasicSet("[i] -> {[l]: i<=l<i+2}")
cell_to_v_face = isl.BasicSet("[i] -> {[l, k]: l=i and 0<=k<3}")

# faces to edges
h_to_h_edge = isl.BasicSet("[l] -> {[o, p]: o=l and 0<=p<3}")

v_to_h_edge = isl.BasicSet("[l, k] -> {[o, p]: l<=o<l+2 and p=k}")
# BasicSet cannot handle `or` -> use Set
# not sure if Set can be used as loop domain for loopy, maybe need to use get_basic_sets()
v_to_v_edge = isl.Set("[l, k] -> {[o, p]: o=l and 0<=p<3 and "
                      "(p = ((k + 1) mod 3) or p = (k mod 3))}")

# vertices
v_to_vertices = isl.BasicSet("[o, p] -> {[u, v]: o<=u<o+2 and v=p}")
h_to_vertices = isl.Set("[o, p] -> {[u, v]: u=o and 0<=v<3 and "
                        "(v = ((p + 1) mod 3) or v = (p mod 3))}")

# show that the constructions using mod are correct
assert realize_params(v_to_v_edge, [0, 0]) == isl.Set("{[o, p]: o=0 and (p=0 or p=1)}")
assert realize_params(v_to_v_edge, [0, 2]) == isl.Set("{[o, p]: o=0 and (p=2 or p=0)}")
assert realize_params(h_to_vertices, [0, 2]) == isl.Set("{[u, v]: u=0 and (v=2 or v=0)}")

all_h_edges_extr = align(h_to_h_edge, cell_to_h_face)
all_h_edges_extr |= align(v_to_h_edge, cell_to_v_face)
print("Cells to horizontal edges:")
print("    ", all_h_edges_extr)

all_v_edges_extr = align(v_to_v_edge, cell_to_v_face)
print("Cells to vertical edges:")
print("    ", all_v_edges_extr)

all_vertices = align(v_to_vertices, all_v_edges_extr)
all_vertices |= align(h_to_vertices, all_h_edges_extr)
print("Cells to vertices:")
print("    ", all_vertices)