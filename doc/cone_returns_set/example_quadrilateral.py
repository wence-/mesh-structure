import islpy as isl
import loopy as lp


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


cell_to_h_face = isl.BasicSet("[i,j] -> {[l,k]: l=i and j<=k<=j+1}")
h_face_to_vertex = isl.BasicSet("[l,k] -> {[o,p] : l<=o<=l+1 and p=k}")
cell_to_vertices_target = isl.BasicSet("[i,j] -> {[o,p] : i<=o<=i+1 and j<=p<=j+1}")

result = align(h_face_to_vertex, cell_to_h_face)
assert result == cell_to_vertices_target

cell_to_v_face = isl.BasicSet("[i,j] -> {[l2,k2]: i<=l2<i+2 and k2=j}")
v_face_to_vertex = isl.BasicSet("[l2,k2] -> {[o2,p2] : o2=l2 and k2<=p2<k2+2}")
result2 = align(v_face_to_vertex, cell_to_v_face)
assert result2 == result

elements = isl.BasicSet("{[i,j] : 0<=i,j<5}")
knl = lp.make_kernel([elements, cell_to_h_face, h_face_to_vertex], "a[o,p] = 1",
                     silenced_warnings=['inferred_iname'], lang_version=(2018,2))
print("\nKernel without aligned sets:")
print(lp.generate_body(knl))

knl = lp.make_kernel([elements, result], "a[o,p] = 1",
                     silenced_warnings=['inferred_iname'], lang_version=(2018,2))
print("\nKernel with aligned sets:")
print(lp.generate_body(knl))
