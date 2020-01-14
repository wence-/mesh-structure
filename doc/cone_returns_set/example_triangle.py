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


# cell to horizontal, vertical, diagonal edge
hface = isl.BasicSet("[i,j] -> {[l,k]: l=i and k=j}")
vface = isl.BasicSet("[i,j] -> {[l,k]: l=i and k=j}")
dface = isl.BasicSet("[i,j] -> {[l,k]: l=i and k=j}")

# edge to vertex
hvertex = isl.BasicSet("[l,k] -> {[o,p]: l<=o<l+2 and p=k}")
vvertex = isl.BasicSet("[l,k] -> {[o,p]: k<=p<k+2 and o=l}")
dvertex = isl.Set("[l,k] -> {[o,p]: (o=l and p=k+1) or (o=l+1 and p=k)}")

h_aligned = align(hvertex, hface)
v_aligned = align(vvertex, vface)
d_aligned = align(dvertex, dface)

u = h_aligned | v_aligned
u = u | d_aligned
print("Cells to vertices:")
print("    ", u)

u = h_aligned | (v_aligned - (h_aligned & v_aligned))
u = u | (d_aligned - (d_aligned & u))
print("Cells to vertices (reduced):")
print("    ", u)

realized = realize_params(u, [0, 0])
print("Vertices of (lower) cell [0, 0]:")
realized.foreach_point(lambda p: print('    ', p))

realized = realize_params(u, [2, 1])
print("Vertices of (lower) cell [2, 1]:")
realized.foreach_point(lambda p: print('    ', p))
