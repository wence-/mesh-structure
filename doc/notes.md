# Notes on Interface for Structured Refinement

- consists of three parts:
  - entity sets
  - topology
  - geometry
- additional data layout

## Entity Sets

- polyhedral set
- property `isl_set`: set of inames and piecewise affine constraints
  - e.g {i, j: 0 <= i, j and i + j < N} or {i, j, k: 0 <= i, j, k < N}
- property `size`: number of points in set
- property `dimension/codim`: dimension/codimension of the entity it describes
- function `linear_index_map(index)` returns flatten index expression
  - e.g. for {i,j: 0<=i,j and i + j<N}: (N * (N-1)) // 2 - ((N - i) * (N - i - 1)) // 2 + i + j
  - e.g. for {i,j,k: 0<=i,j,k<N}: N**2 * k + N * j + i
- (maybe: knows the points on the boundary/interior of the polyhedral set)

**Note:**
- the pair `(index, entity_set)` is refered to as multiindex since an index needs to be
  associated with an entity set

## Topology

- property `base`: topology of base topology, e.g.
  - quadrilateral, simplex, etc
  - or topology with structured refinement
  - -> allows for extrusion of blockstructured mesh or vice versa
- property `refinement_rule`: description of local refinement
  - possible unnecessary since this information is also in its
	`entity_sets`
- property `dimension`: topological dimension
- function `entity_sets(dim/codim)` returns a list of entity sets with tags
  - tags describe different variants of entities with same reference element
	- **Note:** should these tags already be part of the entity sets?
  - e.g. one tag for vertical and one for horizontal edges of a 2d blockstructured quadrilateral
- function `cone(multiindex)` returns a list of multiindices with `codim + 1`
  - e.g. in case of a blockstructured quadrilateral grid
  `cone(((i,j), cells)) = [((i,j), hface), ((i,j+1), hface),
						   ((i,j), vface), ((i+1,j), vface)`
- function `support(multiindex)` returns a list of multiindices with `codim - 1` together with
  the local index of the input's entity
  - e.g. in case of extrusion of triangular grid with dune local numbering
  `support(((i,), hfaces)) = [(((i,), cells), 3), (((i - 1,), cells), 4)]`
- function `index_relation(multiindex, target_entity_set)` returns a list of multiindices
  with the same entity set as `target_entity_set` which are connected to `multiindex`
  - this has a generic implementation recursivly applying `cone` or `support`

### Generating cone and support relations

The points that are returned by the cone and support relations must be
consistently ordered with some concrete _reference element_ (for
example provided by FIAT, or DUNE).

The proposal is to do this like so:

Each entity set is augmented with a `reference_element` object which
provides integer offsets from a "zero" multiindex and the cone
(or supports). Applying the operation to a general multiindex then
consists of applying these integer offsets to the multiindex. We
additionally require a rule that maps a given multiindex on this
reference element to an entity number on the external package's
reference element. For example, a flattening rule for multiindices
into a single integer on simplices.

#### Generic rules for mapping from multiindices to a flat index

On hypercubes, we number multiindices lexicographically on the
`reference_element`, this gives us a consistent way of mapping into
the DUNE/FIAT element.

On simplices, we construct refinement patterns by performing Kuhn
refinement of a hypercube and discarding those parts of the hypercube
which are outside the simplex sum(i, j, k, ...) >= N.

To map the different types of simplex to a single oriented simplex
(transformation is different, orientation is easy), we order the
simplex vertices lexicographically by the hypercube multiindex.

So for example:

```
(0, 1)----(1, 1)
  |  \      |
  |   \     |
  |    \    |
  |     \   |
  |      \  |
(0, 0)----(0, 1)
```

Numbers the "left" triangle
```
  2
  |\
  | \
  |  \
  |   \
  |    \
  0-----1
```

And the "right" triangle
```
  1-----2
   \    |
    \   |
     \  |
      \ |
       \|
        0
```

Orienting edges from low to high results in consistently oriented
edges.

The same idea works for tets (and higher-dimensional simplices too,
but I can't draw them).

If, conversely, we want to keep the transformations the same (i.e. the
element tensor is a permutation), we split the triangles into two
types "left" and "right". "left" triangles are numbered
lexicographically in increasing multiindex order, "right" triangles
are numbered in reverse lexicographic order.

The same approach again works for the different classes of tets that
are produced in Kuhn refinement: they pair up and one is numbered
forwards, the other one backwards.

## Geometry

- property `topology`: underlying topology
- property `geometric_dimension`: dimension of the coordinate space
- property `topological_dimension`: dimension of topology
- property `finite_element`: defines the space of transformations from the reference element into the
  coordinate space
- function `geometry_dofs(multiindex)`: returns a 2d array of shape
  `(num_basis_functions, geometric_dimension)` with the coordinates of the dofs in the
  closure of the subentity described by `multiindex`
  - the closure is defined as the recursive application of `cone(multiindex)`
  - example blockstructured refinement in 2d: `w_k` corners of macro element,
	then the coordinates of the dofs for a cell with index `i,j` are
	`v_l = sum_k w_k phi_k(y_l)` with `y_l = (e_l + (i,j))/N` and `e_l` the
	corners of the reference element
- function `spatial_coordinate(multiindex, qp)`: returns the coordinates for a
  quadrature point given in the local coordinate system of the subentity
  described by `multiindex`
  - this can also be computed by using `geometry_dofs` but there are more
	efficient ways
  - example blockstructured refinement in 2d: `w_k` corners of macro element,
	then for a cell `i,j` the global coordinates of the quadrature point are
	`x = sum_k w_k phi_k((qp + (i,j)) / N)`
	- computation using `v_l = geometry_dofs(multiindex)`:
	  `x = sum_l v_l phi_l(qp)` + computation of `v_l`

## Data Layout

- knows underlying topology
- knows the #dofs per entity <- not clear if this should be exported
