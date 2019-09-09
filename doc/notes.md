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

## Geometry

- property `topology`: underlying topology
- property `geometric_dimension`: dimension of the coordinate space
- property `topological_dimension`: dimension of topology
- property `finite_element`: defines the space of transformations from the reference element into the
  coordinate space
- function `geometry_dofs(multiindex)` returns a 2d array of shape
  `(num_basis_functions, geometric_dimension)` with the coordinates of the dofs in the
  closure of the subentity described by `multiindex`
  - the closure is defined as the recursive application of `cone(multiindex)`

## Data Layout

- knows underlying topology
- knows the #dofs per entity <- not clear if this should be exported
