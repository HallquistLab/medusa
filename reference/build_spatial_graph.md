# Build a sparse spatial neighborhood graph for bilateral HRF and noise pooling

Constructs a compressed-sparse-row (CSR) spatial neighborhood graph
(`"medusa_spatial_graph"`) for
[`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)
when `pooling = "graph"`.

## Usage

``` r
build_spatial_graph(
  x,
  spatial_fwhm = NULL,
  max_neighbors = 32L,
  max_radius = NULL,
  spatial_groups = NULL,
  parcel_boundary = c("soft", "hard", "ignore"),
  metric = c("auto", "geodesic", "euclidean"),
  n_voxels = NULL,
  voxel_ids = NULL,
  x_type = c("auto", "coords", "distance")
)
```

## Arguments

- x:

  Spatial input: a `"medusa_nifti_ts"` object, a list containing
  `voxel_map` (and optional `geometry`), a `voxel_map` `data.frame` with
  lattice columns `i, j, k` (or world coordinates `x, y, z`), a numeric
  \\V \times D\\ coordinate matrix, an edge `data.frame` with columns
  `from`, `to` (and optional positive `dist` together with `n_voxels`),
  a [`dist`](https://rdrr.io/r/stats/dist.html) object, a square numeric
  or sparse/dense `Matrix` distance matrix (positive finite off-diagonal
  entries giving neighbor distances; `0` or `Inf` denoting disconnected
  pairs), or an existing `"medusa_spatial_graph"`.

- spatial_fwhm:

  Optional positive numeric Gaussian spatial kernel FWHM (in mm or
  coordinate units). Defaults to `2.5 * median(spacing)` for lattice
  inputs and `2.5` times the median nearest-neighbor distance for
  coordinate and distance matrices.

- max_neighbors:

  Maximum number of non-self spatial neighbors retained per node
  (default `32L`). Every node always retains itself at index 1 of its
  CSR row with distance `0`.

- max_radius:

  Optional positive search radius cutoff. Defaults to
  `max(2.75 * sigma, 1.80 * step_scale)` where
  `sigma = spatial_fwhm / sqrt(8 * log(2))`.

- spatial_groups:

  Optional character or integer vector of parcel/group labels of length
  \\V\\. Only restricts graph construction when
  `parcel_boundary = "hard"`; under `"soft"` or `"ignore"`, cross-border
  neighbors remain in the graph so bilateral pooling can modulate
  cross-border weights cleanly.

- parcel_boundary:

  Boundary handling when `spatial_groups` is present: `"soft"`
  (default), `"hard"` (omit all cross-parcel edges and geodesic paths),
  or `"ignore"`.

- metric:

  Distance metric for volumetric/coordinate inputs: `"auto"` (uses
  26-connected in-mask lattice geodesic distance whenever unique integer
  lattice indices `i, j, k` are available, so out-of-mask CSF/sulcal
  gaps act as natural barriers, and falls back to 3D Euclidean
  radius/\\k\\-NN otherwise), `"geodesic"`, or `"euclidean"`.

- n_voxels:

  Optional total node count \\V\\ when `x` is an edge `data.frame` that
  may contain isolated nodes.

- voxel_ids:

  Optional character vector of node IDs of length \\V\\.

- x_type:

  Disambiguation for matrix inputs `x`: `"auto"` (default; treats `dist`
  and `Matrix` inputs, as well as square \\V \times V\\ numeric matrices
  with \\V \> 3\\, zero diagonal, and nonnegative entries as pairwise
  distance matrices, and other numeric matrices as \\V \times D\\
  coordinate matrices), `"coords"` (force interpretation as a \\V \times
  D\\ coordinate matrix), or `"distance"` (force interpretation as a
  square \\V \times V\\ distance matrix, including small \\2 \times 2\\
  or \\3 \times 3\\ distance matrices).

## Value

An S3 object of class `"medusa_spatial_graph"` containing 1-based CSR
arrays `row_ptr` (length \\V + 1\\), `col_idx`, and `dist`, along with
`n_voxels`, `voxel_ids`, `spatial_fwhm`, `sigma`, `max_radius`,
`max_neighbors`, `metric`, and `parcel_boundary`.
