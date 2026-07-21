# LoopyCuts Reproduction Modifications

This is an independent research reproduction snapshot based on:

https://github.com/mlivesu/LoopyCuts

Upstream baseline:

`c36b81154a03e79208f83725b9f4542f30ee4285`

Modified reproduction state:

`aa29c5327e6a7a47c4a0c1fa1753cf46a8b62627`

Modification date:

2026-07-17

This snapshot includes Linux/WSL compatibility changes, robustness
fixes, SHARP and loop-processing fixes, volumetric cutter changes,
build adjustments, and batch-processing tools.

This repository is not the official LoopyCuts repository.
Original copyright and license notices are retained.

## 2026-07-21 runtime robustness update

The following robustness fixes were added and tested on the `dancer`
and `des6` author models:

- Avoid querying an empty feature-edge octree during mesh smoothing.
  When the target surface contains no marked feature edges, feature
  vertices are reprojected onto the target surface instead.
- Preserve standard midpoint subdivision when a surface UV patch is
  unavailable.
- Preserve standard midpoint subdivision when a patch corner mapping
  is unavailable.
- Retain the CinoLib GCC compatibility fixes for ambiguous AABB
  initializer-list construction and renamed quality-update functions.

The subdivision fallback preserves the topology produced by midpoint
subdivision, but may locally differ from UV-patch-based surface
repositioning. Fallback events are explicitly reported with
`[SUBDIVISION_PATCH_FALLBACK]`.

Tested results:

- `dancer`: completed with 1090 hexahedra and 32 other polyhedra.
- `des6`: completed with 760 hexahedra and 8 other polyhedra.

These are robustness modifications to the reproduction snapshot and
are not part of the original upstream implementation.
