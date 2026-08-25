# LoopyCuts V4 Stage-2 Runtime Provenance

This record freezes the known source state corresponding to the
LoopyCuts Stage-2 RL runtime used before V5 development.

## Parent LoopyCuts base

- branch at freeze: `loopycuts-reproduction-fixes`
- pre-freeze parent HEAD:
  `03bcff8b09ce6c03b75a1ee59e43bc5aeac3730c`

## Frozen V4 executable

Path at freeze:

`volumetric_cutter/volumetric_cutter`

SHA256:

`920ff358d7840dfa13dc53cc0b603721347cc77e8951afb14d322de14f84a565`

Recorded size:

`2123976 bytes`

Recorded mtime:

`2026-08-13 09:52:47.619529300 +0800`

## CinoLib

Parent gitlink:

`3d53718818480c99306d3d6830f49c3604d33ebd`

The V4 runtime source tree contains local CinoLib compatibility/runtime
fixes on top of that exact upstream commit.

The patch is versioned as:

`provenance/v4/cinolib_v4_runtime.patch`

Expected patch SHA256:

`39043fd06bc0f425588792e4e108f2c555f6fcb42c5918c42439d73192573d93`

Apply/verify it with:

`provenance/v4/apply_cinolib_v4_runtime_patch.sh`

Modified CinoLib source files:

- `include/cinolib/meshes/hexmesh.cpp`
- `include/cinolib/octree.cpp`
- `include/cinolib/smoother.cpp`

All three current modified source files have mtimes earlier than the
frozen executable.

This temporal evidence is consistent with their participation in the
V4 build, but mtime evidence alone is not claimed as a cryptographic
source-to-binary proof.

## Stage-2 V4 runtime source

Tracked runtime modifications:

- `volumetric_cutter/batch.cpp`
- `volumetric_cutter/batch.h`
- `volumetric_cutter/cut.cpp`
- `volumetric_cutter/finalization.cpp`
- `volumetric_cutter/main.cpp`
- `volumetric_cutter/mesh_extractor.h`
- `volumetric_cutter/subdivision_helper.cpp`
- `volumetric_cutter/volumetric_cutter.pro`

New runtime source:

- `volumetric_cutter/cut_step.cpp`
- `volumetric_cutter/cut_step.h`
- `volumetric_cutter/rl_server.cpp`
- `volumetric_cutter/rl_server.h`

## Scope

This freeze commit deliberately excludes:

- Stage-1 loop-distribution reproduction changes
- generated AntTweakBar objects/libraries
- build directories and logs
- experimental patch scripts
- VS Code configuration
- generated test meshes
- temporary `_tet.mesh` / `_hex.mesh`
- CinoLib backup files

Stage-1 reproduction changes are to be frozen separately.

## External full-worktree snapshot

Before this freeze commit, the complete working tree excluding `.git`
was archived externally.

Snapshot SHA256:

`5391b133a13d3e59a95bd04af5270e56e3403d28975ecb7c0437b49e9d278d58`

The snapshot is a recovery artifact, not a substitute for the
versioned source provenance above.
