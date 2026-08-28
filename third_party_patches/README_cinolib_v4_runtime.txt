LoopyCuts V5 inherited cinolib runtime patch
===========================================

Submodule:
    lib/cinolib

Required base commit:
    3d53718818480c99306d3d6830f49c3604d33ebd

Patch:
    third_party_patches/cinolib_v4_runtime.patch

Patch SHA256:
    39043fd06bc0f425588792e4e108f2c555f6fcb42c5918c42439d73192573d93

Patched files:
    include/cinolib/meshes/hexmesh.cpp
    include/cinolib/octree.cpp
    include/cinolib/smoother.cpp

This patch is inherited byte-for-byte from the frozen V4 runtime.
It is not a V5 algorithm change.

Apply from repository root:

    git -C lib/cinolib checkout 3d53718818480c99306d3d6830f49c3604d33ebd
    git -C lib/cinolib apply ../../third_party_patches/cinolib_v4_runtime.patch

The M14-validated V5 executable SHA256 is:

    0adfbb90e86d1166f6b85bd5c21f1b8bf39d1b808cf3fb80a1b4092e35bcb846
