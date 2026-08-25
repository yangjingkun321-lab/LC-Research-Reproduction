#!/usr/bin/env bash

set -euo pipefail

ROOT="$(
    cd "$(
        dirname "${BASH_SOURCE[0]}"
    )/../.."
    pwd
)"

SUBMODULE="$ROOT/lib/cinolib"

PATCH="$ROOT/provenance/v4/cinolib_v4_runtime.patch"

EXPECTED_CINOLIB_HEAD="3d53718818480c99306d3d6830f49c3604d33ebd"

EXPECTED_PATCH_SHA256="39043fd06bc0f425588792e4e108f2c555f6fcb42c5918c42439d73192573d93"


observed_head="$(
    git -C "$SUBMODULE" rev-parse HEAD
)"

if test "$observed_head" != "$EXPECTED_CINOLIB_HEAD"; then

    echo \
        "ERROR: unexpected CinoLib HEAD:" \
        "$observed_head" \
        >&2

    exit 1
fi


observed_patch_sha="$(
    sha256sum "$PATCH" |
    awk '{print $1}'
)"

if test "$observed_patch_sha" != "$EXPECTED_PATCH_SHA256"; then

    echo \
        "ERROR: CinoLib patch SHA256 mismatch:" \
        "$observed_patch_sha" \
        >&2

    exit 1
fi


if git -C "$SUBMODULE" apply \
    --reverse \
    --check \
    "$PATCH" \
    >/dev/null 2>&1
then

    echo "CINOLIB V4 PATCH ALREADY APPLIED"

    exit 0
fi


if ! git -C "$SUBMODULE" diff --quiet; then

    echo \
        "ERROR: CinoLib has non-matching tracked modifications" \
        >&2

    git -C "$SUBMODULE" status --short >&2

    exit 1
fi


git -C "$SUBMODULE" apply \
    --check \
    "$PATCH"

git -C "$SUBMODULE" apply \
    "$PATCH"


if ! git -C "$SUBMODULE" apply \
    --reverse \
    --check \
    "$PATCH" \
    >/dev/null 2>&1
then

    echo \
        "ERROR: patch application verification failed" \
        >&2

    exit 1
fi


echo "CINOLIB V4 PATCH APPLY PASS"
