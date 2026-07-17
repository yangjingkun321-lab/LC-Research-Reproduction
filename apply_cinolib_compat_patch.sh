#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CINOLIB="$ROOT/lib/cinolib"
PATCH="$ROOT/patches/cinolib_gcc11_compat.patch"

if [[ ! -e "$CINOLIB/.git" ]]; then
    echo "[FAIL] CinoLib 子模块不存在。"
    echo "请先执行：git submodule update --init --recursive"
    exit 1
fi

if [[ ! -s "$PATCH" ]]; then
    echo "[FAIL] 补丁不存在或为空：$PATCH"
    exit 1
fi

if git -C "$CINOLIB" apply --reverse --check "$PATCH" \
    >/dev/null 2>&1
then
    echo "[PASS] CinoLib 兼容性补丁已经应用"
    exit 0
fi

git -C "$CINOLIB" apply --check "$PATCH"
git -C "$CINOLIB" apply "$PATCH"

echo "[PASS] CinoLib GCC 11 兼容性补丁已应用"
