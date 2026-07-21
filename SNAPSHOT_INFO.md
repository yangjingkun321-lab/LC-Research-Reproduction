# LoopyCuts Research Reproduction Snapshot

本仓库是用于科研复现、兼容性调试、代码检查和问题定位的
LoopyCuts 独立源码快照。

本仓库不是原作者的官方仓库，也不是 GitHub Fork。

## 原作者项目

仓库：

https://github.com/mlivesu/LoopyCuts

快照所采用的原作者版本：

`c36b81154a03e79208f83725b9f4542f30ee4285`

## 当前复现版本

本地工作分支：

`loopycuts-reproduction-fixes`

当前完整提交：

`aa29c5327e6a7a47c4a0c1fa1753cf46a8b62627`

主要稳定化提交：

`e5476986db0d98eac98f5f3ab32df5c9c608cddd`

批处理报告修正提交：

`aa29c5327e6a7a47c4a0c1fa1753cf46a8b62627`

## 子模块处理

本快照没有保留原项目和子模块的 Git 历史。

以下子模块在其记录提交处导出，并作为普通源码目录嵌入：

```text
 3d53718818480c99306d3d6830f49c3604d33ebd lib/cinolib (v1.0~205)
 36b95962756c1fce8e29b1f8bc45967f30773c00 lib/cinolib/external/eigen (3.3.0-1993-g36b959627)
 66376566852b704a0e57bf49dcac74ee5210ff18 lib/cinolib/external/graph_cut (heads/master)
 fb5258c157f9a9c40cbf4cc914dc1a7622865207 lib/vcg (v1.0.0-661-gfb5258c1)
```

因此，读取本快照不需要执行：

```bash
git submodule update --init --recursive
```
