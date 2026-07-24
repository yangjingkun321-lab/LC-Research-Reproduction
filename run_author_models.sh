#!/usr/bin/env bash
set -u

ROOT="${LOOPYCUTS_ROOT:-$HOME/codes/LoopyCuts}"
DATA="$ROOT/test_data"
RUNROOT="${LOOPYCUTS_RUNROOT:-$HOME/loopycuts_runs}"

LD_BIN="$ROOT/loop_distribution/loop_distributor"
VC_BIN="$ROOT/volumetric_cutter/volumetric_cutter"

mkdir -p "$RUNROOT"

die()
{
    echo "[FAIL] $*" >&2
    exit 1
}

test -d "$DATA" || die "找不到作者数据目录：$DATA"
test -x "$LD_BIN" || die "找不到：$LD_BIN"
test -x "$VC_BIN" || die "找不到：$VC_BIN"


find_stem()
{
    local src="$1"
    local setup
    local stem
    local -a valid=()

    while IFS= read -r -d '' setup
    do
        stem="${setup%_setup.txt}"

        if [[ \
            -s "${stem}.obj" &&
            -s "${stem}.rosy" &&
            -s "${stem}.sharp"
        ]]
        then
            valid+=("$(basename "$stem")")
        fi
    done < <(
        find "$src" \
          -maxdepth 1 \
          -type f \
          -name '*_setup.txt' \
          -print0 \
        | sort -z
    )

    if (( ${#valid[@]} == 0 ))
    then
        return 1
    fi

    if (( ${#valid[@]} > 1 ))
    then
        echo \
          "[WARN] $(basename "$src") 有多个正式 stem，默认使用 ${valid[0]}：${valid[*]}" \
          >&2
    fi

    printf '%s\n' "${valid[0]}"
}


list_cases()
{
    local src
    local case_name
    local stem
    local state

    printf 'CASE\tSTEM\tRUN_STATE\tSOURCE\n'

    while IFS= read -r -d '' src
    do
        case_name="$(basename "$src")"

        if ! stem="$(find_stem "$src")"
        then
            printf \
              '%s\t-\tNO_FORMAL_INPUT\t%s\n' \
              "$case_name" \
              "$src"

            continue
        fi

        if [[ -d "$RUNROOT/${case_name}_clean" ]]
        then
            state="EXISTS"
        else
            state="NOT_RUN"
        fi

        printf \
          '%s\t%s\t%s\t%s\n' \
          "$case_name" \
          "$stem" \
          "$state" \
          "$src"

    done < <(
        find "$DATA" \
          -mindepth 1 \
          -maxdepth 1 \
          -type d \
          -print0 \
        | sort -z
    )
}


write_report()
{
    local run="$1"
    local status="$2"
    local ld_exit="$3"
    local vc_exit="$4"
    local loops="$5"
    local stats="$6"
    local kind="$7"
    local infnan="$8"
    local maxrss="$9"
    local zero="${10}"
    local dead="${11}"
    local revert="${12}"
    local recovery="${13}"

    cat > "$run/run_report.txt" <<EOF
CASE=$CASE
SOURCE=$SRC
STEM=$STEM
FINAL_STATUS=$status
LOOP_DISTRIBUTOR_EXIT=$ld_exit
VOLUMETRIC_CUTTER_EXIT=$vc_exit
CUTTING_LOOPS=$loops
FINAL_CELL_COUNTS=$stats
OUTPUT_KIND=$kind
INF_NAN_FOUND=$infnan
MAX_RSS_KB=$maxrss
ZERO_VECTOR_WARNINGS=$zero
DEAD_END_WARNINGS=$dead
REVERT_CUT_COUNT=$revert
UV_DEGENERATE_RECOVERIES=$recovery
EOF
}


run_case()
{
    CASE="$1"
    local force="${2:-0}"

    SRC="$DATA/$CASE"

    if [[ ! -d "$SRC" ]]
    then
        echo "[FAIL] 不存在作者模型目录：$SRC"
        return 0
    fi

    if ! STEM="$(find_stem "$SRC")"
    then
        echo \
          "[FAIL_INPUT] $CASE 没有完整的 .obj/.rosy/.sharp/_setup.txt 同前缀输入"

        return 0
    fi

    local run="$RUNROOT/${CASE}_clean"
    local stamp

    if [[ -e "$run" ]]
    then
        if [[ "$force" == "1" ]]
        then
            stamp="$(date +%Y%m%d_%H%M%S)"

            mv \
              "$run" \
              "${run}.previous_${stamp}"

            echo "[BACKUP] ${run}.previous_${stamp}"
        else
            echo "[SKIP_EXISTING] $run"
            return 0
        fi
    fi

    echo
    echo "================================================================"
    echo "[START] CASE=$CASE STEM=$STEM"
    echo "================================================================"

    mkdir -p "$run"

    cp -a \
      "$SRC/$STEM.obj" \
      "$SRC/$STEM.rosy" \
      "$SRC/$STEM.sharp" \
      "$SRC/${STEM}_setup.txt" \
      "$run/"

    local obj="$run/$STEM.obj"
    local rosy="$run/$STEM.rosy"
    local sharp="$run/$STEM.sharp"
    local setup="$run/${STEM}_setup.txt"

    local loop="$run/${STEM}_loop.txt"
    local split="$run/${STEM}_splitted.obj"

    local out="$run/output"
    local ldlog="$run/loop_distributor_batch.log"
    local vclog="$run/volumetric_cutter_batch.log"

    # test_data 可能包含作者已经生成的输出。
    # 必须删除运行副本中的旧输出，防止第一阶段失败后误用旧文件。
    rm -f \
      "$loop" \
      "$split" \
      "$run/${STEM}_console.txt" \
      "$ldlog" \
      "$vclog"

    rm -rf "$out"
    mkdir -p "$out"

    local faces
    local rosy_faces
    local sharp_header
    local sharp_records

    faces="$(
        grep -cE \
          '^[[:space:]]*f[[:space:]]' \
          "$obj" \
        || true
    )"

    rosy_faces="$(
        awk \
          'NF && $1 !~ /^#/ {print $1; exit}' \
          "$rosy"
    )"

    sharp_header="$(
        awk \
          'NF && $1 !~ /^#/ {print $1; exit}' \
          "$sharp"
    )"

    sharp_records="$(
        grep -cve \
          '^[[:space:]]*$' \
          "$sharp"
    )"

    sharp_records=$((sharp_records - 1))

    {
        echo "OBJ=$obj"
        echo "OBJ_FACES=$faces"
        echo "ROSY_HEADER=$rosy_faces"
        echo "SHARP_HEADER=$sharp_header"
        echo "SHARP_RECORDS=$sharp_records"
        echo "SETUP=$setup"
    } | tee "$run/input_check.txt"

    if [[ \
        ! "$faces" =~ ^[0-9]+$ ||
        ! "$rosy_faces" =~ ^[0-9]+$ ||
        "$faces" -ne "$rosy_faces"
    ]]
    then
        echo "[FAIL_INPUT] OBJ 面数与 ROSY 头部不一致"

        write_report \
          "$run" \
          "FAIL_INPUT_ROSY" \
          "-" "-" "-" "-" \
          "NONE" "-" "-" "-" "-" "-" "-"

        return 0
    fi

    if [[ \
        ! "$sharp_header" =~ ^[0-9]+$ ||
        "$sharp_header" -ne "$sharp_records"
    ]]
    then
        echo "[FAIL_INPUT] SHARP 头部与记录数不一致"

        write_report \
          "$run" \
          "FAIL_INPUT_SHARP" \
          "-" "-" "-" "-" \
          "NONE" "-" "-" "-" "-" "-" "-"

        return 0
    fi

    echo "[STAGE 1] loop_distribution"

    cd "$ROOT/loop_distribution" || return 1

    set +e
    set -o pipefail

    /usr/bin/time -v \
      "$LD_BIN" \
      "$obj" \
      batch \
      2>&1 | tee "$ldlog"

    local ld_exit=${PIPESTATUS[0]}

    set +o pipefail

    echo "loop_distributor exit=$ld_exit" \
      | tee -a "$ldlog"

    if [[ \
        "$ld_exit" -ne 0 ||
        ! -s "$loop" ||
        ! -s "$split"
    ]]
    then
        echo "[FAIL_STAGE1] $CASE exit=$ld_exit"

        write_report \
          "$run" \
          "FAIL_LOOP_DISTRIBUTION" \
          "$ld_exit" "-" "-" "-" \
          "NONE" "-" "-" "-" "-" "-" "-"

        return 0
    fi

    local loops

    loops="$(
        awk \
          'NF {print $1; exit}' \
          "$loop"
    )"

    if [[ ! "$loops" =~ ^[0-9]+$ ]]
    then
        loops="INVALID"
    fi

    echo "[STAGE 1 PASS] cutting_loops=$loops"
    echo "[STAGE 2] volumetric_cutter"

    cd "$ROOT/volumetric_cutter" || return 1

    set +e
    set -o pipefail

    /usr/bin/time -v \
      "$VC_BIN" \
      "$split" \
      "$loop" \
      -batch-mode "$out" \
      2>&1 | tee "$vclog"

    local vc_exit=${PIPESTATUS[0]}

    set +o pipefail

    echo "volumetric_cutter exit=$vc_exit" \
      | tee -a "$vclog"

    local smoothed
    local hexmesh

    smoothed="$(
        find "$out" \
          -maxdepth 1 \
          -type f \
          -name '*_mm_subdivided_smoothed.hedra' \
          -size +0c \
          -print \
        | head -n 1
    )"

    hexmesh="$(
        find "$out" \
          -maxdepth 1 \
          -type f \
          -name '*_hex.mesh' \
          -size +0c \
          -print \
        | head -n 1
    )"

    local infnan=0
    local file

    while IFS= read -r -d '' file
    do
        if grep -qiE \
          '(^|[[:space:]])[-+]?(inf|nan)([[:space:]]|$)' \
          "$file"
        then
            infnan=1
            break
        fi
    done < <(
        find "$out" \
          -maxdepth 1 \
          -type f \
          -name '*.hedra' \
          -print0
    )

    local stats

    stats="$(
        python3 - "$vclog" <<'PY_CELL_COUNTS'
from pathlib import Path
import re
import sys

text = Path(sys.argv[1]).read_text(
    encoding="utf-8",
    errors="replace",
)

# 删除终端颜色和其他 ANSI CSI 控制序列。
text = re.sub(
    r"\x1b\[[0-?]*[ -/]*[@-~]",
    "",
    text,
)

matches = re.findall(
    r"\d+\s+hexa\s+-\s+"
    r"\d+\s+prisms\s+-\s+"
    r"\d+\s+hexable with midpoint\s+-\s+"
    r"\d+\s+others",
    text,
)

print(matches[-1] if matches else "NOT_FOUND")
PY_CELL_COUNTS
    )"

    local maxrss

    maxrss="$(
        sed -nE \
          's/^[[:space:]]*Maximum resident set size \(kbytes\):[[:space:]]*([0-9]+).*$/\1/p' \
          "$vclog" \
        | tail -n 1
    )"

    if [[ -z "$maxrss" ]]
    then
        maxrss="-"
    fi

    local zero
    local dead
    local revert
    local recovery

    zero="$(
        grep -cF \
          'WARNING: normalization of zero length vector!' \
          "$vclog" \
        || true
    )"

    dead="$(
        grep -cF \
          'WARNING: removing a dead end has generated a valence 2 vertex. FixMe' \
          "$vclog" \
        || true
    )"

    revert="$(
        grep -cF \
          'Revert cut' \
          "$vclog" \
        || true
    )"

    recovery="$(
        grep -cF \
          'SUBDIVISION_DEGENERATE_UV_TRIANGLE_RECOVERED' \
          "$vclog" \
        || true
    )"

    local status="FAIL_VOLUMETRIC_CUTTER"
    local kind="NONE"

    if [[ \
        "$vc_exit" -eq 0 &&
        -n "$smoothed" &&
        "$infnan" -eq 0
    ]]
    then
        if [[ -n "$hexmesh" ]]
        then
            status="PASS_HEX"
            kind="PURE_HEX_MESH"
        else
            status="PASS_HYBRID"
            kind="HYBRID_HEDRA"
        fi
    elif [[ "$infnan" -ne 0 ]]
    then
        status="FAIL_INF_NAN"
    fi

    write_report \
      "$run" \
      "$status" \
      "$ld_exit" \
      "$vc_exit" \
      "$loops" \
      "$stats" \
      "$kind" \
      "$infnan" \
      "$maxrss" \
      "$zero" \
      "$dead" \
      "$revert" \
      "$recovery"

    echo "[END] $CASE -> $status"
    cat "$run/run_report.txt"
}


show_status()
{
    local report

    printf \
      'CASE\tSTATUS\tLOOPS\tFINAL_CELL_COUNTS\tOUTPUT_KIND\tMAX_RSS_KB\n'

    while IFS= read -r -d '' report
    do
        awk -F= '
            /^CASE=/              {case_name=$2}
            /^FINAL_STATUS=/      {status=$2}
            /^CUTTING_LOOPS=/     {loops=$2}
            /^FINAL_CELL_COUNTS=/ {counts=$2}
            /^OUTPUT_KIND=/       {kind=$2}
            /^MAX_RSS_KB=/        {rss=$2}

            END
            {
                printf \
                    "%s\t%s\t%s\t%s\t%s\t%s\n", \
                    case_name, status, loops, counts, kind, rss
            }
        ' "$report"

    done < <(
        find "$RUNROOT" \
          -mindepth 2 \
          -maxdepth 2 \
          -type f \
          -name run_report.txt \
          ! -path '*.previous_*/*' \
          -print0 \
        | sort -z
    )
}


run_all()
{
    local src
    local case_name
    local stem

    while IFS= read -r -d '' src
    do
        case_name="$(basename "$src")"

        if stem="$(find_stem "$src")"
        then
            run_case "$case_name" 0
        else
            echo "[SKIP_NO_FORMAL_INPUT] $case_name"
        fi

    done < <(
        find "$DATA" \
          -mindepth 1 \
          -maxdepth 1 \
          -type d \
          -print0 \
        | sort -z
    )
}


case "${1:-}" in

    list)
        list_cases
        ;;

    status)
        show_status \
          | tee "$RUNROOT/author_models_summary.tsv"
        ;;

    run)
        [[ $# -eq 2 ]] \
          || die "用法：$0 run <CASE>"

        run_case "$2" 0
        ;;

    rerun)
        [[ $# -eq 2 ]] \
          || die "用法：$0 rerun <CASE>"

        run_case "$2" 1
        ;;

    all)
        run_all

        show_status \
          | tee "$RUNROOT/author_models_summary.tsv"
        ;;

    *)
        cat <<EOF
用法：
  $0 list
  $0 run <CASE>
  $0 rerun <CASE>
  $0 all
  $0 status
EOF
        exit 2
        ;;

esac
