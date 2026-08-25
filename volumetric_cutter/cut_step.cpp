#include "cut_step.h"

#include "cut.h"
#include "mesh_extractor.h"

#include <cassert>
#include <chrono>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

const char *cut_step_status_string(CutStepStatus status)
{
    switch (status)
    {
    case CutStepStatus::NOOP_USED:
        return "NOOP_USED";

    case CutStepStatus::NOOP_REVERTED:
        return "NOOP_REVERTED";

    case CutStepStatus::NOOP_NICO_BUG:
        return "NOOP_NICO_BUG";

    case CutStepStatus::COMMITTED:
        return "COMMITTED";

    case CutStepStatus::REVERTED:
        return "REVERTED";

    case CutStepStatus::UNKNOWN:
    default:
        return "UNKNOWN";
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CutStepResult execute_cut_step(GlobalState &state,
                               uint loop_id)
{
    assert(loop_id < state.loops.size());

    CutStepResult result;

    result.loop_id = loop_id;

    Loop &loop = state.loops.at(loop_id);

    //
    // State before calling the original LoopyCuts cut().
    //
    const bool pre_used =
        loop.used;

    const bool pre_reverted =
        loop.reverted;

    const bool pre_nico_bug =
        loop.Nico_bug;

    result.verts_before =
        state.m_vol.num_verts();

    result.tets_before =
        state.m_vol.num_polys();

    //
    // Preserve the classification semantics currently used by
    // run_external_order().
    //
    if (pre_used)
    {
        result.status =
            CutStepStatus::NOOP_USED;
    }
    else if (pre_reverted)
    {
        result.status =
            CutStepStatus::NOOP_REVERTED;
    }
    else if (pre_nico_bug)
    {
        result.status =
            CutStepStatus::NOOP_NICO_BUG;
    }
    else
    {
        result.attempted = true;
    }

    const auto step_begin =
        std::chrono::steady_clock::now();

    //
    // ORIGINAL LoopyCuts geometry operation.
    //
    cut(state, loop_id);

    //
    // ORIGINAL Stage-2 meta-mesh extraction and convergence test.
    //
    state.profiler.push("Mesh Extractor");

    MeshExtractor me(state.m_vol);

    state.m_poly = me.mm;

    result.converged =
        me.converged(state.m_vol);

    //
    // MeshExtractor has already computed these values.
    // Reading them here adds no second extraction pass.
    //
    result.nonmanifold_polys =
        me.nonmanifold_poly_count();

    result.high_genus_polys =
        me.high_genus_poly_count();

    result.buggy_chains =
        me.buggy_chain_count();

    state.profiler.pop("\n");

    const auto step_end =
        std::chrono::steady_clock::now();

    result.step_time =
        std::chrono::duration<double>(
            step_end - step_begin)
            .count();

    //
    // Determine the result of a real cut attempt.
    //
    if (result.attempted)
    {
        if (loop.reverted)
        {
            result.reverted = true;

            result.status =
                CutStepStatus::REVERTED;
        }
        else if (loop.used)
        {
            result.committed = true;

            result.status =
                CutStepStatus::COMMITTED;
        }
        else
        {
            result.status =
                CutStepStatus::UNKNOWN;
        }
    }

    result.verts_after =
        state.m_vol.num_verts();

    result.tets_after =
        state.m_vol.num_polys();

    return result;
}