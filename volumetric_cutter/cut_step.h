#ifndef CUT_STEP_H
#define CUT_STEP_H

#include "state.h"

// Result category of one Stage-2 loop action.
enum class CutStepStatus
{
    NOOP_USED,
    NOOP_REVERTED,
    NOOP_NICO_BUG,
    COMMITTED,
    REVERTED,
    UNKNOWN
};

struct CutStepResult
{
    uint loop_id = 0;

    CutStepStatus status = CutStepStatus::UNKNOWN;

    bool attempted = false;
    bool committed = false;
    bool reverted = false;
    bool converged = false;

    uint verts_before = 0;
    uint verts_after = 0;

    uint tets_before = 0;
    uint tets_after = 0;

    //
    // Diagnostics of the MeshExtractor built after this step.
    //
    uint nonmanifold_polys = 0;
    uint high_genus_polys = 0;
    uint buggy_chains = 0;

    double step_time = 0.0;
};

// Execute exactly one LoopyCuts Stage-2 loop action.
//
// This function:
//   1. inspects the current loop state,
//   2. calls the original cut(),
//   3. rebuilds the current meta mesh,
//   4. tests convergence,
//   5. returns compact step information.
//
// It does NOT:
//   - choose which loop to execute,
//   - finalize the decomposition,
//   - save output meshes,
//   - maintain RL replay/history.
CutStepResult execute_cut_step(GlobalState &state,
                               uint loop_id);

const char *cut_step_status_string(CutStepStatus status);

#endif // CUT_STEP_H