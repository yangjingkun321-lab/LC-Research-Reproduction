#include "rl_server.h"

#include "cut_step.h"
#include "finalization.h"
#include "subdivision_helper.h"
#include "mesh_smoother.h"
#include "polyhedral_decomposition.h"

#include <cinolib/export_hexahedra.h>
#include <cinolib/string_utilities.h>

#include <iostream>
#include <sstream>
#include <string>

extern std::string model_name;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Convert LoopyCuts loop type to a stable string for the RL protocol.
//
static const char *rl_loop_type_string(const Loop &loop)
{
    if (loop.type == CONCAVE)
    {
        return "CONCAVE";
    }

    if (loop.type == CONVEX)
    {
        return "CONVEX";
    }

    //
    // Current LoopyCuts loop files use:
    // CONCAVE / REGULAR / CONVEX.
    //
    return "REGULAR";
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Base validity of a volumetric loop action.
//
// This does NOT consider convergence or the RL phase.
// It only checks whether the loop itself is still usable.
//
static bool is_base_action_valid(const GlobalState &state,
                                 const uint loop_id)
{
    if (loop_id >= state.loops.size())
    {
        return false;
    }

    const Loop &loop =
        state.loops.at(loop_id);

    if (loop.type == CONVEX)
    {
        return false;
    }

    if (loop.used)
    {
        return false;
    }

    if (loop.reverted)
    {
        return false;
    }

    if (loop.Nico_bug)
    {
        return false;
    }

    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// RL-V1 Stage-2 action semantics.
//
// Before the first convergence event:
//
//     valid CONCAVE + REGULAR loops are selectable.
//
// Once convergence has been reached at least once:
//
//     regular_phase_closed = true
//
// and this flag NEVER becomes false again.
//
// Therefore:
//
//     REGULAR loops are permanently closed.
//
//     Only still-valid CONCAVE loops may remain selectable.
//
// IMPORTANT:
//
// "converged" and "regular_phase_closed" are different:
//
//     converged
//         = whether the CURRENT volumetric state is converged.
//
//     regular_phase_closed
//         = whether convergence has EVER been reached and therefore
//           the REGULAR-selection phase has permanently ended.
//
// A later CONCAVE cut may produce:
//
//     converged: 1 -> 0
//
// but:
//
//     regular_phase_closed remains true.
//
// Hence previously skipped REGULAR loops do NOT reopen.
//
static bool is_rl_action_legal(
    const GlobalState &state,
    const uint loop_id,
    const bool regular_phase_closed)
{
    if (!is_base_action_valid(
            state,
            loop_id))
    {
        return false;
    }

    const Loop &loop =
        state.loops.at(loop_id);

    //
    // After the first convergence event,
    // only remaining CONCAVE loops are legal.
    //
    if (regular_phase_closed &&
        loop.type != CONCAVE)
    {
        return false;
    }

    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static uint count_available_actions(
    const GlobalState &state,
    const bool regular_phase_closed)
{
    uint count = 0;

    for (uint i = 0;
         i < state.loops.size();
         ++i)
    {
        if (is_rl_action_legal(
                state,
                i,
                regular_phase_closed))
        {
            ++count;
        }
    }

    return count;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// An RL loop-selection episode terminates when there are no legal
// Stage-2 loop actions left.
//
// Note:
//
//     terminal != selection_success
//
// It is possible to have:
//
//     terminal = true
//     converged = false
//
// This happens, for example, when:
//
//     1. convergence was reached once,
//     2. REGULAR phase was permanently closed,
//     3. a later mandatory CONCAVE cut destroyed convergence,
//     4. all remaining CONCAVE loops were processed,
//     5. convergence was not recovered.
//
// Such an episode is a valid terminal loop-selection state, but is a
// selection failure. The original LoopyCuts finalization stage may still
// subsequently attempt CUTS TO DO / CUTS TO UNDO recovery.
//
static bool rl_terminal(
    const GlobalState &state,
    const bool regular_phase_closed)
{
    return count_available_actions(
               state,
               regular_phase_closed) == 0;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static void print_rl_state(
    const GlobalState &state,
    const bool converged,
    const bool regular_phase_closed,
    const uint step_count,
    const bool finalized,
    const bool diagnostics_valid,
    const uint nonmanifold_polys,
    const uint high_genus_polys,
    const uint buggy_chains)
{
    const uint available =
        finalized
            ? 0
            : count_available_actions(
                  state,
                  regular_phase_closed);

    const bool terminal =
        finalized
            ? true
            : (available == 0);

    const bool selection_success =
        terminal && converged;

    std::cout
        << "[RL] STATE"
        << " step=" << step_count
        << " loops=" << state.loops.size()
        << " available=" << available
        << " verts=" << state.m_vol.num_verts()
        << " tets=" << state.m_vol.num_polys()
        << " mm_verts=" << state.m_poly.num_verts()
        << " mm_edges=" << state.m_poly.num_edges()
        << " mm_faces=" << state.m_poly.num_faces()
        << " mm_polys=" << state.m_poly.num_polys()
        << " converged=" << converged
        << " regular_phase_closed="
        << regular_phase_closed
        << " terminal=" << terminal
        << " selection_success="
        << selection_success
        << " finalized=" << finalized
        << " diagnostics_valid="
        << diagnostics_valid
        << " nonmanifold_polys="
        << nonmanifold_polys
        << " high_genus_polys="
        << high_genus_polys
        << " buggy_chains="
        << buggy_chains
        << std::endl;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Dynamic per-loop state.
//
// These lines contain ORIGINAL LoopyCuts loop IDs.
// Python will later convert them into fixed-size boolean arrays.
//
// IMPORTANT:
//
//     EXISTS != ACTION LEGALITY
//
// A loop may exist but be unavailable because it is:
//
//     used,
//     reverted,
//     Nico_bug,
//     TOP_RELEVANT,
//     or closed by the RL phase semantics.
//
// [RL] ACTIONS remains the authoritative legality mask.
//
static void print_loop_status(
    const GlobalState &state)
{
    std::cout
        << "[RL] USED";

    for (uint i = 0;
         i < state.loops.size();
         ++i)
    {
        if (state.loops.at(i).used)
        {
            std::cout
                << " "
                << i;
        }
    }

    std::cout
        << std::endl;

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    std::cout
        << "[RL] REVERTED";

    for (uint i = 0;
         i < state.loops.size();
         ++i)
    {
        if (state.loops.at(i).reverted)
        {
            std::cout
                << " "
                << i;
        }
    }

    std::cout
        << std::endl;

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    std::cout
        << "[RL] NICO_BUG";

    for (uint i = 0;
         i < state.loops.size();
         ++i)
    {
        if (state.loops.at(i).Nico_bug)
        {
            std::cout
                << " "
                << i;
        }
    }

    std::cout
        << std::endl;

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    std::cout
        << "[RL] TOP_RELEVANT";

    for (uint i = 0;
         i < state.loops.size();
         ++i)
    {
        if (state.loops.at(i).type ==
            TOP_RELEVANT)
        {
            std::cout
                << " "
                << i;
        }
    }

    std::cout
        << std::endl;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static void print_available_actions(
    const GlobalState &state,
    const bool regular_phase_closed,
    const bool finalized)
{
    std::cout
        << "[RL] ACTIONS";

    if (!finalized)
    {
        for (uint i = 0;
             i < state.loops.size();
             ++i)
        {
            if (is_rl_action_legal(
                    state,
                    i,
                    regular_phase_closed))
            {
                std::cout
                    << " "
                    << i;
            }
        }
    }

    std::cout
        << std::endl;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Perform the same post-processing pipeline used by batch mode.
//
// IMPORTANT:
//
// This mutates state.m_poly / state.m_vol and must only be called once
// after the loop-selection episode has terminated.
//
static void finalize_rl_episode(
    GlobalState &state,
    const std::string &output_dir,
    const bool save_outputs)
{
    const std::string model =
        get_file_name(
            model_name,
            false);

    if (save_outputs)
    {
        std::cout
            << "[RL] FINALIZE_BEGIN"
            << " output="
            << output_dir
            << std::endl;
    }
    else
    {
        std::cout
            << "[RL] FINALIZE_EVAL_BEGIN"
            << std::endl;
    }

    //
    // IMPORTANT:
    //
    // FINALIZE and FINALIZE_EVAL intentionally execute the exact
    // same geometric post-processing pipeline.
    //
    // FINALIZE_EVAL differs ONLY by skipping file serialization.
    //
    finalize_block_decomposition(
        state);

    if (save_outputs)
    {
        state.m_poly.save(
            (output_dir + "/" +
             model +
             "_mm.hedra")
                .c_str());
    }

    SubdivisionHelper sh(
        state.m_vol,
        state.m_poly);

    state.m_poly =
        sh.subdivide();

    if (save_outputs)
    {
        state.m_poly.save(
            (output_dir + "/" +
             model +
             "_mm_subdivided.hedra")
                .c_str());
    }

    smoother(
        state.m_poly,
        state.m_srf);

    if (save_outputs)
    {
        state.m_poly.save(
            (output_dir + "/" +
             model +
             "_mm_subdivided_smoothed.hedra")
                .c_str());
    }

    classify_polyhedra(
        state.m_poly);

    //
    // MAKE HEXMESH
    //
    Hexmesh<MM, MV, ME, MF, MP> hm;

    export_hexahedra(
        state.m_poly,
        hm);

    const uint nh =
        hm.num_polys();

    const uint np =
        state.m_poly.num_polys();

    const bool full_hex =
        (nh == np);

    if (full_hex)
    {
        //
        // This MUST run in both FINALIZE and FINALIZE_EVAL.
        //
        // A failure here is part of the real finalization outcome
        // and must not be hidden merely because no file is saved.
        //
        hm.poly_fix_orientation();

        if (save_outputs)
        {
            hm.save(
                (output_dir + "/" +
                 model +
                 "_hex.mesh")
                    .c_str());
        }
    }

    std::cout
        << "[RL] FINAL_RESULT"
        << " hex=" << nh
        << " total_polys=" << np
        << " full_hex=" << full_hex
        << std::endl;

    if (save_outputs)
    {
        std::cout
            << "[RL] FINALIZE_END"
            << std::endl;
    }
    else
    {
        std::cout
            << "[RL] FINALIZE_EVAL_END"
            << std::endl;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int run_rl_server(
    GlobalState &state)
{
    //
    // Current geometric convergence state.
    //
    bool converged = false;

    //
    // RL-V1 phase flag.
    //
    // false:
    //     CONCAVE + REGULAR may be selected.
    //
    // true:
    //     REGULAR loops are permanently closed.
    //     Only remaining CONCAVE loops may be selected.
    //
    // This flag changes only:
    //
    //     false -> true
    //
    // and never returns to false.
    //
    bool regular_phase_closed = false;

    bool finalized = false;

    uint step_count = 0;

    //
    // MeshExtractor diagnostics are not defined as an RL observation
    // until execute_cut_step() has built the first post-step extractor.
    //
    bool diagnostics_valid = false;

    uint nonmanifold_polys = 0;
    uint high_genus_polys = 0;
    uint buggy_chains = 0;


    std::cout
        << std::endl;

    std::cout
        << "=========================================="
        << std::endl;

    std::cout
        << " LOOPYCUTS RL SERVER"
        << std::endl;

    std::cout
        << "=========================================="
        << std::endl;

    std::cout
        << "[RL] READY"
        << " loops="
        << state.loops.size()
        << " verts="
        << state.m_vol.num_verts()
        << " tets="
        << state.m_vol.num_polys()
        << std::endl;

    print_rl_state(
        state,
        converged,
        regular_phase_closed,
        step_count,
        finalized,
        diagnostics_valid,
        nonmanifold_polys,
        high_genus_polys,
        buggy_chains);

    print_loop_status(
        state);

    print_available_actions(
        state,
        regular_phase_closed,
        finalized);

    std::string line;

    while (std::getline(
        std::cin,
        line))
    {
        if (line.empty())
        {
            continue;
        }

        std::istringstream iss(
            line);

        std::string command;

        iss >> command;

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        if (command == "STATE")
        {
            print_rl_state(
                state,
                converged,
                regular_phase_closed,
                step_count,
                finalized,
                diagnostics_valid,
                nonmanifold_polys,
                high_genus_polys,
                buggy_chains);

            print_loop_status(
                state);

            print_available_actions(
                state,
                regular_phase_closed,
                finalized);

            continue;
        }

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        if (command == "STEP")
        {
            if (finalized)
            {
                std::cout
                    << "[RL] ERROR"
                    << " reason=ALREADY_FINALIZED"
                    << std::endl;

                continue;
            }

            long long input_id = -1;

            if (!(iss >> input_id))
            {
                std::cout
                    << "[RL] ERROR"
                    << " reason=MISSING_LOOP_ID"
                    << std::endl;

                continue;
            }

            if (input_id < 0 ||
                input_id >=
                    static_cast<long long>(
                        state.loops.size()))
            {
                std::cout
                    << "[RL] ERROR"
                    << " reason=INVALID_LOOP_ID"
                    << " loop_id="
                    << input_id
                    << std::endl;

                continue;
            }

            const uint loop_id =
                static_cast<uint>(
                    input_id);

            const Loop &loop =
                state.loops.at(
                    loop_id);

            //
            // Detailed base-validity errors.
            //
            if (loop.type == CONVEX)
            {
                std::cout
                    << "[RL] ERROR"
                    << " reason=CONVEX_LOOP"
                    << " loop_id="
                    << loop_id
                    << std::endl;

                continue;
            }

            if (loop.used)
            {
                std::cout
                    << "[RL] ERROR"
                    << " reason=ALREADY_USED"
                    << " loop_id="
                    << loop_id
                    << std::endl;

                continue;
            }

            if (loop.reverted)
            {
                std::cout
                    << "[RL] ERROR"
                    << " reason=ALREADY_REVERTED"
                    << " loop_id="
                    << loop_id
                    << std::endl;

                continue;
            }

            if (loop.Nico_bug)
            {
                std::cout
                    << "[RL] ERROR"
                    << " reason=NICO_BUG"
                    << " loop_id="
                    << loop_id
                    << std::endl;

                continue;
            }

            //
            // RL-V1 phase legality.
            //
            // Once the first convergence event has closed the
            // REGULAR-selection phase, REGULAR loops can never
            // become legal again, even if a later CONCAVE cut
            // changes convergence from 1 to 0.
            //
            if (regular_phase_closed &&
                loop.type != CONCAVE)
            {
                std::cout
                    << "[RL] ERROR"
                    << " reason=REGULAR_PHASE_CLOSED"
                    << " loop_id="
                    << loop_id
                    << std::endl;

                continue;
            }

            //
            // Record semantic state BEFORE the action.
            //
            const bool converged_before =
                converged;

            const bool
                regular_phase_closed_before =
                    regular_phase_closed;

            const char *loop_type_name =
                rl_loop_type_string(
                    loop);

            //
            // Execute one real LoopyCuts Stage-2 action.
            //
            CutStepResult result =
                execute_cut_step(
                    state,
                    loop_id);

            diagnostics_valid = true;

            nonmanifold_polys =
                result.nonmanifold_polys;

            high_genus_polys =
                result.high_genus_polys;

            buggy_chains =
                result.buggy_chains;

            //
            // Update CURRENT geometric convergence.
            //
            converged =
                result.converged;

            //
            // RL-V1:
            //
            // The first convergence event permanently closes
            // the REGULAR-selection phase.
            //
            // This is intentionally monotonic.
            //
            if (converged)
            {
                regular_phase_closed =
                    true;
            }

            ++step_count;

            std::cout
                << "[RL] STEP_RESULT"
                << " step="
                << step_count
                << " loop_id="
                << loop_id
                << " loop_type="
                << loop_type_name
                << " converged_before="
                << converged_before
                << " regular_phase_closed_before="
                << regular_phase_closed_before
                << " status="
                << cut_step_status_string(
                       result.status)
                << " attempted="
                << result.attempted
                << " committed="
                << result.committed
                << " reverted="
                << result.reverted
                << " verts_before="
                << result.verts_before
                << " verts_after="
                << result.verts_after
                << " tets_before="
                << result.tets_before
                << " tets_after="
                << result.tets_after
                << " nonmanifold_polys="
                << result.nonmanifold_polys
                << " high_genus_polys="
                << result.high_genus_polys
                << " buggy_chains="
                << result.buggy_chains
                << " step_time="
                << result.step_time
                << " converged="
                << converged
                << " regular_phase_closed="
                << regular_phase_closed
                << std::endl;

            print_rl_state(
                state,
                converged,
                regular_phase_closed,
                step_count,
                finalized,
                diagnostics_valid,
                nonmanifold_polys,
                high_genus_polys,
                buggy_chains);

            print_loop_status(
                state);

            print_available_actions(
                state,
                regular_phase_closed,
                finalized);

            continue;
        }

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        if (command == "FINALIZE" ||
            command == "FINALIZE_EVAL")
        {
            if (finalized)
            {
                std::cout
                    << "[RL] ERROR"
                    << " reason=ALREADY_FINALIZED"
                    << std::endl;

                continue;
            }

            const bool save_outputs =
                (command == "FINALIZE");

            std::string output_dir;

            if (save_outputs)
            {
                if (!(iss >> output_dir))
                {
                    std::cout
                        << "[RL] ERROR"
                        << " reason=MISSING_OUTPUT_DIR"
                        << std::endl;

                    continue;
                }
            }
            else
            {
                //
                // FINALIZE_EVAL is deliberately argument-free.
                //
                std::string unexpected_argument;

                if (iss >> unexpected_argument)
                {
                    std::cout
                        << "[RL] ERROR"
                        << " reason=UNEXPECTED_ARGUMENT"
                        << " command=FINALIZE_EVAL"
                        << std::endl;

                    continue;
                }
            }

            //
            // Do not finalize while there are still legal
            // Stage-2 loop actions.
            //
            // Note that terminal does NOT require
            // converged == true.
            //
            // A failed loop-selection trajectory may still
            // hand control to LoopyCuts finalization.
            //
            if (!rl_terminal(
                    state,
                    regular_phase_closed))
            {
                std::cout
                    << "[RL] ERROR"
                    << " reason=EPISODE_NOT_TERMINAL"
                    << " available="
                    << count_available_actions(
                           state,
                           regular_phase_closed)
                    << std::endl;

                continue;
            }

            finalize_rl_episode(
                state,
                output_dir,
                save_outputs);

            finalized = true;

            //
            // Finalization mutates the decomposition.
            // The diagnostics stored from the final RL selection step
            // no longer describe this new state.
            //
            diagnostics_valid = false;

            nonmanifold_polys = 0;
            high_genus_polys = 0;
            buggy_chains = 0;

            print_rl_state(
                state,
                converged,
                regular_phase_closed,
                step_count,
                finalized,
                diagnostics_valid,
                nonmanifold_polys,
                high_genus_polys,
                buggy_chains);

            print_loop_status(
                state);

            print_available_actions(
                state,
                regular_phase_closed,
                finalized);

            continue;
        }

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        if (command == "QUIT")
        {
            std::cout
                << "[RL] BYE"
                << std::endl;

            return 0;
        }

        //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        std::cout
            << "[RL] ERROR"
            << " reason=UNKNOWN_COMMAND"
            << " command="
            << command
            << std::endl;
    }

    return 0;
}