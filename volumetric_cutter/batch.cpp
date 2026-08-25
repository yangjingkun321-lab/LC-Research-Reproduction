/***************************************************************************/
/* Copyright(C) 2020

 Marco Livesu
 Italian National Research Council

 and

 Nico Pietroni
 University Of Technology Sydney

 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

#include "batch.h"
#include "cut.h"
#include "cut_step.h"
#include "polyhedral_decomposition.h"
#include "subdivision_helper.h"
#include "mesh_extractor.h"
#include "finalization.h"
#include "mesh_smoother.h"
#include <cinolib/export_hexahedra.h>
#include <cinolib/string_utilities.h>

#include <fstream>
#include <vector>
#include <iostream>


extern std::string output_name;
extern std::string model_name;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void run_batch(GlobalState & state)
{
    state.profiler.push("LoopyCuts");

    uint cut_count = 0;
    bool converged = false;

    for(uint i=0; i<state.loops.size(); ++i)
    {
        if(state.loops.at(i).type==CONCAVE || !converged)
        {
            if(state.loops.at(i).type==CONVEX) continue;

            CutStepResult result =
                execute_cut_step(state, i);

            converged = result.converged;

            ++cut_count;
        }
    }

    auto model = get_file_name(model_name,false);

    std::cout << cut_count << " cuts performed" << std::endl;
    finalize_block_decomposition(state);
    state.m_poly.save((output_name + "/" + model + "_mm.hedra").c_str());
    SubdivisionHelper sh(state.m_vol, state.m_poly);
    state.m_poly = sh.subdivide();
    state.m_poly.save((output_name + "/" + model + "_mm_subdivided.hedra").c_str());
    smoother(state.m_poly, state.m_srf);
    state.m_poly.save((output_name + "/" + model + "_mm_subdivided_smoothed.hedra").c_str());
    classify_polyhedra(state.m_poly);
    state.profiler.pop();

    // MAKE HEXMESH
    Hexmesh<MM,MV,ME,MF,MP> hm;
    export_hexahedra(state.m_poly, hm);
    uint  nh = hm.num_polys();
    uint  np = state.m_poly.num_polys();
    if(nh==np)
    {
        hm.poly_fix_orientation();
        hm.save((output_name + "/" + model + "_hex.mesh").c_str());
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// External-order mode
//
// The order file must contain every loop ID exactly once.
// Example for 5 loops:
//
// 0
// 1
// 2
// 3
// 4
//
// IDs may also be separated by spaces.
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static bool load_external_order(const std::string &filename,
                                const uint num_loops,
                                std::vector<uint> &order)
{
    std::ifstream in(filename.c_str());

    if (!in.is_open())
    {
        std::cerr << "[EXTERNAL_ORDER_ERROR] Cannot open order file: "
                  << filename << std::endl;
        return false;
    }

    order.clear();

    std::vector<bool> seen(num_loops, false);

    long long id = -1;

    while (in >> id)
    {
        if (id < 0 || id >= static_cast<long long>(num_loops))
        {
            std::cerr
                << "[EXTERNAL_ORDER_ERROR] Invalid loop id "
                << id
                << ". Valid range is [0, "
                << (num_loops == 0 ? 0 : num_loops - 1)
                << "]."
                << std::endl;

            return false;
        }

        const uint uid = static_cast<uint>(id);

        if (seen.at(uid))
        {
            std::cerr
                << "[EXTERNAL_ORDER_ERROR] Loop id "
                << uid
                << " appears more than once."
                << std::endl;

            return false;
        }

        seen.at(uid) = true;
        order.push_back(uid);
    }

    // Detect malformed contents such as letters.
    if (!in.eof())
    {
        std::cerr
            << "[EXTERNAL_ORDER_ERROR] Order file must contain "
            << "integer loop IDs only."
            << std::endl;

        return false;
    }

    if (order.size() != num_loops)
    {
        std::cerr
            << "[EXTERNAL_ORDER_ERROR] Order file contains "
            << order.size()
            << " loop IDs, but the model contains "
            << num_loops
            << " loops."
            << std::endl;

        for (uint i = 0; i < num_loops; ++i)
        {
            if (!seen.at(i))
            {
                std::cerr
                    << "[EXTERNAL_ORDER_ERROR] Missing loop id "
                    << i
                    << std::endl;
            }
        }

        return false;
    }

    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool run_external_order(GlobalState &state,
                        const std::string &order_file)
{
    std::vector<uint> order;

    if (!load_external_order(order_file,
                             state.loops.size(),
                             order))
    {
        return false;
    }

    std::cout << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << " EXTERNAL LOOP ORDER MODE" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "order file : " << order_file << std::endl;
    std::cout << "loop count : " << order.size() << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;

    state.profiler.push("LoopyCuts");

    // Keep this counter for compatibility with the original experimental logs.
    // It counts calls entering the cut-processing branch, not necessarily
    // successfully committed geometric cuts.
    uint cut_count = 0;

    bool converged = false;

    // More precise statistics for External / future RL mode.
    uint selected_count = 0;
    uint attempted_count = 0;
    uint committed_count = 0;
    uint reverted_count = 0;

    uint noop_used_count = 0;
    uint noop_reverted_count = 0;
    uint noop_nico_bug_count = 0;

    uint skipped_convex_count = 0;
    uint skipped_converged_count = 0;

    uint max_verts = state.m_vol.num_verts();
    uint max_tets = state.m_vol.num_polys();

    for (uint pos = 0; pos < order.size(); ++pos)
    {
        const uint i = order.at(pos);

        std::cout
            << "[EXTERNAL_ORDER] position="
            << pos
            << " loop_id="
            << i
            << std::endl;

        //
        // Preserve exactly the same CONCAVE / convergence logic
        // as the original run_batch().
        //
        if (state.loops.at(i).type == CONCAVE || !converged)
        {
            //
            // CONVEX loops are feature loops and are not volumetric cuts.
            //
            if (state.loops.at(i).type == CONVEX)
            {
                ++skipped_convex_count;

                std::cout
                    << "[EXTERNAL_STEP] "
                    << "position=" << pos
                    << " loop_id=" << i
                    << " status=SKIP_CONVEX"
                    << std::endl;

                continue;
            }

            ++selected_count;

            //
            // Execute exactly one Stage-2 action using the shared core.
            //
            CutStepResult result =
                execute_cut_step(state, i);

            ++cut_count;

            //
            // Update External Order statistics.
            //
            if (result.attempted)
            {
                ++attempted_count;
            }

            if (result.committed)
            {
                ++committed_count;
            }

            if (result.reverted)
            {
                ++reverted_count;
            }

            if (result.status == CutStepStatus::NOOP_USED)
            {
                ++noop_used_count;
            }
            else if (result.status == CutStepStatus::NOOP_REVERTED)
            {
                ++noop_reverted_count;
            }
            else if (result.status == CutStepStatus::NOOP_NICO_BUG)
            {
                ++noop_nico_bug_count;
            }

            converged = result.converged;

            if (result.verts_after > max_verts)
            {
                max_verts = result.verts_after;
            }

            if (result.tets_after > max_tets)
            {
                max_tets = result.tets_after;
            }

            //
            // Machine-readable per-step statistics.
            //
            std::cout
                << "[EXTERNAL_STEP] "
                << "position=" << pos
                << " loop_id=" << i
                << " status="
                << cut_step_status_string(result.status)
                << " verts_before=" << result.verts_before
                << " verts_after=" << result.verts_after
                << " tets_before=" << result.tets_before
                << " tets_after=" << result.tets_after
                << " step_time=" << result.step_time
                << " converged=" << result.converged
                << std::endl;
        }
        else
        {
            ++skipped_converged_count;

            std::cout
                << "[EXTERNAL_STEP] "
                << "position=" << pos
                << " loop_id=" << i
                << " status=SKIP_CONVERGED"
                << std::endl;
        }
    }

    //
    // Everything below is intentionally identical to run_batch().
    // We duplicate it for this first experiment so that the original
    // batch implementation remains untouched.
    //

    auto model = get_file_name(model_name, false);

    std::cout << std::endl;

    std::cout
        << "[EXTERNAL_SUMMARY] "
        << "selected=" << selected_count
        << " attempted=" << attempted_count
        << " committed=" << committed_count
        << " reverted=" << reverted_count
        << " noop_used=" << noop_used_count
        << " noop_reverted=" << noop_reverted_count
        << " noop_nico_bug=" << noop_nico_bug_count
        << " skipped_convex=" << skipped_convex_count
        << " skipped_converged=" << skipped_converged_count
        << " max_verts=" << max_verts
        << " max_tets=" << max_tets
        << " converged=" << converged
        << std::endl;

    std::cout << cut_count << " cuts performed" << std::endl;

    finalize_block_decomposition(state);

    state.m_poly.save(
        (output_name + "/" + model + "_mm.hedra").c_str());

    SubdivisionHelper sh(state.m_vol, state.m_poly);

    state.m_poly = sh.subdivide();

    state.m_poly.save(
        (output_name + "/" + model + "_mm_subdivided.hedra").c_str());

    smoother(state.m_poly, state.m_srf);

    state.m_poly.save(
        (output_name + "/" + model + "_mm_subdivided_smoothed.hedra").c_str());

    classify_polyhedra(state.m_poly);

    state.profiler.pop();

    // MAKE HEXMESH
    Hexmesh<MM, MV, ME, MF, MP> hm;

    export_hexahedra(state.m_poly, hm);

    uint nh = hm.num_polys();
    uint np = state.m_poly.num_polys();

    if (nh == np)
    {
        hm.poly_fix_orientation();

        hm.save(
            (output_name + "/" + model + "_hex.mesh").c_str());
    }

    return true;
}