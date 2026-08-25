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

#include "subdivision_helper.h"
#include <cmath>
#include <limits>
#include <cinolib/geometry/n_sided_poygon.h>
#include <cinolib/sampling.h>
#include <cinolib/harmonic_map.h>
#include <cinolib/subdivision_midpoint.h>
#include <algorithm>

SubdivisionHelper::SubdivisionHelper(const TetMesh & m, const MetaMesh & mm)
: m(m)
, mm(mm)
{
    extract_surface_patches();
    make_uv_maps();
    subdivide();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void SubdivisionHelper::extract_surface_patches()
{
    /* meta mesh extrctor labels all facets in the tetmesh and associates
     * a unique label to each meta mesh face, creating a connection between
     * the two face sets. Here I exploit this connection to retrieve the patches
     * to be uv-maped
    */

    for(uint mm_fid=0; mm_fid<mm.num_faces(); ++mm_fid)
    {
        if(!mm.face_is_on_srf(mm_fid)) continue;
        int l = mm.face_data(mm_fid).label;
        if(l<0) continue; // from the second run on, no uv maps available (will do standard midpoint)

        std::cout << "extracting patch for MM surface face " << mm_fid << ". Face label on M is " << l << std::endl;

        Patch p;
        std::unordered_map<uint,uint> vmap;
        for(uint m_fid=0; m_fid<m.num_faces(); ++m_fid)
        {
            if(m.face_data(m_fid).label==l)
            {
                uint tri[3];
                for(int i=0; i<3; ++i)
                {
                    uint vid = m.face_vert_id(m_fid,i);
                    auto query = vmap.find(vid);
                    if(query==vmap.end())
                    {
                        uint fresh_id = p.m.vert_add(m.vert(vid));
                        vmap[vid] = fresh_id;
                        tri[i] = fresh_id;
                    }
                    else tri[i] = query->second;
                }
                p.m.poly_add(tri[0], tri[1], tri[2]);
            }
        }

        for(uint mm_corner : mm.adj_f2v(mm_fid))
        {
            uint m_corner = mm.vert_data(mm_corner).vid_on_m;
            uint p_corner = vmap.at(m_corner);
            p.corners.push_back(p_corner);
            p.mm2corners[mm_corner] = p_corner;
        }
        patches[l] = p;
        //patches[l].m.save(("/Users/cino/Desktop/" + std::to_string(l) + ".obj").c_str());
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void SubdivisionHelper::make_uv_maps()
{
    for(auto & p : patches)
    {
        std::cout << "uv mapping patch " << p.first << "(corners " << p.second.corners.size() << ")" << std::endl;
        if(map_to_polygon(p.second.m, p.second.corners)) p.second.octree.build_from_mesh_polys(p.second.m);
        else
        {
            p.second.flawed = true;

            const std::string dump_path =
                "/tmp/loopycuts_uv_patch_" +
                std::to_string(p.first) +
                ".obj";

            std::cout
                << "[UV_PATCH_FALLBACK]"
                << " patch=" << p.first
                << " mode=standard_midpoint"
                << " dump=" << dump_path
                << std::endl;

            p.second.m.save(dump_path.c_str());
        }
        //for(auto c : p.second.corners_uv) std::cout << "\t" << c.first << " => " << c.second << std::endl;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

bool SubdivisionHelper::map_to_polygon(Trimesh<> & m, const std::vector<uint> & corners)
{
    // make sure it's a topological disk
    //if(m.Euler_characteristic()!=1)
    //{
    //    std::cout << "WARNING: patch is not a topological disk!" << std::endl;
    //    return false;
    //}

    // get the ordered list of boundary vertices, starting from a corner
    m.vert_unmark_all();
    for(uint vid : corners) m.vert_data(vid).flags[MARKED] = true;
    std::vector<uint> border = m.get_ordered_boundary_vertices();

    if(border.empty())
    {
        std::cout
            << "[UV_PATCH_INVALID_BOUNDARY]"
            << " reason=empty_or_non_orientable_boundary"
            << " verts=" << m.num_verts()
            << " edges=" << m.num_edges()
            << " faces=" << m.num_polys()
            << " corners=" << corners.size()
            << " euler=" << m.Euler_characteristic()
            << std::endl;

        return false;
    }

    if(corners.empty())
    {
        std::cout
            << "[UV_PATCH_INVALID_BOUNDARY]"
            << " reason=no_corners"
            << " verts=" << m.num_verts()
            << " edges=" << m.num_edges()
            << " faces=" << m.num_polys()
            << " border_vertices=" << border.size()
            << " euler=" << m.Euler_characteristic()
            << std::endl;

        return false;
    }

    // Every meta-mesh corner must lie on the ordered patch boundary.
    for(uint corner : corners)
    {
        const auto it = std::find(border.begin(), border.end(), corner);

        if(it == border.end())
        {
            std::cout
                << "[UV_PATCH_INVALID_BOUNDARY]"
                << " reason=corner_not_on_boundary"
                << " missing_corner=" << corner
                << " verts=" << m.num_verts()
                << " edges=" << m.num_edges()
                << " faces=" << m.num_polys()
                << " border_vertices=" << border.size()
                << " corners=" << corners.size()
                << " euler=" << m.Euler_characteristic()
                << std::endl;

            return false;
        }
    }

    CIRCULAR_SHIFT_VEC(border, corners.front());

    // split the boundary into n edges, with n = #corners
    std::vector<std::vector<uint>> edges;
    for(uint i=0; i<border.size(); ++i)
    {
        std::vector<uint> e = { border.at(i) };
        for(uint j=i+1; j<border.size() && !m.vert_data(border.at(j)).flags[MARKED]; ++j,++i)
        {
            e.push_back(border.at(j));
        }
        e.push_back(border.at((i+1)%border.size()));
        edges.push_back(e);
    }

    if(edges.size() != corners.size())
    {
        std::cout
            << "[UV_PATCH_INVALID_BOUNDARY]"
            << " reason=boundary_segment_count_mismatch"
            << " segments=" << edges.size()
            << " corners=" << corners.size()
            << " border_vertices=" << border.size()
            << " verts=" << m.num_verts()
            << " edges=" << m.num_edges()
            << " faces=" << m.num_polys()
            << " euler=" << m.Euler_characteristic()
            << std::endl;

        return false;
    }

    std::vector<vec3d> poly = n_sided_polygon(vec3d(0,0,0), corners.size(), 1.0);
    std::map<uint,vec3d> dirichlet_bcs;
    for(uint i=0; i<poly.size(); ++i)
    {
        std::vector<vec3d> e_bcs = sample_within_interval(poly.at(i), poly.at((i+1)%poly.size()), edges.at(i).size());
        for(uint j=0; j<e_bcs.size(); ++j) dirichlet_bcs[edges.at(i).at(j)] = e_bcs.at(j);
    }

    // map the interior vertices
    m.copy_xyz_to_uvw(UVW_param); // save XYZ coordinates in UVW
    m.vector_verts() = harmonic_map_3d(m, dirichlet_bcs, 1, COTANGENT); // map to disk

    return true;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

MetaMesh SubdivisionHelper::subdivide()
{
    // do standard midpoint subdivision
    std::unordered_map<uint,uint> edge_verts;
    std::unordered_map<uint,uint> face_verts;
    std::unordered_map<uint,uint> poly_verts;
    subdivision_midpoint(mm, mm_refined, edge_verts, face_verts, poly_verts);

    // restore sharp creases;
    for(auto e : edge_verts)
    {
        uint mm_eid = e.first;
        bool crease = mm.edge_data(mm_eid).flags[MARKED];
        if(crease)
        {
            uint vid   = e.second;
            uint mm_v0 = mm.edge_vert_id(mm_eid,0);
            uint mm_v1 = mm.edge_vert_id(mm_eid,1);
            int  e0    = mm_refined.edge_id(mm_v0,vid); assert(e0>=0);
            int  e1    = mm_refined.edge_id(mm_v1,vid); assert(e1>=0);
            mm_refined.edge_data(e0).flags[MARKED]= true;
            mm_refined.edge_data(e1).flags[MARKED]= true;
        }
    }

    // if per patch UV maps are available, reposition surface midpoints using the map
    if(!patches.empty())
    {
        // reposition face midpoints using the uv map of its associated MM face
        for(auto f : face_verts)
        {
            uint mm_fid = f.first;
            if(!mm.face_is_on_srf(mm_fid)) continue;
            uint vid      = f.second;
            int patch_id = mm.face_data(mm_fid).label;
            auto patch_it = patches.find(patch_id);

            if(patch_it == patches.end())
            {
                std::cerr
                    << "[SUBDIVISION_PATCH_FALLBACK]"
                    << " type=face"
                    << " mm_fid=" << mm_fid
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " reason=missing_patch"
                    << std::endl;

                // subdivision_midpoint() has already assigned this vertex.
                // Keep the standard midpoint when no UV patch is available.
                mm_refined.vert_data(vid).color = Color::RED();
                continue;
            }

            Patch & p = patch_it->second;
            if(p.flawed)
            {
                mm_refined.vert_data(vid).color = Color::RED();
                continue;
            }

            vec3d query(0,0,0);
            for(auto c : p.corners) query += p.m.vert(c);
            query /= static_cast<double>(p.corners.size());

            uint   pid;
            vec3d  pos;
            double dist;
            const bool SubdivisionQueryFinite =
                std::isfinite(query.x()) &&
                std::isfinite(query.y()) &&
                std::isfinite(query.z());

            if (!SubdivisionQueryFinite)
            {
                std::cerr
                    << "[SUBDIVISION_NONFINITE_QUERY]"
                    << " type=face"
                    << " mm_fid=" << mm_fid
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " patch_corners="
                    << p.corners.size()
                    << " query=" << query
                    << std::endl;

                assert(false);
            }

            p.octree.closest_point(query, pid, pos, dist);

            const bool SubdivisionClosestFinite =
                pid < p.m.num_polys() &&
                std::isfinite(pos.x()) &&
                std::isfinite(pos.y()) &&
                std::isfinite(pos.z()) &&
                std::isfinite(dist);

            if (!SubdivisionClosestFinite)
            {
                std::cerr
                    << "[SUBDIVISION_NONFINITE_CLOSEST_POINT]"
                    << " type=face"
                    << " mm_fid=" << mm_fid
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " patch_corners="
                    << p.corners.size()
                    << " patch_polys="
                    << p.m.num_polys()
                    << " pid=" << pid
                    << " query=" << query
                    << " pos=" << pos
                    << " dist=" << dist
                    << std::endl;

                assert(false);
            }
            //std::cout << "found face midpoint in uv: " << pos << "\t" << dist << std::endl;
            double bary[3];
            p.m.poly_bary_coords(pid, pos, bary);

            const bool SubdivisionBaryFinite =
                std::isfinite(bary[0]) &&
                std::isfinite(bary[1]) &&
                std::isfinite(bary[2]);

            if (!SubdivisionBaryFinite)
            {
                const uint OriginalPid = pid;
                const vec3d OriginalPos = pos;
                const double OriginalDist = dist;

                const vec3d OriginalA =
                    p.m.poly_vert(OriginalPid,0);

                const vec3d OriginalB =
                    p.m.poly_vert(OriginalPid,1);

                const vec3d OriginalC =
                    p.m.poly_vert(OriginalPid,2);

                const vec3d OriginalU =
                    OriginalB - OriginalA;

                const vec3d OriginalV =
                    OriginalC - OriginalA;

                const double OriginalD00 =
                    OriginalU.dot(OriginalU);

                const double OriginalD01 =
                    OriginalU.dot(OriginalV);

                const double OriginalD11 =
                    OriginalV.dot(OriginalV);

                const double OriginalDen =
                    OriginalD00 * OriginalD11 -
                    OriginalD01 * OriginalD01;

                uint BestPid =
                    std::numeric_limits<uint>::max();

                vec3d BestPos(0,0,0);

                double BestDist =
                    std::numeric_limits<double>::infinity();

                double BestDen = 0.0;

                uint DegenerateTriangleCount = 0;
                uint NonFiniteTriangleCount = 0;

                for (
                    uint CandidatePid=0;
                    CandidatePid<p.m.num_polys();
                    ++CandidatePid
                )
                {
                    const vec3d A =
                        p.m.poly_vert(CandidatePid,0);

                    const vec3d B =
                        p.m.poly_vert(CandidatePid,1);

                    const vec3d C =
                        p.m.poly_vert(CandidatePid,2);

                    const bool VerticesFinite =
                        std::isfinite(A.x()) &&
                        std::isfinite(A.y()) &&
                        std::isfinite(A.z()) &&
                        std::isfinite(B.x()) &&
                        std::isfinite(B.y()) &&
                        std::isfinite(B.z()) &&
                        std::isfinite(C.x()) &&
                        std::isfinite(C.y()) &&
                        std::isfinite(C.z());

                    if (!VerticesFinite)
                    {
                        ++NonFiniteTriangleCount;
                        continue;
                    }

                    const vec3d EdgeU = B - A;
                    const vec3d EdgeV = C - A;

                    const double D00 =
                        EdgeU.dot(EdgeU);

                    const double D01 =
                        EdgeU.dot(EdgeV);

                    const double D11 =
                        EdgeV.dot(EdgeV);

                    const double Den =
                        D00 * D11 -
                        D01 * D01;

                    double Scale = D00 * D11;

                    if (Scale < 1e-30)
                        Scale = 1e-30;

                    // Den / Scale is approximately sin(angle)^2.
                    // Reject exact and numerically near-degenerate triangles.
                    if (
                        !std::isfinite(Den) ||
                        std::fabs(Den) <=
                            1e-12 * Scale
                    )
                    {
                        ++DegenerateTriangleCount;
                        continue;
                    }

                    const vec3d CandidatePos =
                        cinolib::triangle_closest_point(
                            query,
                            A,
                            B,
                            C
                        );

                    const double CandidateDist =
                        (CandidatePos-query).length();

                    const bool CandidateFinite =
                        std::isfinite(CandidatePos.x()) &&
                        std::isfinite(CandidatePos.y()) &&
                        std::isfinite(CandidatePos.z()) &&
                        std::isfinite(CandidateDist);

                    if (!CandidateFinite)
                    {
                        ++NonFiniteTriangleCount;
                        continue;
                    }

                    if (CandidateDist < BestDist)
                    {
                        BestPid = CandidatePid;
                        BestPos = CandidatePos;
                        BestDist = CandidateDist;
                        BestDen = Den;
                    }
                }

                if (
                    BestPid ==
                    std::numeric_limits<uint>::max()
                )
                {
                    std::cerr
                        << "[SUBDIVISION_DEGENERATE_UV_RECOVERY_FAILED]"
                    << " type=face"
                    << " mm_fid=" << mm_fid
                    << " refined_vid=" << vid
                        << " patch_id=" << patch_id
                        << " patch_corners="
                        << p.corners.size()
                        << " original_pid="
                        << OriginalPid
                        << " original_pos="
                        << OriginalPos
                        << " original_dist="
                        << OriginalDist
                        << " original_den="
                        << OriginalDen
                        << " patch_polys="
                        << p.m.num_polys()
                        << " degenerate_triangles="
                        << DegenerateTriangleCount
                        << " nonfinite_triangles="
                        << NonFiniteTriangleCount
                        << std::endl;

                    assert(false);
                }

                pid = BestPid;
                pos = BestPos;
                dist = BestDist;

                p.m.poly_bary_coords(
                    pid,
                    pos,
                    bary
                );

                const bool RecoveredBaryFinite =
                    std::isfinite(bary[0]) &&
                    std::isfinite(bary[1]) &&
                    std::isfinite(bary[2]);

                std::cerr
                    << "[SUBDIVISION_DEGENERATE_UV_TRIANGLE_RECOVERED]"
                    << " type=face"
                    << " mm_fid=" << mm_fid
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " patch_corners="
                    << p.corners.size()
                    << " original_pid="
                    << OriginalPid
                    << " replacement_pid="
                    << pid
                    << " original_den="
                    << OriginalDen
                    << " replacement_den="
                    << BestDen
                    << " original_dist="
                    << OriginalDist
                    << " replacement_dist="
                    << BestDist
                    << " degenerate_triangles="
                    << DegenerateTriangleCount
                    << " bary=("
                    << bary[0] << ","
                    << bary[1] << ","
                    << bary[2] << ")"
                    << std::endl;

                if (!RecoveredBaryFinite)
                {
                    std::cerr
                        << "[SUBDIVISION_DEGENERATE_UV_RECOVERY_NONFINITE_BARY]"
                        << " patch_id=" << patch_id
                        << " replacement_pid="
                        << pid
                        << std::endl;

                    assert(false);
                }
            }
            double u = p.m.poly_sample_param_at(pid, bary, U_param);
            double v = p.m.poly_sample_param_at(pid, bary, V_param);
            double w = p.m.poly_sample_param_at(pid, bary, W_param);

            const bool SubdivisionParamFinite =
                std::isfinite(u) &&
                std::isfinite(v) &&
                std::isfinite(w);

            if (!SubdivisionParamFinite)
            {
                std::cerr
                    << "[SUBDIVISION_NONFINITE_PARAM]"
                    << " type=face"
                    << " mm_fid=" << mm_fid
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " patch_corners="
                    << p.corners.size()
                    << " pid=" << pid
                    << " query=" << query
                    << " pos=" << pos
                    << " dist=" << dist
                    << " bary=("
                    << bary[0] << ","
                    << bary[1] << ","
                    << bary[2] << ")"
                    << " uvw=("
                    << u << ","
                    << v << ","
                    << w << ")"
                    << std::endl;

                for (
                    uint pvid :
                    p.m.adj_p2v(pid)
                )
                {
                    std::cerr
                        << "  [SUBDIVISION_BAD_TRI_VERTEX]"
                        << " pvid=" << pvid
                        << " uv=" << p.m.vert(pvid)
                        << std::endl;
                }

                assert(false);
            }
            vec3d new_pos(u,v,w);
            mm_refined.vert(vid) = new_pos;
        }

        // reposition edge midpoints using one of the uv map of its incident MM faces
        for(auto e : edge_verts)
        {
            uint mm_eid = e.first;
            if(!mm.edge_is_on_srf(mm_eid)) continue;
            uint vid      = e.second;
            assert(mm.edge_adj_srf_faces(mm_eid).size()>0);
            uint mm_fid   = mm.edge_adj_srf_faces(mm_eid).front();
            int patch_id = mm.face_data(mm_fid).label;
            auto patch_it = patches.find(patch_id);

            if(patch_it == patches.end())
            {
                std::cerr
                    << "[SUBDIVISION_PATCH_FALLBACK]"
                    << " type=edge"
                    << " mm_eid=" << mm_eid
                    << " mm_fid=" << mm_fid
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " reason=missing_patch"
                    << std::endl;

                // Preserve the standard midpoint.
                mm_refined.vert_data(vid).color = Color::RED();
                continue;
            }

            Patch & p = patch_it->second;
            if(p.flawed)
            {
                mm_refined.vert_data(vid).color = Color::RED();
                continue;
            }

            uint mm_v0    = mm.edge_vert_id(mm_eid,0);
            uint mm_v1    = mm.edge_vert_id(mm_eid,1);

            auto c0_it = p.mm2corners.find(mm_v0);
            auto c1_it = p.mm2corners.find(mm_v1);

            if(c0_it == p.mm2corners.end() ||
               c1_it == p.mm2corners.end())
            {
                std::cerr
                    << "[SUBDIVISION_PATCH_FALLBACK]"
                    << " type=edge"
                    << " mm_eid=" << mm_eid
                    << " mm_fid=" << mm_fid
                    << " mm_v0=" << mm_v0
                    << " mm_v1=" << mm_v1
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " reason=missing_corner_mapping"
                    << std::endl;

                // Preserve the standard midpoint.
                mm_refined.vert_data(vid).color = Color::RED();
                continue;
            }

            uint c0 = c0_it->second;
            uint c1 = c1_it->second;
            vec3d query = (p.m.vert(c0) + p.m.vert(c1))*0.5;

            uint   pid;
            vec3d  pos;
            double dist;
            const bool SubdivisionQueryFinite =
                std::isfinite(query.x()) &&
                std::isfinite(query.y()) &&
                std::isfinite(query.z());

            if (!SubdivisionQueryFinite)
            {
                std::cerr
                    << "[SUBDIVISION_NONFINITE_QUERY]"
                    << " type=edge"
                    << " mm_eid=" << mm_eid
                    << " mm_fid=" << mm_fid
                    << " mm_v0=" << mm_v0
                    << " mm_v1=" << mm_v1
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " patch_corners="
                    << p.corners.size()
                    << " query=" << query
                    << std::endl;

                assert(false);
            }

            p.octree.closest_point(query, pid, pos, dist);

            const bool SubdivisionClosestFinite =
                pid < p.m.num_polys() &&
                std::isfinite(pos.x()) &&
                std::isfinite(pos.y()) &&
                std::isfinite(pos.z()) &&
                std::isfinite(dist);

            if (!SubdivisionClosestFinite)
            {
                std::cerr
                    << "[SUBDIVISION_NONFINITE_CLOSEST_POINT]"
                    << " type=edge"
                    << " mm_eid=" << mm_eid
                    << " mm_fid=" << mm_fid
                    << " mm_v0=" << mm_v0
                    << " mm_v1=" << mm_v1
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " patch_corners="
                    << p.corners.size()
                    << " patch_polys="
                    << p.m.num_polys()
                    << " pid=" << pid
                    << " query=" << query
                    << " pos=" << pos
                    << " dist=" << dist
                    << std::endl;

                assert(false);
            }
            //std::cout << "found edge midpoint in uv: " << pos << "\tpid: " << pid << "\tdist: " << dist << std::endl;
            double bary[3];
            p.m.poly_bary_coords(pid, pos, bary);

            const bool SubdivisionBaryFinite =
                std::isfinite(bary[0]) &&
                std::isfinite(bary[1]) &&
                std::isfinite(bary[2]);

            if (!SubdivisionBaryFinite)
            {
                const uint OriginalPid = pid;
                const vec3d OriginalPos = pos;
                const double OriginalDist = dist;

                const vec3d OriginalA =
                    p.m.poly_vert(OriginalPid,0);

                const vec3d OriginalB =
                    p.m.poly_vert(OriginalPid,1);

                const vec3d OriginalC =
                    p.m.poly_vert(OriginalPid,2);

                const vec3d OriginalU =
                    OriginalB - OriginalA;

                const vec3d OriginalV =
                    OriginalC - OriginalA;

                const double OriginalD00 =
                    OriginalU.dot(OriginalU);

                const double OriginalD01 =
                    OriginalU.dot(OriginalV);

                const double OriginalD11 =
                    OriginalV.dot(OriginalV);

                const double OriginalDen =
                    OriginalD00 * OriginalD11 -
                    OriginalD01 * OriginalD01;

                uint BestPid =
                    std::numeric_limits<uint>::max();

                vec3d BestPos(0,0,0);

                double BestDist =
                    std::numeric_limits<double>::infinity();

                double BestDen = 0.0;

                uint DegenerateTriangleCount = 0;
                uint NonFiniteTriangleCount = 0;

                for (
                    uint CandidatePid=0;
                    CandidatePid<p.m.num_polys();
                    ++CandidatePid
                )
                {
                    const vec3d A =
                        p.m.poly_vert(CandidatePid,0);

                    const vec3d B =
                        p.m.poly_vert(CandidatePid,1);

                    const vec3d C =
                        p.m.poly_vert(CandidatePid,2);

                    const bool VerticesFinite =
                        std::isfinite(A.x()) &&
                        std::isfinite(A.y()) &&
                        std::isfinite(A.z()) &&
                        std::isfinite(B.x()) &&
                        std::isfinite(B.y()) &&
                        std::isfinite(B.z()) &&
                        std::isfinite(C.x()) &&
                        std::isfinite(C.y()) &&
                        std::isfinite(C.z());

                    if (!VerticesFinite)
                    {
                        ++NonFiniteTriangleCount;
                        continue;
                    }

                    const vec3d EdgeU = B - A;
                    const vec3d EdgeV = C - A;

                    const double D00 =
                        EdgeU.dot(EdgeU);

                    const double D01 =
                        EdgeU.dot(EdgeV);

                    const double D11 =
                        EdgeV.dot(EdgeV);

                    const double Den =
                        D00 * D11 -
                        D01 * D01;

                    double Scale = D00 * D11;

                    if (Scale < 1e-30)
                        Scale = 1e-30;

                    // Den / Scale is approximately sin(angle)^2.
                    // Reject exact and numerically near-degenerate triangles.
                    if (
                        !std::isfinite(Den) ||
                        std::fabs(Den) <=
                            1e-12 * Scale
                    )
                    {
                        ++DegenerateTriangleCount;
                        continue;
                    }

                    const vec3d CandidatePos =
                        cinolib::triangle_closest_point(
                            query,
                            A,
                            B,
                            C
                        );

                    const double CandidateDist =
                        (CandidatePos-query).length();

                    const bool CandidateFinite =
                        std::isfinite(CandidatePos.x()) &&
                        std::isfinite(CandidatePos.y()) &&
                        std::isfinite(CandidatePos.z()) &&
                        std::isfinite(CandidateDist);

                    if (!CandidateFinite)
                    {
                        ++NonFiniteTriangleCount;
                        continue;
                    }

                    if (CandidateDist < BestDist)
                    {
                        BestPid = CandidatePid;
                        BestPos = CandidatePos;
                        BestDist = CandidateDist;
                        BestDen = Den;
                    }
                }

                if (
                    BestPid ==
                    std::numeric_limits<uint>::max()
                )
                {
                    std::cerr
                        << "[SUBDIVISION_DEGENERATE_UV_RECOVERY_FAILED]"
                    << " type=edge"
                    << " mm_eid=" << mm_eid
                    << " mm_fid=" << mm_fid
                    << " mm_v0=" << mm_v0
                    << " mm_v1=" << mm_v1
                    << " refined_vid=" << vid
                        << " patch_id=" << patch_id
                        << " patch_corners="
                        << p.corners.size()
                        << " original_pid="
                        << OriginalPid
                        << " original_pos="
                        << OriginalPos
                        << " original_dist="
                        << OriginalDist
                        << " original_den="
                        << OriginalDen
                        << " patch_polys="
                        << p.m.num_polys()
                        << " degenerate_triangles="
                        << DegenerateTriangleCount
                        << " nonfinite_triangles="
                        << NonFiniteTriangleCount
                        << std::endl;

                    assert(false);
                }

                pid = BestPid;
                pos = BestPos;
                dist = BestDist;

                p.m.poly_bary_coords(
                    pid,
                    pos,
                    bary
                );

                const bool RecoveredBaryFinite =
                    std::isfinite(bary[0]) &&
                    std::isfinite(bary[1]) &&
                    std::isfinite(bary[2]);

                std::cerr
                    << "[SUBDIVISION_DEGENERATE_UV_TRIANGLE_RECOVERED]"
                    << " type=edge"
                    << " mm_eid=" << mm_eid
                    << " mm_fid=" << mm_fid
                    << " mm_v0=" << mm_v0
                    << " mm_v1=" << mm_v1
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " patch_corners="
                    << p.corners.size()
                    << " original_pid="
                    << OriginalPid
                    << " replacement_pid="
                    << pid
                    << " original_den="
                    << OriginalDen
                    << " replacement_den="
                    << BestDen
                    << " original_dist="
                    << OriginalDist
                    << " replacement_dist="
                    << BestDist
                    << " degenerate_triangles="
                    << DegenerateTriangleCount
                    << " bary=("
                    << bary[0] << ","
                    << bary[1] << ","
                    << bary[2] << ")"
                    << std::endl;

                if (!RecoveredBaryFinite)
                {
                    std::cerr
                        << "[SUBDIVISION_DEGENERATE_UV_RECOVERY_NONFINITE_BARY]"
                        << " patch_id=" << patch_id
                        << " replacement_pid="
                        << pid
                        << std::endl;

                    assert(false);
                }
            }
            double u = p.m.poly_sample_param_at(pid, bary, U_param);
            double v = p.m.poly_sample_param_at(pid, bary, V_param);
            double w = p.m.poly_sample_param_at(pid, bary, W_param);

            const bool SubdivisionParamFinite =
                std::isfinite(u) &&
                std::isfinite(v) &&
                std::isfinite(w);

            if (!SubdivisionParamFinite)
            {
                std::cerr
                    << "[SUBDIVISION_NONFINITE_PARAM]"
                    << " type=edge"
                    << " mm_eid=" << mm_eid
                    << " mm_fid=" << mm_fid
                    << " mm_v0=" << mm_v0
                    << " mm_v1=" << mm_v1
                    << " refined_vid=" << vid
                    << " patch_id=" << patch_id
                    << " patch_corners="
                    << p.corners.size()
                    << " pid=" << pid
                    << " query=" << query
                    << " pos=" << pos
                    << " dist=" << dist
                    << " bary=("
                    << bary[0] << ","
                    << bary[1] << ","
                    << bary[2] << ")"
                    << " uvw=("
                    << u << ","
                    << v << ","
                    << w << ")"
                    << std::endl;

                for (
                    uint pvid :
                    p.m.adj_p2v(pid)
                )
                {
                    std::cerr
                        << "  [SUBDIVISION_BAD_TRI_VERTEX]"
                        << " pvid=" << pvid
                        << " uv=" << p.m.vert(pvid)
                        << std::endl;
                }

                assert(false);
            }
            vec3d new_pos(u,v,w);
            mm_refined.vert(vid) = new_pos;
        }
    }

    return mm_refined;
}
