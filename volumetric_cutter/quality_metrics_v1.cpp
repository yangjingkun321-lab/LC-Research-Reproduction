#include "quality_metrics_v1.h"

#include <cinolib/octree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace loopycuts_quality_v1
{

using cinolib::Octree;


// ================================================================
// Common validation.
// ================================================================

static void require_finite(
    const double value,
    const char *name)
{
    if (!std::isfinite(value))
    {
        throw std::runtime_error(
            std::string(
                "quality metric non-finite: ")
            +
            name);
    }
}


// ================================================================
// D_C V1.
// ================================================================

DcMetricsV1 compute_dc_v1(
    const std::uint64_t n_total,
    const std::uint64_t n_hex)
{
    if (n_total == 0)
    {
        throw std::runtime_error(
            "D_C V1: N_total must be positive");
    }

    if (n_hex > n_total)
    {
        throw std::runtime_error(
            "D_C V1: N_hex exceeds N_total");
    }


    const std::uint64_t n_nonhex =
        n_total - n_hex;


    //
    // IMPORTANT:
    //
    // Arithmetic order is frozen.
    //
    const double h =
        static_cast<double>(
            n_hex)
        /
        static_cast<double>(
            n_total);


    const double log1p_n_nonhex =
        std::log1p(
            static_cast<double>(
                n_nonhex));


    const double log1p_n_total =
        std::log1p(
            static_cast<double>(
                n_total));


    const double log_ratio =
        log1p_n_nonhex
        /
        log1p_n_total;


    const double d_nonhex =
        1.0
        -
        log_ratio;


    const double h_times_d_nonhex =
        h
        *
        d_nonhex;


    const double d_c =
        std::sqrt(
            h_times_d_nonhex);


    require_finite(
        h,
        "D_C.H");

    require_finite(
        log1p_n_nonhex,
        "D_C.log1p_n_nonhex");

    require_finite(
        log1p_n_total,
        "D_C.log1p_n_total");

    require_finite(
        log_ratio,
        "D_C.log_ratio");

    require_finite(
        d_nonhex,
        "D_C.D_nonhex");

    require_finite(
        h_times_d_nonhex,
        "D_C.product");

    require_finite(
        d_c,
        "D_C");


    DcMetricsV1 out;

    out.n_total =
        n_total;

    out.n_hex =
        n_hex;

    out.n_nonhex =
        n_nonhex;

    out.h =
        h;

    out.log1p_n_nonhex =
        log1p_n_nonhex;

    out.log1p_n_total =
        log1p_n_total;

    out.log_ratio =
        log_ratio;

    out.d_nonhex =
        d_nonhex;

    out.h_times_d_nonhex =
        h_times_d_nonhex;

    out.d_c =
        d_c;

    return out;
}


// ================================================================
// Frozen geometry primitives.
// ================================================================

struct MetricTriangleV1
{
    unsigned int v0;
    unsigned int v1;
    unsigned int v2;

    double area;
    double local_h;
};


static double edge_length_v1(
    const vec3d &a,
    const vec3d &b)
{
    const double dx =
        a.x() - b.x();

    const double dy =
        a.y() - b.y();

    const double dz =
        a.z() - b.z();


    return std::sqrt(
        dx * dx
        +
        dy * dy
        +
        dz * dz);
}


static double triangle_area2_v1(
    const vec3d &p0,
    const vec3d &p1,
    const vec3d &p2)
{
    const double ux =
        p1.x() - p0.x();

    const double uy =
        p1.y() - p0.y();

    const double uz =
        p1.z() - p0.z();


    const double vx =
        p2.x() - p0.x();

    const double vy =
        p2.y() - p0.y();

    const double vz =
        p2.z() - p0.z();


    const double cx =
        uy * vz
        -
        uz * vy;

    const double cy =
        uz * vx
        -
        ux * vz;

    const double cz =
        ux * vy
        -
        uy * vx;


    //
    // Frozen final-triangle runtime arithmetic.
    //
    return std::sqrt(
        cx * cx
        +
        cy * cy
        +
        cz * cz);
}


static double triangle_local_scale_v1(
    const vec3d &p0,
    const vec3d &p1,
    const vec3d &p2)
{
    std::array<double, 3>
        lengths = {{
            edge_length_v1(
                p0,
                p1),

            edge_length_v1(
                p1,
                p2),

            edge_length_v1(
                p2,
                p0),
        }};


    std::sort(
        lengths.begin(),
        lengths.end());


    const double h =
        lengths[1];


    if (
        !std::isfinite(h)
        ||
        h <= 1e-18)
    {
        throw std::runtime_error(
            "invalid Stage2 local triangle scale");
    }


    return h;
}


// ================================================================
// Stage2 metric surface.
//
// Frozen Train49 contract:
//
//     one Stage2 OBJ face == one triangle
//     polygon ID == metric triangle ID
//
// We fail closed rather than silently changing triangle IDs.
// ================================================================

static std::vector<MetricTriangleV1>
build_stage2_triangles_v1(
    SrfMesh &mesh)
{
    if (
        mesh.num_verts() == 0
        ||
        mesh.num_polys() == 0)
    {
        throw std::runtime_error(
            "empty Stage2 metric surface");
    }


    std::vector<MetricTriangleV1>
        triangles;

    triangles.reserve(
        mesh.num_polys());


    for (
        unsigned int pid = 0;
        pid < mesh.num_polys();
        ++pid)
    {
        const std::vector<unsigned int>
            vids =
                mesh.poly_verts_id(
                    pid);


        if (vids.size() != 3)
        {
            throw std::runtime_error(
                "Stage2 metric face is not triangular");
        }


        const unsigned int a =
            vids[0];

        const unsigned int b =
            vids[1];

        const unsigned int c =
            vids[2];


        if (
            a == b
            ||
            b == c
            ||
            a == c)
        {
            throw std::runtime_error(
                "Stage2 metric face has repeated vertex ID");
        }


        const vec3d &p0 =
            mesh.vert(a);

        const vec3d &p1 =
            mesh.vert(b);

        const vec3d &p2 =
            mesh.vert(c);


        const double area2 =
            triangle_area2_v1(
                p0,
                p1,
                p2);


        if (
            !std::isfinite(
                area2)
            ||
            area2 <= 1e-18)
        {
            throw std::runtime_error(
                "invalid Stage2 metric triangle area");
        }


        MetricTriangleV1 t;

        t.v0 = a;
        t.v1 = b;
        t.v2 = c;

        t.area =
            0.5
            *
            area2;

        t.local_h =
            triangle_local_scale_v1(
                p0,
                p1,
                p2);


        triangles.push_back(
            t);
    }


    if (
        triangles.size()
        !=
        mesh.num_polys())
    {
        throw std::runtime_error(
            "Stage2 metric triangle ID mapping changed");
    }


    return triangles;
}


// ================================================================
// Final metric surface.
//
// Frozen final topology contract:
//
//     get_surface_faces() order
//     face_verts_id(fid, false)
//     fan:
//         (face[0], face[i], face[i+1])
//
// repeated IDs are skipped;
// area2 <= 1e-18 is skipped.
//
// Final area:
//
//     0.5 * sqrt(cx*cx + cy*cy + cz*cz)
//
// This source form is the one proven exact on Available47.
// ================================================================

static std::vector<MetricTriangleV1>
build_final_triangles_v1(
    MetaMesh &mesh)
{
    const std::vector<unsigned int>
        boundary_fids =
            mesh.get_surface_faces();


    if (boundary_fids.empty())
    {
        throw std::runtime_error(
            "final metric surface has no boundary faces");
    }


    std::vector<MetricTriangleV1>
        triangles;


    for (
        std::size_t bid = 0;
        bid < boundary_fids.size();
        ++bid)
    {
        const unsigned int fid =
            boundary_fids[
                bid
            ];


        const std::vector<unsigned int>
            face =
                mesh.face_verts_id(
                    fid,
                    false);


        if (face.size() < 3)
        {
            throw std::runtime_error(
                "final boundary face has fewer than 3 vertices");
        }


        for (
            std::size_t i = 1;
            i + 1 < face.size();
            ++i)
        {
            const unsigned int v0 =
                face[0];

            const unsigned int v1 =
                face[i];

            const unsigned int v2 =
                face[i + 1];


            if (
                v0 == v1
                ||
                v1 == v2
                ||
                v2 == v0)
            {
                continue;
            }


            const vec3d &p0 =
                mesh.vert(v0);

            const vec3d &p1 =
                mesh.vert(v1);

            const vec3d &p2 =
                mesh.vert(v2);


            const double area2 =
                triangle_area2_v1(
                    p0,
                    p1,
                    p2);


            if (!std::isfinite(area2))
            {
                throw std::runtime_error(
                    "non-finite final triangle area2");
            }


            if (area2 <= 1e-18)
            {
                continue;
            }


            MetricTriangleV1 t;

            t.v0 = v0;
            t.v1 = v1;
            t.v2 = v2;

            t.area =
                0.5
                *
                area2;

            t.local_h =
                0.0;


            if (
                !std::isfinite(
                    t.area)
                ||
                t.area <= 0.0)
            {
                throw std::runtime_error(
                    "invalid final sampling triangle area");
            }


            triangles.push_back(
                t);
        }
    }


    if (triangles.empty())
    {
        throw std::runtime_error(
            "no valid final metric triangles");
    }


    return triangles;
}


// ================================================================
// NumPy 1.24.4 np.sum(float64) structural mirror.
//
// Frozen:
//     PW_BLOCKSIZE = 128
//     ufunc outer buffer = 8192
// ================================================================

static const std::size_t
    NUMPY_PW_BLOCKSIZE_V1 =
        128;


static const std::size_t
    NUMPY_UFUNC_BUFSIZE_V1 =
        8192;


static double numpy_double_pairwise_sum_v1(
    const double *a,
    const std::size_t n)
{
    if (n < 8)
    {
        //
        // NumPy intentionally starts from -0.0.
        //
        double res =
            -0.0;


        for (
            std::size_t i = 0;
            i < n;
            ++i)
        {
            res +=
                a[i];
        }


        return res;
    }


    if (
        n
        <=
        NUMPY_PW_BLOCKSIZE_V1)
    {
        std::size_t i;


        double r[8];

        r[0] = a[0];
        r[1] = a[1];
        r[2] = a[2];
        r[3] = a[3];
        r[4] = a[4];
        r[5] = a[5];
        r[6] = a[6];
        r[7] = a[7];


        const std::size_t stop =
            n
            -
            (
                n
                %
                8
            );


        for (
            i = 8;
            i < stop;
            i += 8)
        {
            r[0] += a[i + 0];
            r[1] += a[i + 1];
            r[2] += a[i + 2];
            r[3] += a[i + 3];
            r[4] += a[i + 4];
            r[5] += a[i + 5];
            r[6] += a[i + 6];
            r[7] += a[i + 7];
        }


        double res =
            (
                (r[0] + r[1])
                +
                (r[2] + r[3])
            )
            +
            (
                (r[4] + r[5])
                +
                (r[6] + r[7])
            );


        for (
            ;
            i < n;
            ++i)
        {
            res +=
                a[i];
        }


        return res;
    }


    std::size_t n2 =
        n / 2;


    n2 -=
        n2
        %
        8;


    return (
        numpy_double_pairwise_sum_v1(
            a,
            n2)
        +
        numpy_double_pairwise_sum_v1(
            a + n2,
            n - n2)
    );
}


static double numpy_sum_v1(
    const std::vector<double> &values)
{
    if (values.empty())
    {
        throw std::runtime_error(
            "NumPy sum input is empty");
    }


    double result =
        0.0;


    for (
        std::size_t start = 0;
        start < values.size();
        start += NUMPY_UFUNC_BUFSIZE_V1)
    {
        const std::size_t remaining =
            values.size()
            -
            start;


        const std::size_t count =
            (
                remaining
                <
                NUMPY_UFUNC_BUFSIZE_V1
            )
            ?
            remaining
            :
            NUMPY_UFUNC_BUFSIZE_V1;


        result +=
            numpy_double_pairwise_sum_v1(
                values.data()
                +
                start,
                count);
    }


    require_finite(
        result,
        "numpy_sum_v1");


    return result;
}


// ================================================================
// Frozen final CDF.
// ================================================================

static std::vector<double>
build_final_cdf_v1(
    const std::vector<MetricTriangleV1> &triangles)
{
    std::vector<double>
        areas(
            triangles.size());


    for (
        std::size_t i = 0;
        i < triangles.size();
        ++i)
    {
        areas[i] =
            triangles[i].area;
    }


    const double total =
        numpy_sum_v1(
            areas);


    if (
        !std::isfinite(total)
        ||
        total <= 0.0)
    {
        throw std::runtime_error(
            "invalid final total area");
    }


    std::vector<double>
        probability(
            areas.size());


    for (
        std::size_t i = 0;
        i < areas.size();
        ++i)
    {
        probability[i] =
            areas[i]
            /
            total;
    }


    //
    // np.cumsum(p, dtype=np.float64)
    //
    std::vector<double>
        cdf(
            probability.size());


    double cumulative =
        0.0;


    for (
        std::size_t i = 0;
        i < probability.size();
        ++i)
    {
        cumulative +=
            probability[i];

        cdf[i] =
            cumulative;
    }


    const double
        pre_normalization_last =
            cdf.back();


    if (
        !std::isfinite(
            pre_normalization_last)
        ||
        pre_normalization_last <= 0.0)
    {
        throw std::runtime_error(
            "invalid final CDF normalization denominator");
    }


    //
    // cdf /= cdf[-1]
    //
    for (
        std::size_t i = 0;
        i < cdf.size();
        ++i)
    {
        cdf[i] /=
            pre_normalization_last;
    }


    if (cdf.back() != 1.0)
    {
        throw std::runtime_error(
            "final CDF does not terminate at exact 1");
    }


    for (
        std::size_t i = 1;
        i < cdf.size();
        ++i)
    {
        if (
            cdf[i]
            <
            cdf[i - 1])
        {
            throw std::runtime_error(
                "final CDF is not monotone");
        }
    }


    return cdf;
}


// ================================================================
// Frozen barycentric final sampling.
// ================================================================

static vec3d reconstruct_final_point_v1(
    MetaMesh &mesh,
    const MetricTriangleV1 &triangle,
    const FinalSurfaceDraw &draw)
{
    const double u =
        draw.u;

    const double v =
        draw.v;


    const double su =
        std::sqrt(
            u);


    const double w0 =
        1.0
        -
        su;


    const double one_minus_v =
        1.0
        -
        v;


    const double w1 =
        su
        *
        one_minus_v;


    const double w2 =
        su
        *
        v;


    const vec3d &p0 =
        mesh.vert(
            triangle.v0);

    const vec3d &p1 =
        mesh.vert(
            triangle.v1);

    const vec3d &p2 =
        mesh.vert(
            triangle.v2);


    //
    // Frozen:
    //
    //     term0 + term1 + term2
    //
    // explicitly left-associated.
    //
    const double t0x =
        w0 * p0.x();

    const double t0y =
        w0 * p0.y();

    const double t0z =
        w0 * p0.z();


    const double t1x =
        w1 * p1.x();

    const double t1y =
        w1 * p1.y();

    const double t1z =
        w1 * p1.z();


    const double t2x =
        w2 * p2.x();

    const double t2y =
        w2 * p2.y();

    const double t2z =
        w2 * p2.z();


    const double x =
        (t0x + t1x)
        +
        t2x;

    const double y =
        (t0y + t1y)
        +
        t2y;

    const double z =
        (t0z + t1z)
        +
        t2z;


    require_finite(
        x,
        "final_sample.x");

    require_finite(
        y,
        "final_sample.y");

    require_finite(
        z,
        "final_sample.z");


    return vec3d(
        x,
        y,
        z);
}


// ================================================================
// NumPy 1.24.4 percentile(method="linear") scalar mirror.
//
// For finite float64 data:
//
// virtual_index = (n - 1) * q / 100
//
// _lerp:
//
//     diff = b - a
//     result = a + diff*t
//
//     if t >= 0.5:
//         result = b - diff*(1-t)
//
// Sorting the complete finite array is value-equivalent to NumPy's
// partition for the selected order statistics.
// ================================================================

static double percentile_linear_v1(
    std::vector<double> values,
    const double percentile)
{
    if (values.empty())
    {
        throw std::runtime_error(
            "percentile input is empty");
    }


    if (
        !std::isfinite(
            percentile)
        ||
        percentile < 0.0
        ||
        percentile > 100.0)
    {
        throw std::runtime_error(
            "invalid percentile");
    }


    for (
        std::size_t i = 0;
        i < values.size();
        ++i)
    {
        if (
            !std::isfinite(
                values[i]))
        {
            throw std::runtime_error(
                "non-finite percentile input");
        }
    }


    std::sort(
        values.begin(),
        values.end());


    if (values.size() == 1)
    {
        return values[0];
    }


    const double virtual_index =
        (
            static_cast<double>(
                values.size() - 1)
            *
            percentile
        )
        /
        100.0;


    const double lower_double =
        std::floor(
            virtual_index);


    const std::size_t lower =
        static_cast<std::size_t>(
            lower_double);


    if (
        lower
        >=
        values.size())
    {
        throw std::runtime_error(
            "percentile lower index out of range");
    }


    if (
        lower
        ==
        values.size() - 1)
    {
        return values[
            lower
        ];
    }


    const std::size_t upper =
        lower + 1;


    const double gamma =
        virtual_index
        -
        lower_double;


    const double a =
        values[
            lower
        ];

    const double b =
        values[
            upper
        ];


    const double diff_b_a =
        b
        -
        a;


    //
    // NumPy 1.24.4 _lerp first expression.
    //
    double result =
        a
        +
        diff_b_a
        *
        gamma;


    //
    // NumPy _lerp precision branch.
    //
    if (gamma >= 0.5)
    {
        result =
            b
            -
            diff_b_a
            *
            (
                1.0
                -
                gamma
            );
    }


    require_finite(
        result,
        "percentile result");


    return result;
}


// ================================================================
// Q_spurious.
// ================================================================

SpuriousMetricsV1 compute_q_spurious_v1(
    SrfMesh &input_surface,
    MetaMesh &final_mesh,
    const QualityRef &ref)
{
    validate_quality_ref_v1(
        ref);


    if (
        ref.final_draws.size()
        !=
        FINAL_DRAW_COUNT_V1)
    {
        throw std::runtime_error(
            "unexpected FINAL_DRAWS count");
    }


    const std::vector<MetricTriangleV1>
        stage2_triangles =
            build_stage2_triangles_v1(
                input_surface);


    const std::vector<MetricTriangleV1>
        final_triangles =
            build_final_triangles_v1(
                final_mesh);


    const std::vector<double>
        cdf =
            build_final_cdf_v1(
                final_triangles);


    //
    // Frozen Stage2 nearest backend:
    //
    // Octree(10, 32)
    //
    Octree tree(
        10,
        32);


    for (
        std::size_t tid = 0;
        tid < stage2_triangles.size();
        ++tid)
    {
        const MetricTriangleV1 &t =
            stage2_triangles[
                tid
            ];


        tree.add_triangle(
            static_cast<unsigned int>(
                tid),
            {
                input_surface.vert(
                    t.v0),

                input_surface.vert(
                    t.v1),

                input_surface.vert(
                    t.v2),
            });
    }


    tree.build();


    std::vector<double>
        normalized;

    normalized.reserve(
        ref.final_draws.size());


    std::size_t coverage_1h_count =
        0;


    for (
        std::size_t i = 0;
        i < ref.final_draws.size();
        ++i)
    {
        const FinalSurfaceDraw &draw =
            ref.final_draws[
                i
            ];


        //
        // np.searchsorted(
        //     cdf,
        //     area_draw,
        //     side="right"
        // )
        //
        const std::vector<double>::const_iterator
            selected_it =
                std::upper_bound(
                    cdf.begin(),
                    cdf.end(),
                    draw.area_draw);


        const std::size_t final_tid =
            static_cast<std::size_t>(
                selected_it
                -
                cdf.begin());


        if (
            final_tid
            >=
            final_triangles.size())
        {
            throw std::runtime_error(
                "final draw triangle selection out of range");
        }


        const vec3d query =
            reconstruct_final_point_v1(
                final_mesh,
                final_triangles[
                    final_tid
                ],
                draw);


        unsigned int nearest_tid =
            0;


        vec3d closest(
            0.0,
            0.0,
            0.0);


        double distance_squared =
            0.0;


        tree.closest_point(
            query,
            nearest_tid,
            closest,
            distance_squared);


        if (
            nearest_tid
            >=
            stage2_triangles.size())
        {
            throw std::runtime_error(
                "Stage2 Octree returned invalid triangle ID");
        }


        if (
            !std::isfinite(
                distance_squared)
            ||
            distance_squared
            <
            -1e-18)
        {
            throw std::runtime_error(
                "Stage2 Octree returned invalid squared distance");
        }


        const double distance =
            std::sqrt(
                std::max(
                    0.0,
                    distance_squared));


        const double h =
            stage2_triangles[
                nearest_tid
            ].local_h;


        const double normalized_distance =
            distance
            /
            h;


        if (
            !std::isfinite(
                normalized_distance)
            ||
            normalized_distance < 0.0)
        {
            throw std::runtime_error(
                "invalid normalized Final->Input distance");
        }


        normalized.push_back(
            normalized_distance);


        //
        // Frozen:
        //
        // np.mean(normalized <= 1.0)
        //
        if (
            normalized_distance
            <=
            1.0)
        {
            ++coverage_1h_count;
        }
    }


    if (normalized.empty())
    {
        throw std::runtime_error(
            "no Final->Input samples");
    }


    const double coverage_1h =
        static_cast<double>(
            coverage_1h_count)
        /
        static_cast<double>(
            normalized.size());


    const double p99_h =
        percentile_linear_v1(
            normalized,
            99.0);


    const double tail_transform =
        1.0
        /
        std::max(
            1.0,
            p99_h);


    const double q_spurious =
        std::sqrt(
            coverage_1h
            *
            tail_transform);


    require_finite(
        coverage_1h,
        "Q_spurious.coverage_1h");

    require_finite(
        p99_h,
        "Q_spurious.p99_h");

    require_finite(
        tail_transform,
        "Q_spurious.tail");

    require_finite(
        q_spurious,
        "Q_spurious");


    SpuriousMetricsV1 out;

    out.stage2_triangle_count =
        stage2_triangles.size();

    out.final_triangle_count =
        final_triangles.size();

    out.sample_count =
        normalized.size();

    out.coverage_1h_count =
        coverage_1h_count;

    out.coverage_1h =
        coverage_1h;

    out.p99_h =
        p99_h;

    out.tail_transform =
        tail_transform;

    out.q_spurious =
        q_spurious;


    return out;
}


// ================================================================
// Q_missing.
//
// Frozen Input -> Final contract:
//
//     query XYZ = QualityRef::input_geometry[i].xyz
//     local h   = QualityRef::input_geometry[i].h
//
//     nearest surface = frozen final boundary fan triangulation
//
//     normalized = distance / input_h
//
//     coverage_1h = mean(normalized <= 1.0)
//
//     Q_missing = sqrt(
//         coverage_1h
//         *
//         (1 / max(1, p99_h))
//     )
// ================================================================

MissingMetricsV1 compute_q_missing_v1(
    MetaMesh &final_mesh,
    const QualityRef &ref)
{
    validate_quality_ref_v1(
        ref);


    if (
        ref.input_geometry.size()
        !=
        GEOMETRY_SAMPLE_COUNT_V1)
    {
        throw std::runtime_error(
            "unexpected INPUT_GEOMETRY count");
    }


    const std::vector<MetricTriangleV1>
        final_triangles =
            build_final_triangles_v1(
                final_mesh);


    //
    // Frozen final-surface nearest backend:
    //
    //     Octree(10, 32)
    //
    Octree tree(
        10,
        32);


    for (
        std::size_t tid = 0;
        tid < final_triangles.size();
        ++tid)
    {
        const MetricTriangleV1 &t =
            final_triangles[
                tid
            ];


        tree.add_triangle(
            static_cast<unsigned int>(
                tid),
            {
                final_mesh.vert(
                    t.v0),

                final_mesh.vert(
                    t.v1),

                final_mesh.vert(
                    t.v2),
            });
    }


    tree.build();


    std::vector<double>
        normalized;

    normalized.reserve(
        ref.input_geometry.size());


    std::size_t coverage_1h_count =
        0;


    for (
        std::size_t i = 0;
        i < ref.input_geometry.size();
        ++i)
    {
        const InputGeometrySample &sample =
            ref.input_geometry[
                i
            ];


        if (
            !std::isfinite(
                sample.x)
            ||
            !std::isfinite(
                sample.y)
            ||
            !std::isfinite(
                sample.z)
            ||
            !std::isfinite(
                sample.h)
            ||
            sample.h <= 1e-18)
        {
            throw std::runtime_error(
                "invalid INPUT_GEOMETRY sample");
        }


        const vec3d query(
            sample.x,
            sample.y,
            sample.z);


        unsigned int nearest_tid =
            0;


        vec3d closest(
            0.0,
            0.0,
            0.0);


        double distance_squared =
            0.0;


        tree.closest_point(
            query,
            nearest_tid,
            closest,
            distance_squared);


        if (
            nearest_tid
            >=
            final_triangles.size())
        {
            throw std::runtime_error(
                "final Octree returned invalid triangle ID");
        }


        if (
            !std::isfinite(
                distance_squared)
            ||
            distance_squared
            <
            -1e-18)
        {
            throw std::runtime_error(
                "final Octree returned invalid squared distance");
        }


        const double distance =
            std::sqrt(
                std::max(
                    0.0,
                    distance_squared));


        const double normalized_distance =
            distance
            /
            sample.h;


        if (
            !std::isfinite(
                normalized_distance)
            ||
            normalized_distance < 0.0)
        {
            throw std::runtime_error(
                "invalid normalized Input->Final distance");
        }


        normalized.push_back(
            normalized_distance);


        //
        // Frozen:
        //
        // np.mean(normalized <= 1.0)
        //
        if (
            normalized_distance
            <=
            1.0)
        {
            ++coverage_1h_count;
        }
    }


    if (normalized.empty())
    {
        throw std::runtime_error(
            "no Input->Final samples");
    }


    const double coverage_1h =
        static_cast<double>(
            coverage_1h_count)
        /
        static_cast<double>(
            normalized.size());


    const double p99_h =
        percentile_linear_v1(
            normalized,
            99.0);


    const double tail_transform =
        1.0
        /
        std::max(
            1.0,
            p99_h);


    const double q_missing =
        std::sqrt(
            coverage_1h
            *
            tail_transform);


    require_finite(
        coverage_1h,
        "Q_missing.coverage_1h");

    require_finite(
        p99_h,
        "Q_missing.p99_h");

    require_finite(
        tail_transform,
        "Q_missing.tail");

    require_finite(
        q_missing,
        "Q_missing");


    MissingMetricsV1 out;

    out.final_triangle_count =
        final_triangles.size();

    out.sample_count =
        normalized.size();

    out.coverage_1h_count =
        coverage_1h_count;

    out.coverage_1h =
        coverage_1h;

    out.p99_h =
        p99_h;

    out.tail_transform =
        tail_transform;

    out.q_missing =
        q_missing;


    return out;
}


// ================================================================
// Active SHARP V1 geometry.
// ================================================================

struct SharpNormalV1
{
    bool valid;

    double x;
    double y;
    double z;
};


struct SharpEdgeRecordV1
{
    unsigned int v0;
    unsigned int v1;

    std::vector<std::size_t>
        incident_faces;
};


struct SharpAngleEdgeV1
{
    unsigned int v0;
    unsigned int v1;

    double angle;
};


static SharpNormalV1 sharp_newell_normal_v1(
    MetaMesh &mesh,
    const std::vector<unsigned int> &face)
{
    if (face.size() < 3)
    {
        throw std::runtime_error(
            "SHARP boundary face has fewer than 3 vertices");
    }


    double nx =
        0.0;

    double ny =
        0.0;

    double nz =
        0.0;


    //
    // Frozen sharp_metric_common_v1.py::face_normal().
    //
    for (
        std::size_t i = 0;
        i < face.size();
        ++i)
    {
        const vec3d &p =
            mesh.vert(
                face[
                    i
                ]);


        const vec3d &q =
            mesh.vert(
                face[
                    (
                        i + 1
                    )
                    %
                    face.size()
                ]);


        nx +=
            (
                p.y()
                -
                q.y()
            )
            *
            (
                p.z()
                +
                q.z()
            );


        ny +=
            (
                p.z()
                -
                q.z()
            )
            *
            (
                p.x()
                +
                q.x()
            );


        nz +=
            (
                p.x()
                -
                q.x()
            )
            *
            (
                p.y()
                +
                q.y()
            );
    }


    const double xx =
        nx * nx;

    const double yy =
        ny * ny;

    const double zz =
        nz * nz;


    //
    // Frozen NumPy 1.24.4 three-element norm:
    //
    //     sqrt((x*x + z*z) + y*y)
    //
    const double xz =
        xx
        +
        zz;


    const double length =
        std::sqrt(
            xz
            +
            yy);


    if (!std::isfinite(length))
    {
        throw std::runtime_error(
            "non-finite SHARP Newell normal length");
    }


    SharpNormalV1 result;

    result.valid =
        false;

    result.x =
        0.0;

    result.y =
        0.0;

    result.z =
        0.0;


    if (length <= 1e-18)
    {
        return result;
    }


    result.valid =
        true;

    result.x =
        nx
        /
        length;

    result.y =
        ny
        /
        length;

    result.z =
        nz
        /
        length;


    return result;
}


static std::vector<SharpAngleEdgeV1>
build_sharp_angle_edges_v1(
    MetaMesh &mesh,
    std::size_t &boundary_face_count,
    std::size_t &non_two_sided_edge_count)
{
    const std::vector<unsigned int>
        boundary_fids =
            mesh.get_surface_faces();


    if (boundary_fids.empty())
    {
        throw std::runtime_error(
            "SHARP final mesh has no boundary faces");
    }


    boundary_face_count =
        boundary_fids.size();


    std::vector<
        std::vector<unsigned int>
    > boundary_faces;


    boundary_faces.reserve(
        boundary_fids.size());


    std::vector<SharpNormalV1>
        normals;


    normals.reserve(
        boundary_fids.size());


    for (
        std::size_t bid = 0;
        bid < boundary_fids.size();
        ++bid)
    {
        const unsigned int fid =
            boundary_fids[
                bid
            ];


        const std::vector<unsigned int>
            face =
                mesh.face_verts_id(
                    fid,
                    false);


        if (face.size() < 3)
        {
            throw std::runtime_error(
                "SHARP boundary face has fewer than 3 vertices");
        }


        boundary_faces.push_back(
            face);


        normals.push_back(
            sharp_newell_normal_v1(
                mesh,
                face));
    }


    //
    // Frozen Python defaultdict insertion order:
    //
    // map performs lookup;
    // vector preserves first encounter order.
    //
    std::map<
        std::pair<
            unsigned int,
            unsigned int
        >,
        std::size_t
    > lookup;


    std::vector<SharpEdgeRecordV1>
        adjacency;


    for (
        std::size_t face_id = 0;
        face_id < boundary_faces.size();
        ++face_id)
    {
        const std::vector<unsigned int> &face =
            boundary_faces[
                face_id
            ];


        for (
            std::size_t i = 0;
            i < face.size();
            ++i)
        {
            const unsigned int a =
                face[
                    i
                ];


            const unsigned int b =
                face[
                    (
                        i + 1
                    )
                    %
                    face.size()
                ];


            if (a == b)
            {
                continue;
            }


            const unsigned int lo =
                std::min(
                    a,
                    b);


            const unsigned int hi =
                std::max(
                    a,
                    b);


            const std::pair<
                unsigned int,
                unsigned int
            > key(
                lo,
                hi);


            std::map<
                std::pair<
                    unsigned int,
                    unsigned int
                >,
                std::size_t
            >::iterator it =
                lookup.find(
                    key);


            if (
                it
                ==
                lookup.end())
            {
                SharpEdgeRecordV1 record;

                record.v0 =
                    lo;

                record.v1 =
                    hi;


                adjacency.push_back(
                    record);


                const std::size_t index =
                    adjacency.size()
                    -
                    1;


                lookup[
                    key
                ] = index;


                it =
                    lookup.find(
                        key);
            }


            adjacency[
                it->second
            ].incident_faces.push_back(
                face_id);
        }
    }


    non_two_sided_edge_count =
        0;


    std::vector<SharpAngleEdgeV1>
        angle_edges;


    //
    // Python math.degrees(x) conversion factor.
    //
    const double rad_to_deg =
        57.295779513082320876798154814105;


    for (
        std::size_t i = 0;
        i < adjacency.size();
        ++i)
    {
        const SharpEdgeRecordV1 &edge =
            adjacency[
                i
            ];


        if (
            edge.incident_faces.size()
            !=
            2)
        {
            ++non_two_sided_edge_count;

            continue;
        }


        const SharpNormalV1 &n0 =
            normals[
                edge.incident_faces[
                    0
                ]
            ];


        const SharpNormalV1 &n1 =
            normals[
                edge.incident_faces[
                    1
                ]
            ];


        if (
            !n0.valid
            ||
            !n1.valid)
        {
            continue;
        }


        const double m0 =
            n0.x
            *
            n1.x;


        const double m1 =
            n0.y
            *
            n1.y;


        const double m2 =
            n0.z
            *
            n1.z;


        //
        // Frozen NumPy 1.24.4 np.dot:
        //
        //     (m0 + m2) + m1
        //
        const double s02 =
            m0
            +
            m2;


        double dot =
            s02
            +
            m1;


        dot =
            std::fabs(
                dot);


        dot =
            std::min(
                1.0,
                std::max(
                    0.0,
                    dot));


        const double radians =
            std::acos(
                dot);


        const double angle =
            radians
            *
            rad_to_deg;


        if (!std::isfinite(angle))
        {
            throw std::runtime_error(
                "non-finite SHARP edge angle");
        }


        SharpAngleEdgeV1 result;

        result.v0 =
            edge.v0;

        result.v1 =
            edge.v1;

        result.angle =
            angle;


        angle_edges.push_back(
            result);
    }


    if (angle_edges.empty())
    {
        throw std::runtime_error(
            "no final two-sided SHARP edges");
    }


    return angle_edges;
}


// ================================================================
// Frozen scalar point-to-segment distance.
//
// NumPy 1.24.4 exact trees established by mechanical02:
//
// denom:
//     (x*x + z*z) + y*y
//
// numerator:
//     (rx*x + rz*z) + ry*y
//
// final norm:
//     sqrt((qx*qx + qy*qy) + qz*qz)
// ================================================================

static double sharp_point_segment_distance_v1(
    const SharpSample &sample,
    const vec3d &p0,
    const vec3d &p1)
{
    const double dx =
        p1.x()
        -
        p0.x();


    const double dy =
        p1.y()
        -
        p0.y();


    const double dz =
        p1.z()
        -
        p0.z();


    const double denom =
        (
            dx * dx
            +
            dz * dz
        )
        +
        dy * dy;


    const double rx =
        sample.x
        -
        p0.x();


    const double ry =
        sample.y
        -
        p0.y();


    const double rz =
        sample.z
        -
        p0.z();


    const double numerator =
        (
            rx * dx
            +
            rz * dz
        )
        +
        ry * dy;


    double t =
        0.0;


    if (denom > 0.0)
    {
        t =
            numerator
            /
            denom;
    }


    t =
        std::min(
            1.0,
            std::max(
                0.0,
                t));


    const double closest_x =
        p0.x()
        +
        t * dx;


    const double closest_y =
        p0.y()
        +
        t * dy;


    const double closest_z =
        p0.z()
        +
        t * dz;


    const double qx =
        closest_x
        -
        sample.x;


    const double qy =
        closest_y
        -
        sample.y;


    const double qz =
        closest_z
        -
        sample.z;


    const double distance_squared =
        (
            qx * qx
            +
            qy * qy
        )
        +
        qz * qz;


    if (
        !std::isfinite(
            distance_squared)
        ||
        distance_squared < 0.0)
    {
        throw std::runtime_error(
            "invalid SHARP point-segment squared distance");
    }


    const double distance =
        std::sqrt(
            distance_squared);


    require_finite(
        distance,
        "Q_sharp.segment_distance");


    return distance;
}


// ================================================================
// Q_sharp.
// ================================================================

SharpMetricsV1 compute_q_sharp_v1(
    MetaMesh &final_mesh,
    const QualityRef &ref)
{
    validate_quality_ref_v1(
        ref);


    if (ref.sharp_samples.empty())
    {
        throw std::runtime_error(
            "Q_sharp requested with no SHARP samples");
    }


    std::size_t boundary_face_count =
        0;


    std::size_t non_two_sided_edge_count =
        0;


    const std::vector<SharpAngleEdgeV1>
        edges =
            build_sharp_angle_edges_v1(
                final_mesh,
                boundary_face_count,
                non_two_sided_edge_count);


    std::vector<double>
        retention;


    retention.reserve(
        ref.sharp_samples.size());


    std::size_t covered75_count =
        0;


    std::size_t sample_edge_evaluations =
        0;


    for (
        std::size_t sample_id = 0;
        sample_id < ref.sharp_samples.size();
        ++sample_id)
    {
        const SharpSample &sample =
            ref.sharp_samples[
                sample_id
            ];


        if (
            !std::isfinite(
                sample.x)
            ||
            !std::isfinite(
                sample.y)
            ||
            !std::isfinite(
                sample.z)
            ||
            !std::isfinite(
                sample.h)
            ||
            !std::isfinite(
                sample.theta)
            ||
            sample.h <= 0.0
            ||
            sample.theta <= 0.0)
        {
            throw std::runtime_error(
                "invalid SHARP reference sample");
        }


        bool any_within =
            false;


        bool covered75 =
            false;


        double strongest =
            0.0;


        const double threshold75 =
            0.75
            *
            sample.theta;


        for (
            std::size_t edge_id = 0;
            edge_id < edges.size();
            ++edge_id)
        {
            const SharpAngleEdgeV1 &edge =
                edges[
                    edge_id
                ];


            const double distance =
                sharp_point_segment_distance_v1(
                    sample,
                    final_mesh.vert(
                        edge.v0),
                    final_mesh.vert(
                        edge.v1));


            ++sample_edge_evaluations;


            if (distance > sample.h)
            {
                continue;
            }


            any_within =
                true;


            if (
                edge.angle
                >
                strongest)
            {
                strongest =
                    edge.angle;
            }


            if (
                edge.angle
                >=
                threshold75)
            {
                covered75 =
                    true;
            }
        }


        double sample_retention =
            0.0;


        if (any_within)
        {
            sample_retention =
                std::min(
                    1.0,
                    strongest
                    /
                    sample.theta);
        }


        retention.push_back(
            sample_retention);


        if (covered75)
        {
            ++covered75_count;
        }
    }


    const std::size_t expected_evaluations =
        ref.sharp_samples.size()
        *
        edges.size();


    if (
        sample_edge_evaluations
        !=
        expected_evaluations)
    {
        throw std::runtime_error(
            "SHARP sample-edge evaluation count mismatch");
    }


    const double coverage75 =
        static_cast<double>(
            covered75_count)
        /
        static_cast<double>(
            ref.sharp_samples.size());


    //
    // Frozen np.mean(float64 retention):
    //
    //     NumPy-compatible sum
    //     /
    //     sample count
    //
    const double retention_sum =
        numpy_sum_v1(
            retention);


    const double retention_mean =
        retention_sum
        /
        static_cast<double>(
            retention.size());


    const double q_sharp =
        std::sqrt(
            coverage75
            *
            retention_mean);


    require_finite(
        coverage75,
        "Q_sharp.coverage75");


    require_finite(
        retention_mean,
        "Q_sharp.retention_mean");


    require_finite(
        q_sharp,
        "Q_sharp");


    SharpMetricsV1 out;

    out.boundary_face_count =
        boundary_face_count;

    out.two_sided_edge_count =
        edges.size();

    out.non_two_sided_edge_count =
        non_two_sided_edge_count;

    out.sharp_sample_count =
        ref.sharp_samples.size();

    out.sample_edge_evaluations =
        sample_edge_evaluations;

    out.covered75_count =
        covered75_count;

    out.coverage75 =
        coverage75;

    out.retention_mean =
        retention_mean;

    out.q_sharp =
        q_sharp;


    return out;
}


// ================================================================
// Q_fidelity V1.
// ================================================================

FidelityMetricsV1 compute_q_fidelity_v1(
    SrfMesh &input_surface,
    MetaMesh &final_mesh,
    const QualityRef &ref)
{
    validate_quality_ref_v1(
        ref);


    FidelityMetricsV1 out =
        {};


    //
    // Frozen geometry components.
    //
    out.spurious =
        compute_q_spurious_v1(
            input_surface,
            final_mesh,
            ref);


    out.missing =
        compute_q_missing_v1(
            final_mesh,
            ref);


    //
    // Frozen shape bottleneck rule:
    //
    //     min(Q_missing, Q_spurious)
    //
    // std::min returns the first operand on an exact tie,
    // matching Python min(a,b) first-item tie behavior.
    //
    out.q_shape =
        std::min(
            out.missing.q_missing,
            out.spurious.q_spurious);


    require_finite(
        out.q_shape,
        "Q_fidelity.q_shape");


    //
    // QUALITY_REF_V1 is the authority for whether
    // SHARP fidelity is active.
    //
    out.sharp_active =
        ref.sharp_present;


    out.sharp_metrics_valid =
        false;


    if (out.sharp_active)
    {
        if (ref.sharp_samples.empty())
        {
            throw std::runtime_error(
                "SHARP_PRESENT=1 requires SHARP samples");
        }


        out.sharp =
            compute_q_sharp_v1(
                final_mesh,
                ref);


        out.sharp_metrics_valid =
            true;


        //
        // Frozen active rule:
        //
        //     Q_fidelity =
        //         Q_shape * Q_sharp
        //
        out.q_fidelity =
            out.q_shape
            *
            out.sharp.q_sharp;
    }
    else
    {
        if (!ref.sharp_samples.empty())
        {
            throw std::runtime_error(
                "SHARP_PRESENT=0 requires zero SHARP samples");
        }


        //
        // Frozen inactive identity rule.
        //
        // Do not call compute_q_sharp_v1().
        //
        out.q_fidelity =
            out.q_shape;
    }


    require_finite(
        out.q_fidelity,
        "Q_fidelity");


    return out;
}

}
