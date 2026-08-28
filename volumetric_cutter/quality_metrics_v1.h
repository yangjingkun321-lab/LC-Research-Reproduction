#ifndef LOOPYCUTS_QUALITY_METRICS_V1_H
#define LOOPYCUTS_QUALITY_METRICS_V1_H

#include "definitions.h"
#include "quality_ref_v1.h"

#include <cstddef>
#include <cstdint>


namespace loopycuts_quality_v1
{

// ================================================================
// D_C V1
// ================================================================

struct DcMetricsV1
{
    std::uint64_t n_total;
    std::uint64_t n_hex;
    std::uint64_t n_nonhex;

    double h;
    double log1p_n_nonhex;
    double log1p_n_total;
    double log_ratio;
    double d_nonhex;
    double h_times_d_nonhex;
    double d_c;
};


DcMetricsV1 compute_dc_v1(
    std::uint64_t n_total,
    std::uint64_t n_hex);


// ================================================================
// Final -> Input / Q_spurious V1
// ================================================================

struct SpuriousMetricsV1
{
    std::size_t stage2_triangle_count;
    std::size_t final_triangle_count;
    std::size_t sample_count;

    std::size_t coverage_1h_count;

    double coverage_1h;
    double p99_h;

    double tail_transform;
    double q_spurious;
};


// Frozen Final -> Input metric.
//
// input_surface:
//     Exact Stage2 surface corresponding to
//     QualityRef::stage2_input_sha256.
//
// final_mesh:
//     Final post-subdivision/post-smoothing MetaMesh,
//     i.e. the geometry used by the frozen final-quality oracle.
//
// ref:
//     Frozen QUALITY_REF_V1 containing 30,000 FINAL_DRAWS.
//
// Contract:
//
//     Q_spurious = sqrt(
//         coverage_1h
//         *
//         (1 / max(1, p99_h))
//     )
//
// No floating-point tolerance is part of this API.
//
SpuriousMetricsV1 compute_q_spurious_v1(
    SrfMesh &input_surface,
    MetaMesh &final_mesh,
    const QualityRef &ref);


// ================================================================
// Input -> Final / Q_missing V1
// ================================================================

struct MissingMetricsV1
{
    std::size_t final_triangle_count;
    std::size_t sample_count;

    std::size_t coverage_1h_count;

    double coverage_1h;
    double p99_h;

    double tail_transform;
    double q_missing;
};


// Frozen Input -> Final metric.
//
// Query XYZ and local h come directly from the immutable
// QUALITY_REF_V1 INPUT_GEOMETRY section.
//
//     normalized_i =
//         distance(
//             input_sample_i,
//             final_boundary_surface)
//         / input_sample_i.h
//
//     Q_missing = sqrt(
//         coverage_1h
//         *
//         (1 / max(1, p99_h))
//     )
//
MissingMetricsV1 compute_q_missing_v1(
    MetaMesh &final_mesh,
    const QualityRef &ref);


// ================================================================
// Active SHARP / Q_sharp V1
// ================================================================

struct SharpMetricsV1
{
    std::size_t boundary_face_count;

    std::size_t two_sided_edge_count;
    std::size_t non_two_sided_edge_count;

    std::size_t sharp_sample_count;
    std::size_t sample_edge_evaluations;

    std::size_t covered75_count;

    double coverage75;
    double retention_mean;
    double q_sharp;
};


// Frozen active-SHARP metric.
//
// For each frozen SHARP sample:
//
//     inspect every valid two-sided final boundary edge
//
//     nearby(edge) iff
//         point_segment_distance <= sample.h
//
//     strongest =
//         max(angle(edge) for nearby edges)
//
//     retention =
//         0                              if no nearby edge
//         min(1, strongest/sample.theta) otherwise
//
//     covered75 iff
//         any nearby edge has
//         angle >= 0.75*sample.theta
//
// Final:
//
//     Q_sharp =
//         sqrt(coverage75 * retention_mean)
//
SharpMetricsV1 compute_q_sharp_v1(
    MetaMesh &final_mesh,
    const QualityRef &ref);


// ================================================================
// Combined shape / SHARP fidelity V1
// ================================================================

struct FidelityMetricsV1
{
    SpuriousMetricsV1 spurious;
    MissingMetricsV1 missing;

    double q_shape;

    bool sharp_active;
    bool sharp_metrics_valid;

    SharpMetricsV1 sharp;

    double q_fidelity;
};


// Frozen fidelity rule:
//
//     q_shape = min(
//         q_missing,
//         q_spurious)
//
// Active SHARP:
//
//     q_fidelity =
//         q_shape * q_sharp
//
// Inactive SHARP:
//
//     q_fidelity =
//         q_shape
//
// QualityRef::sharp_present is the authoritative
// active/inactive branch flag.
//
// sharp_metrics_valid is false for inactive SHARP;
// in that case the zero-initialized SharpMetricsV1
// payload is not semantically valid and must not be
// interpreted as a measured Q_sharp.
//
FidelityMetricsV1 compute_q_fidelity_v1(
    SrfMesh &input_surface,
    MetaMesh &final_mesh,
    const QualityRef &ref);

}

#endif
