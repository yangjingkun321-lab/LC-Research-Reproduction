#ifndef LOOPYCUTS_QUALITY_REF_V1_H
#define LOOPYCUTS_QUALITY_REF_V1_H

#include <cstdint>
#include <string>
#include <vector>

namespace loopycuts_quality_v1
{

struct InputGeometrySample
{
    double x;
    double y;
    double z;
    double h;
};

struct FinalSurfaceDraw
{
    double area_draw;
    double u;
    double v;
};

struct SharpSample
{
    double x;
    double y;
    double z;
    double h;
    double theta;
};

struct QualityRef
{
    std::string model;

    std::string metric_contract_sha256;
    std::string stage2_input_sha256;

    bool sharp_present;

    bool has_sharp_declared_count;
    std::uint64_t sharp_declared_count;

    bool has_sharp_file_sha256;
    std::string sharp_file_sha256;

    bool has_sharp_source_obj_sha256;
    std::string sharp_source_obj_sha256;

    std::uint64_t input_sample_seed_u64;
    std::uint64_t final_draw_seed_u64;

    std::vector<InputGeometrySample>
        input_geometry;

    std::vector<FinalSurfaceDraw>
        final_draws;

    std::vector<SharpSample>
        sharp_samples;

    QualityRef();
};

extern const char * const
    QUALITY_REF_V1_MAGIC;

extern const char * const
    QUALITY_REF_V1_END;

extern const char * const
    METRIC_CONTRACT_V3_SHA256;

extern const std::size_t
    GEOMETRY_SAMPLE_COUNT_V1;

extern const std::size_t
    FINAL_DRAW_COUNT_V1;


// Validate semantic V1 invariants.
void validate_quality_ref_v1(
    const QualityRef &ref);


// Serialize using the same canonical layout used by
// frozen Python quality_ref_v1.
std::string quality_ref_v1_to_text(
    const QualityRef &ref);


// Read, validate, canonicalize, and require exact
// byte equality with the canonical representation.
QualityRef read_quality_ref_v1(
    const std::string &path);

}

#endif
