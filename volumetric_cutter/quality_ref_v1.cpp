#include "quality_ref_v1.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace loopycuts_quality_v1
{

const char * const
    QUALITY_REF_V1_MAGIC =
        "LOOPYCUTS_QUALITY_REF_V1";

const char * const
    QUALITY_REF_V1_END =
        "END_LOOPYCUTS_QUALITY_REF_V1";

const char * const
    METRIC_CONTRACT_V3_SHA256 =
        "060d7f40293303aaa8fbdcbbe3999a303"
        "b14ee6b9dc0fd8bb93f21abed6d731d";

const std::size_t
    GEOMETRY_SAMPLE_COUNT_V1 =
        30000;

const std::size_t
    FINAL_DRAW_COUNT_V1 =
        30000;


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

QualityRef::QualityRef()
    : sharp_present(false)
    , has_sharp_declared_count(false)
    , sharp_declared_count(0)
    , has_sharp_file_sha256(false)
    , has_sharp_source_obj_sha256(false)
    , input_sample_seed_u64(0)
    , final_draw_seed_u64(0)
{
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static void fail(
    const std::string &message)
{
    throw std::runtime_error(
        message);
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Python reads the reference as UTF-8 before parsing.
//
// Canonical Train49 files are ASCII, but preserving the UTF-8 check
// keeps the C++ parser fail-closed with respect to the frozen contract.
//

static bool is_valid_utf8(
    const std::string &text)
{
    const unsigned char *s =
        reinterpret_cast<
            const unsigned char *>(
                text.data());

    const std::size_t n =
        text.size();

    std::size_t i = 0;

    while (i < n)
    {
        const unsigned char c =
            s[i];

        if (c <= 0x7f)
        {
            ++i;
            continue;
        }

        if (
            c >= 0xc2 &&
            c <= 0xdf)
        {
            if (
                i + 1 >= n ||
                s[i + 1] < 0x80 ||
                s[i + 1] > 0xbf)
            {
                return false;
            }

            i += 2;
            continue;
        }

        if (c == 0xe0)
        {
            if (
                i + 2 >= n ||
                s[i + 1] < 0xa0 ||
                s[i + 1] > 0xbf ||
                s[i + 2] < 0x80 ||
                s[i + 2] > 0xbf)
            {
                return false;
            }

            i += 3;
            continue;
        }

        if (
            c >= 0xe1 &&
            c <= 0xec)
        {
            if (
                i + 2 >= n ||
                s[i + 1] < 0x80 ||
                s[i + 1] > 0xbf ||
                s[i + 2] < 0x80 ||
                s[i + 2] > 0xbf)
            {
                return false;
            }

            i += 3;
            continue;
        }

        if (c == 0xed)
        {
            if (
                i + 2 >= n ||
                s[i + 1] < 0x80 ||
                s[i + 1] > 0x9f ||
                s[i + 2] < 0x80 ||
                s[i + 2] > 0xbf)
            {
                return false;
            }

            i += 3;
            continue;
        }

        if (
            c >= 0xee &&
            c <= 0xef)
        {
            if (
                i + 2 >= n ||
                s[i + 1] < 0x80 ||
                s[i + 1] > 0xbf ||
                s[i + 2] < 0x80 ||
                s[i + 2] > 0xbf)
            {
                return false;
            }

            i += 3;
            continue;
        }

        if (c == 0xf0)
        {
            if (
                i + 3 >= n ||
                s[i + 1] < 0x90 ||
                s[i + 1] > 0xbf ||
                s[i + 2] < 0x80 ||
                s[i + 2] > 0xbf ||
                s[i + 3] < 0x80 ||
                s[i + 3] > 0xbf)
            {
                return false;
            }

            i += 4;
            continue;
        }

        if (
            c >= 0xf1 &&
            c <= 0xf3)
        {
            if (
                i + 3 >= n ||
                s[i + 1] < 0x80 ||
                s[i + 1] > 0xbf ||
                s[i + 2] < 0x80 ||
                s[i + 2] > 0xbf ||
                s[i + 3] < 0x80 ||
                s[i + 3] > 0xbf)
            {
                return false;
            }

            i += 4;
            continue;
        }

        if (c == 0xf4)
        {
            if (
                i + 3 >= n ||
                s[i + 1] < 0x80 ||
                s[i + 1] > 0x8f ||
                s[i + 2] < 0x80 ||
                s[i + 2] > 0xbf ||
                s[i + 3] < 0x80 ||
                s[i + 3] > 0xbf)
            {
                return false;
            }

            i += 4;
            continue;
        }

        return false;
    }

    return true;
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static std::string read_file_bytes(
    const std::string &path)
{
    std::ifstream in(
        path.c_str(),
        std::ios::in |
        std::ios::binary);

    if (!in)
    {
        fail(
            "Cannot open quality reference: "
            +
            path);
    }

    std::ostringstream buffer;

    buffer
        << in.rdbuf();

    if (
        in.bad())
    {
        fail(
            "Failed while reading quality reference: "
            +
            path);
    }

    return buffer.str();
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Mimics Python text.splitlines() for LF-only canonical V1 data.
//
// A single final newline does not produce an additional logical line.
//

static std::vector<std::string>
split_lines_lf(
    const std::string &raw)
{
    std::vector<std::string>
        lines;

    std::size_t start = 0;

    while (
        start < raw.size())
    {
        const std::size_t pos =
            raw.find(
                '\n',
                start);

        if (
            pos ==
            std::string::npos)
        {
            lines.push_back(
                raw.substr(
                    start));

            break;
        }

        lines.push_back(
            raw.substr(
                start,
                pos - start));

        start =
            pos + 1;
    }

    return lines;
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static bool is_lower_sha256(
    const std::string &text)
{
    if (
        text.size() != 64)
    {
        return false;
    }

    for (
        std::size_t i = 0;
        i < text.size();
        ++i)
    {
        const char c =
            text[i];

        const bool digit =
            (
                c >= '0' &&
                c <= '9'
            );

        const bool lower_hex =
            (
                c >= 'a' &&
                c <= 'f'
            );

        if (
            !digit &&
            !lower_hex)
        {
            return false;
        }
    }

    return true;
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static void validate_sha256(
    const std::string &value,
    const std::string &name)
{
    if (
        !is_lower_sha256(
            value))
    {
        fail(
            name
            +
            " must be lowercase SHA256");
    }
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static bool
is_python_unicode_whitespace(
    const std::uint32_t cp)
{
    //
    // Python 3.11 str.isspace() whitespace code points.
    //
    // This includes:
    //
    //   U+0009..U+000D
    //   U+001C..U+001F
    //   U+0020
    //   U+0085
    //   U+00A0
    //   U+1680
    //   U+2000..U+200A
    //   U+2028
    //   U+2029
    //   U+202F
    //   U+205F
    //   U+3000
    //
    return
        (
            cp >= 0x0009 &&
            cp <= 0x000d
        )
        ||
        (
            cp >= 0x001c &&
            cp <= 0x001f
        )
        ||
        cp == 0x0020
        ||
        cp == 0x0085
        ||
        cp == 0x00a0
        ||
        cp == 0x1680
        ||
        (
            cp >= 0x2000 &&
            cp <= 0x200a
        )
        ||
        cp == 0x2028
        ||
        cp == 0x2029
        ||
        cp == 0x202f
        ||
        cp == 0x205f
        ||
        cp == 0x3000;
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static std::uint32_t
decode_valid_utf8_codepoint(
    const std::string &text,
    std::size_t &offset)
{
    const unsigned char *s =
        reinterpret_cast<
            const unsigned char *>(
                text.data());

    const unsigned char c =
        s[offset];

    if (c <= 0x7f)
    {
        ++offset;

        return
            static_cast<
                std::uint32_t>(
                    c);
    }

    if (c <= 0xdf)
    {
        const std::uint32_t cp =
            (
                static_cast<
                    std::uint32_t>(
                        c & 0x1f)
                << 6
            )
            |
            static_cast<
                std::uint32_t>(
                    s[offset + 1]
                    &
                    0x3f);

        offset += 2;

        return cp;
    }

    if (c <= 0xef)
    {
        const std::uint32_t cp =
            (
                static_cast<
                    std::uint32_t>(
                        c & 0x0f)
                << 12
            )
            |
            (
                static_cast<
                    std::uint32_t>(
                        s[offset + 1]
                        &
                        0x3f)
                << 6
            )
            |
            static_cast<
                std::uint32_t>(
                    s[offset + 2]
                    &
                    0x3f);

        offset += 3;

        return cp;
    }

    const std::uint32_t cp =
        (
            static_cast<
                std::uint32_t>(
                    c & 0x07)
            << 18
        )
        |
        (
            static_cast<
                std::uint32_t>(
                    s[offset + 1]
                    &
                    0x3f)
            << 12
        )
        |
        (
            static_cast<
                std::uint32_t>(
                    s[offset + 2]
                    &
                    0x3f)
            << 6
        )
        |
        static_cast<
            std::uint32_t>(
                s[offset + 3]
                &
                0x3f);

    offset += 4;

    return cp;
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static void validate_model(
    const std::string &model)
{
    if (
        model.empty())
    {
        fail(
            "model must not be empty");
    }

    //
    // Python QualityRefV1 validation uses:
    //
    //     any(ch.isspace() for ch in model)
    //
    // Preserve that semantics for UTF-8 model names.
    //
    if (
        !is_valid_utf8(
            model))
    {
        fail(
            "model must be valid UTF-8");
    }

    std::size_t offset = 0;

    while (
        offset <
        model.size())
    {
        const std::uint32_t cp =
            decode_valid_utf8_codepoint(
                model,
                offset);

        if (
            is_python_unicode_whitespace(
                cp))
        {
            fail(
                "model must not contain whitespace");
        }
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static void require_finite(
    const double value,
    const std::string &name)
{
    if (
        !std::isfinite(
            value))
    {
        fail(
            name
            +
            " must be finite");
    }
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Canonical Python V1 uses:
//
//     format(value, ".17g")
//
// libstdc++ defaultfloat + precision(17) follows the same %g-style
// representation for IEEE double on the target Linux environment.
//
// We verify this assumption against all 49 frozen references before
// this parser is ever integrated into the runtime.
//

static std::string
format_double_canonical(
    const double value)
{
    require_finite(
        value,
        "canonical float");

    if (
        value == 0.0)
    {
        return "0";
    }

    std::ostringstream out;

    out.imbue(
        std::locale::classic());

    out
        << std::setprecision(17)
        << std::defaultfloat
        << value;

    return out.str();
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static std::uint64_t
parse_uint64_token(
    const std::string &token,
    const std::string &name)
{
    if (
        token.empty())
    {
        fail(
            name
            +
            " is empty");
    }

    if (
        token.size() > 1 &&
        token[0] == '0')
    {
        fail(
            name
            +
            " is not canonical decimal");
    }

    std::uint64_t value = 0;

    const std::uint64_t max_value =
        std::numeric_limits<
            std::uint64_t>::max();

    for (
        std::size_t i = 0;
        i < token.size();
        ++i)
    {
        const char c =
            token[i];

        if (
            c < '0' ||
            c > '9')
        {
            fail(
                name
                +
                " must be unsigned decimal");
        }

        const std::uint64_t digit =
            static_cast<
                std::uint64_t>(
                    c - '0');

        if (
            value >
            (
                max_value - digit
            ) / 10)
        {
            fail(
                name
                +
                " exceeds uint64 range");
        }

        value =
            value * 10
            +
            digit;
    }

    return value;
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static std::size_t
parse_size_token(
    const std::string &token,
    const std::string &name)
{
    const std::uint64_t value =
        parse_uint64_token(
            token,
            name);

    if (
        value >
        static_cast<
            std::uint64_t>(
                std::numeric_limits<
                    std::size_t>::max()))
    {
        fail(
            name
            +
            " exceeds size_t range");
    }

    return
        static_cast<
            std::size_t>(
                value);
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static double parse_double_token(
    const std::string &token,
    const std::string &name)
{
    if (
        token.empty())
    {
        fail(
            name
            +
            " is empty");
    }

    std::istringstream in(
        token);

    in.imbue(
        std::locale::classic());

    double value = 0.0;

    in
        >> std::noskipws
        >> value;

    if (
        !in ||
        !in.eof())
    {
        fail(
            name
            +
            " is not a valid double");
    }

    require_finite(
        value,
        name);

    return value;
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// Canonical rows contain exactly one ASCII space between tokens.
//

static std::vector<std::string>
split_single_space(
    const std::string &line,
    const std::size_t expected,
    const std::string &context)
{
    if (
        line.empty())
    {
        fail(
            context
            +
            ": empty row");
    }

    std::vector<std::string>
        fields;

    std::size_t start = 0;

    while (true)
    {
        const std::size_t pos =
            line.find(
                ' ',
                start);

        if (
            pos ==
            std::string::npos)
        {
            fields.push_back(
                line.substr(
                    start));

            break;
        }

        if (
            pos == start)
        {
            fail(
                context
                +
                ": non-canonical spacing");
        }

        fields.push_back(
            line.substr(
                start,
                pos - start));

        start =
            pos + 1;

        if (
            start ==
            line.size())
        {
            fail(
                context
                +
                ": trailing space");
        }
    }

    if (
        fields.size() !=
        expected)
    {
        std::ostringstream out;

        out
            << context
            << ": expected "
            << expected
            << " fields; got "
            << fields.size();

        fail(
            out.str());
    }

    return fields;
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

static std::string
parse_header_value(
    const std::string &line,
    const std::string &key)
{
    const std::string prefix =
        key
        +
        " ";

    if (
        line.size() <=
            prefix.size() ||
        line.compare(
            0,
            prefix.size(),
            prefix) != 0)
    {
        fail(
            "Expected header "
            +
            key
            +
            "; got: "
            +
            line);
    }

    return
        line.substr(
            prefix.size());
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void validate_quality_ref_v1(
    const QualityRef &ref)
{
    validate_model(
        ref.model);

    validate_sha256(
        ref.metric_contract_sha256,
        "metric_contract_sha256");

    validate_sha256(
        ref.stage2_input_sha256,
        "stage2_input_sha256");

    if (
        ref.metric_contract_sha256
        !=
        METRIC_CONTRACT_V3_SHA256)
    {
        fail(
            "Unsupported metric-contract SHA256: "
            +
            ref.metric_contract_sha256);
    }


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // SHARP provenance semantics.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (
        !ref.has_sharp_declared_count)
    {
        if (
            ref.sharp_present)
        {
            fail(
                "sharp_present=True requires "
                "an explicit SHARP source");
        }

        if (
            ref.has_sharp_file_sha256)
        {
            fail(
                "No SHARP source requires "
                "sharp_file_sha256=NONE");
        }

        if (
            ref.has_sharp_source_obj_sha256)
        {
            fail(
                "No SHARP source requires "
                "sharp_source_obj_sha256=NONE");
        }

        if (
            !ref.sharp_samples.empty())
        {
            fail(
                "No SHARP source requires "
                "zero SHARP samples");
        }
    }
    else
    {
        if (
            !ref.has_sharp_file_sha256)
        {
            fail(
                "Explicit SHARP source requires "
                "sharp_file_sha256");
        }

        if (
            !ref.has_sharp_source_obj_sha256)
        {
            fail(
                "Explicit SHARP source requires "
                "sharp_source_obj_sha256");
        }

        validate_sha256(
            ref.sharp_file_sha256,
            "sharp_file_sha256");

        validate_sha256(
            ref.sharp_source_obj_sha256,
            "sharp_source_obj_sha256");

        if (
            ref.sharp_declared_count == 0)
        {
            if (
                ref.sharp_present)
            {
                fail(
                    "SHARP_PRESENT must be false "
                    "when declared count is zero");
            }

            if (
                !ref.sharp_samples.empty())
            {
                fail(
                    "Zero-feature SHARP source "
                    "requires zero SHARP samples");
            }
        }
        else
        {
            if (
                !ref.sharp_present)
            {
                fail(
                    "Positive SHARP declared count "
                    "requires SHARP_PRESENT=1");
            }

            if (
                ref.sharp_samples.empty())
            {
                fail(
                    "Positive SHARP declared count "
                    "requires SHARP samples");
            }
        }
    }


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // V1 fixed geometry-draw counts.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (
        ref.input_geometry.size()
        !=
        GEOMETRY_SAMPLE_COUNT_V1)
    {
        fail(
            "V1 requires exactly 30000 "
            "Input geometry samples");
    }

    if (
        ref.final_draws.size()
        !=
        FINAL_DRAW_COUNT_V1)
    {
        fail(
            "V1 requires exactly 30000 "
            "Final surface draws");
    }


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Input geometry.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    for (
        std::size_t i = 0;
        i < ref.input_geometry.size();
        ++i)
    {
        const InputGeometrySample &s =
            ref.input_geometry[i];

        require_finite(
            s.x,
            "input_geometry.x");

        require_finite(
            s.y,
            "input_geometry.y");

        require_finite(
            s.z,
            "input_geometry.z");

        require_finite(
            s.h,
            "input_geometry.h");

        if (
            s.h <= 0.0)
        {
            fail(
                "input_geometry.h "
                "must be positive");
        }
    }


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Final draws.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    for (
        std::size_t i = 0;
        i < ref.final_draws.size();
        ++i)
    {
        const FinalSurfaceDraw &d =
            ref.final_draws[i];

        require_finite(
            d.area_draw,
            "final_draw.area_draw");

        require_finite(
            d.u,
            "final_draw.u");

        require_finite(
            d.v,
            "final_draw.v");

        if (
            d.area_draw < 0.0 ||
            d.area_draw >= 1.0 ||
            d.u < 0.0 ||
            d.u >= 1.0 ||
            d.v < 0.0 ||
            d.v >= 1.0)
        {
            fail(
                "Final draw values must satisfy "
                "0 <= x < 1");
        }
    }


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Continuous SHARP samples.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    for (
        std::size_t i = 0;
        i < ref.sharp_samples.size();
        ++i)
    {
        const SharpSample &s =
            ref.sharp_samples[i];

        require_finite(
            s.x,
            "sharp.x");

        require_finite(
            s.y,
            "sharp.y");

        require_finite(
            s.z,
            "sharp.z");

        require_finite(
            s.h,
            "sharp.h");

        require_finite(
            s.theta,
            "sharp.theta");

        if (
            s.h <= 0.0)
        {
            fail(
                "sharp.h must be positive");
        }

        if (
            s.theta <= 0.0 ||
            s.theta > 90.0)
        {
            fail(
                "sharp.theta must satisfy "
                "0 < theta <= 90");
        }
    }
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::string quality_ref_v1_to_text(
    const QualityRef &ref)
{
    validate_quality_ref_v1(
        ref);

    std::ostringstream out;

    out.imbue(
        std::locale::classic());

    out
        << QUALITY_REF_V1_MAGIC
        << "\n"
        << "\n"
        << "MODEL "
        << ref.model
        << "\n"
        << "METRIC_CONTRACT_SHA256 "
        << ref.metric_contract_sha256
        << "\n"
        << "STAGE2_INPUT_SHA256 "
        << ref.stage2_input_sha256
        << "\n"
        << "SHARP_PRESENT "
        << (
            ref.sharp_present
            ?
            "1"
            :
            "0"
        )
        << "\n"
        << "SHARP_DECLARED_COUNT ";

    if (
        ref.has_sharp_declared_count)
    {
        out
            << ref.sharp_declared_count;
    }
    else
    {
        out
            << "NONE";
    }

    out
        << "\n"
        << "SHARP_FILE_SHA256 ";

    if (
        ref.has_sharp_file_sha256)
    {
        out
            << ref.sharp_file_sha256;
    }
    else
    {
        out
            << "NONE";
    }

    out
        << "\n"
        << "SHARP_SOURCE_OBJ_SHA256 ";

    if (
        ref.has_sharp_source_obj_sha256)
    {
        out
            << ref.sharp_source_obj_sha256;
    }
    else
    {
        out
            << "NONE";
    }

    out
        << "\n"
        << "INPUT_SAMPLE_SEED_U64 "
        << ref.input_sample_seed_u64
        << "\n"
        << "FINAL_DRAW_SEED_U64 "
        << ref.final_draw_seed_u64
        << "\n"
        << "GEOMETRY_SAMPLE_COUNT "
        << ref.input_geometry.size()
        << "\n"
        << "FINAL_DRAW_COUNT "
        << ref.final_draws.size()
        << "\n"
        << "SHARP_SAMPLE_COUNT "
        << ref.sharp_samples.size()
        << "\n"
        << "\n"
        << "BEGIN_INPUT_GEOMETRY"
        << "\n";

    for (
        std::size_t i = 0;
        i < ref.input_geometry.size();
        ++i)
    {
        const InputGeometrySample &s =
            ref.input_geometry[i];

        out
            << format_double_canonical(
                   s.x)
            << " "
            << format_double_canonical(
                   s.y)
            << " "
            << format_double_canonical(
                   s.z)
            << " "
            << format_double_canonical(
                   s.h)
            << "\n";
    }

    out
        << "END_INPUT_GEOMETRY"
        << "\n"
        << "\n"
        << "BEGIN_FINAL_DRAWS"
        << "\n";

    for (
        std::size_t i = 0;
        i < ref.final_draws.size();
        ++i)
    {
        const FinalSurfaceDraw &d =
            ref.final_draws[i];

        out
            << format_double_canonical(
                   d.area_draw)
            << " "
            << format_double_canonical(
                   d.u)
            << " "
            << format_double_canonical(
                   d.v)
            << "\n";
    }

    out
        << "END_FINAL_DRAWS"
        << "\n"
        << "\n"
        << "BEGIN_SHARP"
        << "\n";

    for (
        std::size_t i = 0;
        i < ref.sharp_samples.size();
        ++i)
    {
        const SharpSample &s =
            ref.sharp_samples[i];

        out
            << format_double_canonical(
                   s.x)
            << " "
            << format_double_canonical(
                   s.y)
            << " "
            << format_double_canonical(
                   s.z)
            << " "
            << format_double_canonical(
                   s.h)
            << " "
            << format_double_canonical(
                   s.theta)
            << "\n";
    }

    out
        << "END_SHARP"
        << "\n"
        << "\n"
        << QUALITY_REF_V1_END
        << "\n";

    return out.str();
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

QualityRef read_quality_ref_v1(
    const std::string &path)
{
    const std::string raw =
        read_file_bytes(
            path);

    if (
        !is_valid_utf8(
            raw))
    {
        fail(
            path
            +
            ": not valid UTF-8");
    }

    if (
        raw.find('\r')
        !=
        std::string::npos)
    {
        fail(
            path
            +
            ": CR bytes are not canonical");
    }

    const std::vector<std::string>
        lines =
            split_lines_lf(
                raw);

    std::size_t cursor = 0;


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Local sequential-reader helpers.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    const auto take =
        [&]() -> std::string
        {
            if (
                cursor >=
                lines.size())
            {
                fail(
                    path
                    +
                    ": unexpected EOF");
            }

            return
                lines[
                    cursor++
                ];
        };

    const auto expect =
        [&](const std::string &value)
        {
            const std::string line =
                take();

            if (
                line !=
                value)
            {
                fail(
                    path
                    +
                    ": expected '"
                    +
                    value
                    +
                    "'; got '"
                    +
                    line
                    +
                    "'");
            }
        };


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Header.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    expect(
        QUALITY_REF_V1_MAGIC);

    expect(
        "");

    QualityRef ref;

    ref.model =
        parse_header_value(
            take(),
            "MODEL");

    ref.metric_contract_sha256 =
        parse_header_value(
            take(),
            "METRIC_CONTRACT_SHA256");

    ref.stage2_input_sha256 =
        parse_header_value(
            take(),
            "STAGE2_INPUT_SHA256");


    const std::string
        sharp_present =
            parse_header_value(
                take(),
                "SHARP_PRESENT");

    if (
        sharp_present == "0")
    {
        ref.sharp_present =
            false;
    }
    else if (
        sharp_present == "1")
    {
        ref.sharp_present =
            true;
    }
    else
    {
        fail(
            "SHARP_PRESENT must be 0 or 1");
    }


    const std::string
        sharp_declared =
            parse_header_value(
                take(),
                "SHARP_DECLARED_COUNT");

    if (
        sharp_declared == "NONE")
    {
        ref.has_sharp_declared_count =
            false;

        ref.sharp_declared_count =
            0;
    }
    else
    {
        ref.has_sharp_declared_count =
            true;

        ref.sharp_declared_count =
            parse_uint64_token(
                sharp_declared,
                "SHARP_DECLARED_COUNT");
    }


    const std::string
        sharp_file_sha =
            parse_header_value(
                take(),
                "SHARP_FILE_SHA256");

    if (
        sharp_file_sha == "NONE")
    {
        ref.has_sharp_file_sha256 =
            false;

        ref.sharp_file_sha256.clear();
    }
    else
    {
        ref.has_sharp_file_sha256 =
            true;

        ref.sharp_file_sha256 =
            sharp_file_sha;
    }


    const std::string
        sharp_source_sha =
            parse_header_value(
                take(),
                "SHARP_SOURCE_OBJ_SHA256");

    if (
        sharp_source_sha == "NONE")
    {
        ref.has_sharp_source_obj_sha256 =
            false;

        ref.sharp_source_obj_sha256.clear();
    }
    else
    {
        ref.has_sharp_source_obj_sha256 =
            true;

        ref.sharp_source_obj_sha256 =
            sharp_source_sha;
    }


    ref.input_sample_seed_u64 =
        parse_uint64_token(
            parse_header_value(
                take(),
                "INPUT_SAMPLE_SEED_U64"),
            "INPUT_SAMPLE_SEED_U64");


    ref.final_draw_seed_u64 =
        parse_uint64_token(
            parse_header_value(
                take(),
                "FINAL_DRAW_SEED_U64"),
            "FINAL_DRAW_SEED_U64");


    const std::size_t
        geometry_sample_count =
            parse_size_token(
                parse_header_value(
                    take(),
                    "GEOMETRY_SAMPLE_COUNT"),
                "GEOMETRY_SAMPLE_COUNT");


    const std::size_t
        final_draw_count =
            parse_size_token(
                parse_header_value(
                    take(),
                    "FINAL_DRAW_COUNT"),
                "FINAL_DRAW_COUNT");


    const std::size_t
        sharp_sample_count =
            parse_size_token(
                parse_header_value(
                    take(),
                    "SHARP_SAMPLE_COUNT"),
                "SHARP_SAMPLE_COUNT");


    //
    // Fail before reading or allocating large bodies.
    //
    if (
        geometry_sample_count !=
        GEOMETRY_SAMPLE_COUNT_V1)
    {
        fail(
            "V1 GEOMETRY_SAMPLE_COUNT "
            "must equal 30000");
    }

    if (
        final_draw_count !=
        FINAL_DRAW_COUNT_V1)
    {
        fail(
            "V1 FINAL_DRAW_COUNT "
            "must equal 30000");
    }


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // INPUT_GEOMETRY.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    expect(
        "");

    expect(
        "BEGIN_INPUT_GEOMETRY");

    ref.input_geometry.reserve(
        geometry_sample_count);

    for (
        std::size_t i = 0;
        i < geometry_sample_count;
        ++i)
    {
        const std::vector<std::string>
            fields =
                split_single_space(
                    take(),
                    4,
                    "INPUT_GEOMETRY");

        InputGeometrySample s;

        s.x =
            parse_double_token(
                fields[0],
                "INPUT_GEOMETRY.x");

        s.y =
            parse_double_token(
                fields[1],
                "INPUT_GEOMETRY.y");

        s.z =
            parse_double_token(
                fields[2],
                "INPUT_GEOMETRY.z");

        s.h =
            parse_double_token(
                fields[3],
                "INPUT_GEOMETRY.h");

        ref.input_geometry.push_back(
            s);
    }

    expect(
        "END_INPUT_GEOMETRY");


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // FINAL_DRAWS.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    expect(
        "");

    expect(
        "BEGIN_FINAL_DRAWS");

    ref.final_draws.reserve(
        final_draw_count);

    for (
        std::size_t i = 0;
        i < final_draw_count;
        ++i)
    {
        const std::vector<std::string>
            fields =
                split_single_space(
                    take(),
                    3,
                    "FINAL_DRAWS");

        FinalSurfaceDraw d;

        d.area_draw =
            parse_double_token(
                fields[0],
                "FINAL_DRAWS.area_draw");

        d.u =
            parse_double_token(
                fields[1],
                "FINAL_DRAWS.u");

        d.v =
            parse_double_token(
                fields[2],
                "FINAL_DRAWS.v");

        ref.final_draws.push_back(
            d);
    }

    expect(
        "END_FINAL_DRAWS");


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // SHARP.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    expect(
        "");

    expect(
        "BEGIN_SHARP");

    for (
        std::size_t i = 0;
        i < sharp_sample_count;
        ++i)
    {
        const std::vector<std::string>
            fields =
                split_single_space(
                    take(),
                    5,
                    "SHARP");

        SharpSample s;

        s.x =
            parse_double_token(
                fields[0],
                "SHARP.x");

        s.y =
            parse_double_token(
                fields[1],
                "SHARP.y");

        s.z =
            parse_double_token(
                fields[2],
                "SHARP.z");

        s.h =
            parse_double_token(
                fields[3],
                "SHARP.h");

        s.theta =
            parse_double_token(
                fields[4],
                "SHARP.theta");

        ref.sharp_samples.push_back(
            s);
    }

    expect(
        "END_SHARP");

    expect(
        "");

    expect(
        QUALITY_REF_V1_END);


    if (
        cursor !=
        lines.size())
    {
        fail(
            path
            +
            ": trailing content after "
            +
            QUALITY_REF_V1_END);
    }


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Semantic validation.
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    validate_quality_ref_v1(
        ref);


    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // Canonical byte equality.
    //
    // This is intentionally stronger than merely accepting equivalent
    // numeric data. It mirrors Python:
    //
    //     quality_ref_v1_to_text(ref).encode("utf-8") == raw
    //
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    const std::string canonical =
        quality_ref_v1_to_text(
            ref);

    if (
        canonical !=
        raw)
    {
        fail(
            path
            +
            ": valid data but bytes are not "
            "canonical V1 encoding");
    }

    return ref;
}

}
