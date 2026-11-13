// How PTM3D's DEM tools turn stored raster values into heights in meters.
//
// C++ mirror of scripts/dem_elevation.py (read its docstring for the rule).
// Used by dtm_to_ply_delaunay and dem_mesh_generator; checked against the same
// case table as the Python version (scripts/dem_elevation_cases.tsv) by
// test_dem_elevation, so the tools cannot disagree.
//
//     height_m = stored_value * scale + offset
//
// No GDAL here: dem_elevation_gdal.h reads the metadata this rule consumes.
#pragma once

#include <cctype>
#include <cmath>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ptm {
namespace dem {

constexpr double REF_RADIUS_MIN_M = 100000.0;  // an offset above this is a radius, not a height
constexpr double REF_RADIUS_TOL_M = 20000.0;   // how far a reference radius may sit from the body radius
constexpr double FILL_ABS_MIN = 1e37;          // float fill values near +-FLT_MAX

struct DemElevationError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct ElevationTransform {
    double scale = 1.0;
    double offset = 0.0;
    std::string scale_source;
    std::string offset_source;
    std::vector<std::string> warnings;

    std::string describe() const {
        std::ostringstream s;
        s << "height (m) = stored value * " << scale << " [" << scale_source << "] + "
          << offset << " [" << offset_source << "]";
        return s.str();
    }
};

inline std::string upperTrim(const std::string& in) {
    size_t a = in.find_first_not_of(" \t\r\n");
    size_t b = in.find_last_not_of(" \t\r\n");
    std::string s = (a == std::string::npos) ? std::string() : in.substr(a, b - a + 1);
    for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

inline double unitFactor(const std::string& unit) {
    const std::string k = upperTrim(unit);
    if (k.empty() || k == "M" || k == "METER" || k == "METERS")
        return 1.0;
    if (k == "KM" || k == "KILOMETER" || k == "KILOMETERS")
        return 1000.0;
    throw DemElevationError("unknown elevation unit '" + unit + "' in the file metadata; "
                            "pass the scale explicitly (meters per stored unit)");
}

// UNIT keyword of a PDS label, from GDAL's json:PDS metadata domain ("" if none).
inline std::string pdsUnitFromJson(const char* text) {
    if (!text) return "";
    static const std::regex re("\"UNIT\"\\s*:\\s*\"([^\"]*)\"");
    std::cmatch m;
    return std::regex_search(text, m, re) ? m[1].str() : "";
}

inline std::string fmt(double v) {
    std::ostringstream s;
    s << v;
    return s.str();
}

// Apply the rule. std::nullopt means "absent". Throws DemElevationError.
inline ElevationTransform resolve(std::optional<double> meta_scale,
                                  std::optional<double> meta_offset,
                                  const std::string& meta_unit, bool is_integer,
                                  std::optional<double> user_scale,
                                  std::optional<double> user_offset,
                                  double body_radius_m) {
    auto u = [&]() { return unitFactor(meta_unit); };  // throws only where the unit is needed
    const std::string unit_up = upperTrim(meta_unit);
    const std::string unit_note = unit_up.empty() ? "" : ", unit " + unit_up;
    if (meta_scale && *meta_scale == 0.0) meta_scale.reset();    // a zero scale is not a scale
    if (meta_offset && *meta_offset == 0.0) meta_offset.reset();
    ElevationTransform t;

    if (user_scale) {
        t.scale = *user_scale;
        t.scale_source = "command line";
    } else if (meta_scale) {
        t.scale = *meta_scale * u();
        t.scale_source = "file metadata" + unit_note;
    } else if (!is_integer) {
        t.scale = 1.0 * u();
        t.scale_source = "default, no scale metadata" + unit_note;
    } else {
        throw DemElevationError(
            "integer raster with no scale metadata: the stored integers may be 0.5 m "
            "steps (LOLA, read without its label) or whole meters. Pass the scale "
            "explicitly (meters per stored unit), or open a PDS product through its .LBL");
    }

    if (user_offset) {
        t.offset = *user_offset;
        t.offset_source = "command line";
        if (std::abs(t.offset) > REF_RADIUS_MIN_M)
            t.warnings.push_back("the offset " + fmt(t.offset) + " m from the command line is "
                                 "added to every height; a value this large is a radius, "
                                 "not a height");
    } else if (meta_offset) {
        const double off_m = *meta_offset * u();
        if (std::abs(off_m) > REF_RADIUS_MIN_M) {
            if (std::abs(off_m - body_radius_m) > REF_RADIUS_TOL_M)
                throw DemElevationError(
                    "the file's offset " + fmt(off_m) + " m is a reference radius, but not "
                    "the body radius " + fmt(body_radius_m) + " m. Pass the offset "
                    "explicitly (meters added to every height), or set the body radius");
            t.offset = off_m - body_radius_m;
            t.offset_source = "file metadata: reference radius " + fmt(off_m) + " m" + unit_note;
        } else if (std::abs(off_m * 1000.0 - body_radius_m) <= REF_RADIUS_TOL_M) {
            throw DemElevationError(
                "the file's offset " + fmt(off_m) + " equals the body radius in kilometers, "
                "so its heights are probably kilometers, but the file does not give a unit. "
                "Pass the scale (1000 for kilometers) and an offset of 0 explicitly");
        } else {
            t.offset = off_m;
            t.offset_source = "file metadata" + unit_note;
        }
    } else {
        t.offset = 0.0;
        t.offset_source = "default, no offset metadata";
    }
    return t;
}

}  // namespace dem
}  // namespace ptm
