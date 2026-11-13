// Read the metadata that dem_elevation.h's rule consumes from a GDAL band.
// Mirror of scripts/dem_elevation.py: gdal_band_metadata().
#pragma once

#include "dem_elevation.h"

#include <gdal_priv.h>

#include <optional>
#include <string>

namespace ptm {
namespace dem {

struct BandElevationMetadata {
    std::optional<double> scale;
    std::optional<double> offset;
    std::string unit;
    bool is_integer = false;
};

// Throws DemElevationError if the PDS label contradicts itself.
inline BandElevationMetadata readBandElevationMetadata(GDALDataset* ds, GDALRasterBand* band) {
    BandElevationMetadata m;
    int has = 0;
    double v = band->GetScale(&has);
    if (has) m.scale = v;
    has = 0;
    v = band->GetOffset(&has);
    if (has) m.offset = v;
    const char* ut = band->GetUnitType();
    m.unit = ut ? ut : "";
    char** pds = ds->GetMetadata("json:PDS");
    if (pds && pds[0]) {
        if (m.unit.empty()) m.unit = pdsUnitFromJson(pds[0]);
        if (!m.scale) m.scale = pdsNumberFromJson(pds[0], "SCALING_FACTOR");
        if (!m.offset) m.offset = pdsNumberFromJson(pds[0], "OFFSET");
    }
    m.is_integer = GDALDataTypeIsInteger(band->GetRasterDataType()) != 0;
    return m;
}

}  // namespace dem
}  // namespace ptm
