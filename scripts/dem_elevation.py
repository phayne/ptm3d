"""How PTM3D's DEM tools turn stored raster values into heights in meters.

One rule, shared by scripts/build_lola_ring.py and mirrored line for line in
src/cpp/tools/dem_elevation.h (dtm_to_ply_delaunay, dem_mesh_generator). Both
are checked against the same case table, scripts/dem_elevation_cases.tsv, so the
tools cannot disagree.

    height_m = stored_value * scale + offset

Scale (meters per stored unit):
  1. the command line, if given (meaning unchanged: meters per stored unit);
  2. else the file's scale metadata times its unit (a PDS label's
     SCALING_FACTOR and UNIT; GDAL's band scale);
  3. else 1 times the unit, for a floating-point raster;
  4. else STOP: an integer raster without scale metadata is ambiguous (raw LOLA
     DN are 0.5 m, other products store whole meters), so the user must say.

Offset (meters added to every height):
  1. the command line, if given (meaning unchanged; a warning if it is larger
     than 100 km, since that is a radius, not a height);
  2. else the file's offset metadata times its unit:
     - larger than 100 km: the radius of the DEM's reference sphere (PDS:
       PLANETARY_RADIUS = DN * SCALING_FACTOR + OFFSET), so the offset is
       (that radius - body radius). STOP if they differ by more than 20 km;
     - equal to the body radius in kilometers (within 20 km): the file stores
       kilometers without saying so. STOP and ask for an explicit scale/offset;
     - otherwise a height offset, applied as is;
  3. else 0.

Missing data: a stored value equal to the declared NoData value, NaN, or with
magnitude above 1e37 (float fill near FLT_MAX) is missing.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

REF_RADIUS_MIN_M = 100000.0   # an offset above this is a radius, not a height
REF_RADIUS_TOL_M = 20000.0    # how far a reference radius may sit from the body radius
FILL_ABS_MIN = 1e37           # float fill values near +-FLT_MAX

_UNITS = {"": 1.0, "M": 1.0, "METER": 1.0, "METERS": 1.0,
          "KM": 1000.0, "KILOMETER": 1000.0, "KILOMETERS": 1000.0}

_PDS_UNIT = re.compile(r'"UNIT"\s*:\s*"([^"]*)"')


class DemElevationError(ValueError):
    """The file's scale/offset cannot be interpreted without a user choice."""


@dataclass
class ElevationTransform:
    scale: float
    offset: float
    scale_source: str
    offset_source: str
    warnings: list = field(default_factory=list)

    def describe(self) -> str:
        return (f"height (m) = stored value * {self.scale:g} [{self.scale_source}] "
                f"+ {self.offset:g} [{self.offset_source}]")


def unit_factor(unit: str | None) -> float:
    key = (unit or "").strip().upper()
    if key not in _UNITS:
        raise DemElevationError(f"unknown elevation unit '{unit}' in the file metadata; "
                                f"pass the scale explicitly (meters per stored unit)")
    return _UNITS[key]


def pds_unit_from_json(text: str | None) -> str:
    """UNIT keyword of a PDS label, from GDAL's json:PDS metadata domain ('' if none)."""
    if not text:
        return ""
    m = _PDS_UNIT.search(text)
    return m.group(1) if m else ""


def resolve(meta_scale, meta_offset, meta_unit, is_integer, user_scale, user_offset,
            body_radius_m) -> ElevationTransform:
    """Apply the rule in the module docstring. None means 'absent'."""
    def u():
        return unit_factor(meta_unit)          # raises only where the unit is needed
    unit_note = f", unit {meta_unit.strip().upper()}" if (meta_unit or "").strip() else ""
    if meta_scale is not None and meta_scale == 0.0:
        meta_scale = None                      # a zero scale is not a scale
    if meta_offset is not None and meta_offset == 0.0:
        meta_offset = None
    warnings = []

    if user_scale is not None:
        scale, scale_src = float(user_scale), "command line"
    elif meta_scale is not None:
        scale, scale_src = float(meta_scale) * u(), f"file metadata{unit_note}"
    elif not is_integer:
        scale, scale_src = 1.0 * u(), f"default, no scale metadata{unit_note}"
    else:
        raise DemElevationError(
            "integer raster with no scale metadata: the stored integers may be 0.5 m "
            "steps (LOLA, read without its label) or whole meters. Pass the scale "
            "explicitly (meters per stored unit), or open a PDS product through its .LBL")

    if user_offset is not None:
        offset, offset_src = float(user_offset), "command line"
        if abs(offset) > REF_RADIUS_MIN_M:
            warnings.append(f"the offset {offset:g} m from the command line is added to every "
                            f"height; a value this large is a radius, not a height")
    elif meta_offset is not None:
        off_m = float(meta_offset) * u()
        if abs(off_m) > REF_RADIUS_MIN_M:
            if abs(off_m - body_radius_m) > REF_RADIUS_TOL_M:
                raise DemElevationError(
                    f"the file's offset {off_m:g} m is a reference radius, but not the body "
                    f"radius {body_radius_m:g} m. Pass the offset explicitly (meters added "
                    f"to every height), or set the body radius")
            offset = off_m - body_radius_m
            offset_src = f"file metadata: reference radius {off_m:g} m{unit_note}"
        elif abs(off_m * 1000.0 - body_radius_m) <= REF_RADIUS_TOL_M:
            raise DemElevationError(
                f"the file's offset {off_m:g} equals the body radius in kilometers, so its "
                f"heights are probably kilometers, but the file does not give a unit. Pass "
                f"the scale (1000 for kilometers) and an offset of 0 explicitly")
        else:
            offset, offset_src = off_m, f"file metadata{unit_note}"
    else:
        offset, offset_src = 0.0, "default, no offset metadata"

    return ElevationTransform(scale, offset, scale_src, offset_src, warnings)


def missing_mask(raw, nodata):
    """Boolean array: True where a stored value is missing (see module docstring)."""
    import numpy as np
    raw = np.asarray(raw, dtype=np.float64)
    miss = ~np.isfinite(raw) | (np.abs(raw) > FILL_ABS_MIN)
    if nodata is not None and math.isfinite(nodata):
        miss |= raw == float(nodata)
    return miss


def gdal_band_metadata(ds, band, gdal):
    """(scale, offset, unit, is_integer, nodata) of a GDAL band; None where absent."""
    unit = band.GetUnitType() or ""
    if not unit and "json:PDS" in (ds.GetMetadataDomainList() or []):
        unit = pds_unit_from_json((ds.GetMetadata_List("json:PDS") or [""])[0])
    int_types = {gdal.GDT_Byte, gdal.GDT_Int16, gdal.GDT_UInt16, gdal.GDT_Int32, gdal.GDT_UInt32}
    for name in ("GDT_Int8", "GDT_Int64", "GDT_UInt64"):
        if hasattr(gdal, name):
            int_types.add(getattr(gdal, name))
    return (band.GetScale(), band.GetOffset(), unit, band.DataType in int_types,
            band.GetNoDataValue())
