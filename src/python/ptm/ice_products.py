"""What an ice product in a run's ``statistics/`` directory means, verified.

``statistics/ice_depth.bin`` has held two different quantities: the depth of
the flux-threshold criterion (engines before the retention criterion), and the
depth of the ice table after the retention age (since). A reader must not
assume which. The engine and ``ptm_merge_chunks`` write, beside the retention
products, ``ice_retention.txt`` with one ``file_crc32: <name> <hex>`` line per
product file; this module establishes the criterion from that record and the
bytes on disk:

``retention``
    the record is present and its checksum matches the file: the file is the
    retention product the record describes.
``legacy_flux``
    ``ice_depth.bin`` with no record beside it: written by an engine before
    the retention criterion, so it is the flux-threshold depth. Only
    ``ice_depth.bin`` can be legacy; the other retention files never existed
    without a record.
``unverified``
    a record without checksums (written before they were added), a checksum
    that does not match (the file is not the one the record describes: a
    stale file from an earlier run), or a retention companion with no record.
    The criterion cannot be established; readers refuse or label it so.

Standard library only.
"""

from __future__ import annotations

import zlib
from dataclasses import dataclass
from pathlib import Path

RETENTION_RECORD = "ice_retention.txt"

# The files the retention criterion writes, all described by the record.
RETENTION_FILES = (
    "ice_depth.bin",
    "ice_table_layer.bin",
    "ice_retention_class.bin",
    "ice_toplayer_lifetime.bin",
    "ice_column_lifetime.bin",
    "ice_exposed_loss_tau.bin",
)

RETENTION, LEGACY_FLUX, UNVERIFIED = "retention", "legacy_flux", "unverified"


class IceProvenanceError(ValueError):
    """An ice product whose criterion cannot be established."""


@dataclass(frozen=True)
class IceProductProvenance:
    """``criterion`` is ``retention``, ``legacy_flux`` or ``unverified``;
    ``reason`` says why, in one sentence a user can act on."""
    name: str
    criterion: str
    reason: str

    @property
    def verified(self) -> bool:
        return self.criterion == RETENTION


def read_record_checksums(stats_dir) -> dict | None:
    """``{file name: crc32}`` from ``ice_retention.txt``; ``None`` when there
    is no record. A record without checksum lines gives ``{}``."""
    path = Path(stats_dir) / RETENTION_RECORD
    if not path.is_file():
        return None
    sums = {}
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("file_crc32:"):
            parts = line.split()
            if len(parts) == 3:
                try:
                    sums[parts[1]] = int(parts[2], 16)
                except ValueError:
                    pass
    return sums


def provenance(stats_dir, name: str = "ice_depth.bin") -> IceProductProvenance | None:
    """The criterion of ``statistics/<name>`` (a file of ``RETENTION_FILES``),
    or ``None`` when the file is absent."""
    stats_dir = Path(stats_dir)
    path = stats_dir / name
    if not path.is_file():
        return None
    sums = read_record_checksums(stats_dir)
    if sums is None:
        if name == "ice_depth.bin":
            return IceProductProvenance(
                name, LEGACY_FLUX,
                "no ice_retention.txt beside it: written by an engine before the retention "
                "criterion, so it is the flux-threshold depth")
        return IceProductProvenance(
            name, UNVERIFIED,
            "no ice_retention.txt beside it: a retention file without its record (left by an "
            "earlier run?)")
    if not sums:
        return IceProductProvenance(
            name, UNVERIFIED,
            "ice_retention.txt carries no file checksums (written before they were added), so "
            "it cannot be tied to this file; re-evaluate the run with ptm_merge_chunks --ice-only")
    if name not in sums:
        return IceProductProvenance(name, UNVERIFIED, "ice_retention.txt does not list this file")
    if zlib.crc32(path.read_bytes()) != sums[name]:
        return IceProductProvenance(
            name, UNVERIFIED,
            "the file does not match the checksum in ice_retention.txt: it is not the file that "
            "record describes (a stale file from an earlier run or criterion?)")
    return IceProductProvenance(name, RETENTION, "matches the checksum in ice_retention.txt")


def require_retention(stats_dir, name: str = "ice_depth.bin") -> IceProductProvenance:
    """``provenance`` for a reader that needs the retention product: raises
    ``IceProvenanceError`` unless it is verified (a legacy flux file, an
    unverifiable file, or no file)."""
    p = provenance(stats_dir, name)
    if p is None:
        raise IceProvenanceError(f"{Path(stats_dir) / name}: not present")
    if not p.verified:
        raise IceProvenanceError(f"{Path(stats_dir) / name}: {p.criterion}: {p.reason}")
    return p


def display_label(stats_dir, name: str, label: str, units: str = "") -> str:
    """The label a viewer shows for ``statistics/<name>``: ``label`` for a
    verified retention file or a file this module does not cover; for the
    legacy ``ice_depth.bin``, what it is (the flux-criterion depth); for an
    unverifiable file, ``label`` marked as such. ``units`` (for example
    ``" (m)"``) is appended after the name in every case."""
    if name not in RETENTION_FILES:
        return label + units
    p = provenance(stats_dir, name)
    if p is None or p.verified:
        return label + units
    if p.criterion == LEGACY_FLUX:
        return "Ice Depth, Flux Criterion, Legacy Run" + units
    return label + units + " [Criterion Unverified]"
