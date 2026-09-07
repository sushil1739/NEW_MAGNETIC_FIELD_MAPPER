"""Utilities for building and inspecting TMS magnetic-field mapper files."""

from .core import (
    build_uniform_mapper,
    discover_plate_groups,
    read_fld,
    read_z_coordinates,
    reflect_axial_from_quarter,
)
from .edep import read_mapper_header, validate_mapper_file

__all__ = [
    "build_uniform_mapper",
    "discover_plate_groups",
    "read_fld",
    "read_z_coordinates",
    "reflect_axial_from_quarter",
    "read_mapper_header",
    "validate_mapper_file",
]
