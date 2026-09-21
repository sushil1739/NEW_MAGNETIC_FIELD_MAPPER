# TMS Magnetic Field Mapper — Simulation Integration

## Integration flow

Field_maps
    |
    v
Mapper.txt
    |
    | field-local -> GDML-global translation
    v
Mapper_global.txt
    |
    | GDML ArbBField
    v
volTMS
    |
    v
edep-sim
    |
    v
Geant4 magnetic-field propagation
    |
    v
TMS simulation and reconstruction

## Validated coordinate transformation

Original local mapper header:

    -3800 -2500 -4000 100 100 10

TMS global origin:

    (0, -2202.23, 14860) mm

Field-to-TMS-local translation:

    (0, 850, 342.5) mm

Field-to-GDML-global translation:

    (0, -1352.23, 15202.5) mm

Resulting global mapper header:

    -3800 -3852.23 11202.5 100 100 10

Only coordinates are translated. The magnetic-field components and grid
spacing remain unchanged.

## Current validated ArbBField configuration

The diagnostic integration removes existing BField/ArbBField auxiliaries
from the six TMS steel logical volumes and from volTMS, then inserts one:

    <auxiliary
        auxtype="ArbBField"
        auxvalue="Mapper_global_DIAGNOSTIC.txt" />

on the parent logical volume:

    volTMS

The exact diagnostic injection script is stored in:

    simulation/inject_mapper_on_tms.py

## edep-sim validation evidence

edep-sim successfully reported:

    Reading Mapper_global_DIAGNOSTIC.txt ...
    Printing values for magnetic field.
    m_filename : Mapper_global_DIAGNOSTIC.txt
    m_offset   : -3800, -3852.23, 11202.5
    m_delta    : 100, 100, 10

Geant4 initialized:

    G4ClassicalRK4

with minimum step:

    0.01 mm

The GDML geometry also completed overlap validation with zero illegal
overlaps/extrusions.

Validated so far:

- GDML loads successfully
- ArbBField is recognized on volTMS
- mapper file opens successfully
- global mapper offset is read correctly
- 100 x 100 x 10 mm grid spacing is read correctly
- Geant4 magnetic-field integration initializes

## Diagnostic muon macro

The corrected smoke-test macro is:

    simulation/mapper_smoke.mac

Correct GPS syntax:

    /gps/position 500 -2202.23 11000 mm

The earlier form:

    /gps/position 500 mm -2202.23 mm 11000 mm

was invalid and caused Geant4 command error 503.

## Remaining physics validation

Field loading is validated.

The next required test is:

    field ON mu- trajectory
        versus
    field OFF mu- trajectory

This will verify the actual bending magnitude and direction produced by the
mapped field.

After that, the mapper-enabled GDML can move to larger edep-sim production
samples and TMS reconstruction.

## Generated files

Large/generated diagnostic files should not be committed:

    Mapper_global_DIAGNOSTIC.txt
    mapper_smoke.root
    mapper_smoke.log
    tms_mapper_MOTHER_DIAGNOSTIC.gdml
