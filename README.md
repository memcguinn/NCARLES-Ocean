# NCARLES-Ocean
National Center for Atmospheric Research Large Eddy Simulation (NCARLES) for Ocean Modelling

## RECENT EDITS
2019 - Skyler J. Kern; University of Colorado Boulder \
2021 - Mary E. McGuinn; University of Colorado Boulder

## ABOUT
NCARLES-Ocean is a pseudo-spectral solver for the wave averaged Boussinesq equations. The horizontal directions are solved using spectral techniques while the vertical direction uses finite difference. The code was original created to study the planetary boundary layer in the atmosphere by researchers at NCAR. The code has since been updated for ocean modeling.

NCARLES-Ocean is used to simulate a surface patch of the open ocean. The lateral boundaries are period while the top boundary has an imposed shear. The model is capable of tracking passive and active tracers. Tracer dynamics are in place to study carbonate chemistry.

## TO RUN
Refer to the "runscript_SAMPLE" compile and run script for Intel Fortran. Note that you will need to have OpenMPI installed.
