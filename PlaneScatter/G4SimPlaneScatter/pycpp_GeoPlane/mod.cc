#include "G4Interfaces/GeoConstructPyExport.hh"
#include "G4SimPlaneScatter/GeoPlane.hh"

PYTHON_MODULE
{
  GeoConstructPyExport::exportGeo<GeoPlane>(PYMOD "GeoPlane");
}
