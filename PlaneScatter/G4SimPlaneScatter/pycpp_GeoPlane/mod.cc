#include "G4Interfaces/GeoConstructPyExport.hh"
#include "G4SimPlaneScatter/GeoPlane.hh"

PYTHON_MODULE3
{
  GeoConstructPyExport::exportGeo<GeoPlane>(mod, "GeoPlane");
}
