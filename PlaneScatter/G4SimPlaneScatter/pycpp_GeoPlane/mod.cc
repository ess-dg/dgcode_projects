#include "G4Interfaces/GeoConstructPyExport.hh"
#include "G4SimPlaneScatter/GeoPlane.hh"

PYTHON_MODULE( mod )
{
  GeoConstructPyExport::exportGeo<GeoPlane>(mod, "GeoPlane");
}
