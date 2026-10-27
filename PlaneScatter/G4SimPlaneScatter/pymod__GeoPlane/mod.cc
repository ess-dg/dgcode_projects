#include "G4Interfaces/GeoConstructToPy.hh"
#include "G4SimPlaneScatter/GeoPlane.hh"

SBLDPY_MODULE( _GeoPlane ) { G4Interfaces::geoToPy<GeoPlane>( m ); }
