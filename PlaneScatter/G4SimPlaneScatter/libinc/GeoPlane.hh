#ifndef G4SimPlaneScatter_GeoPlane_hh
#define G4SimPlaneScatter_GeoPlane_hh

#include "G4Interfaces/GeoConstructBase.hh"

//A simple plane parallel to x-y, running from z=0 to z=thickness.

class GeoPlane : public G4Interfaces::GeoConstructBase
{
public:
  GeoPlane();
  virtual ~GeoPlane(){}
  G4VPhysicalVolume* Construct();
};

#endif
