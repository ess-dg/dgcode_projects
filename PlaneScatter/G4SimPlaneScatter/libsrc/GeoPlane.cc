#include "G4SimPlaneScatter/GeoPlane.hh"
#include "G4Box.hh"

GeoPlane::GeoPlane()
  : GeoConstructBase("G4SimPlaneScatter/GeoPlane")
{
  addParameterDouble("thickness_mm",10.0);
  addParameterDouble("transverse_extent_m",1.0);
  addParameterString("material","Al");
}

G4VPhysicalVolume* GeoPlane::Construct()
{
  //Parameters:
  const double dz  = getParameterDouble("thickness_mm")*Units::mm;
  const double dxy = getParameterDouble("transverse_extent_m")*Units::m;
  auto mat = getParameterMaterial("material");
  auto mat_vacuum = getMaterial("Vacuum");

  //Note by TK: KK added double offset = 10*Units::cm here. Next time, just
  //introdue the offset as a parameter (default=0), and modify it in a copy of
  //the sim script. That way you won't change the behaviour when you hijack :-)
  const double offset = 0;

  const double dxy_world = dxy*1.001+offset;
  const double dz_world = 2*dz*1.001+offset;

  //World volume:
  auto worldvols = place(new G4Box("World", dxy_world, dxy_world, dz_world), mat_vacuum,0,0,0,0,INVISIBLE);
  auto lvWorld = worldvols.logvol;
  auto pvWorld = worldvols.physvol;

  //The plane (including detection planes):
  place(new G4Box("BackScatterPlane",   dxy,dxy,0.5*dz), mat_vacuum, 0, 0, -.5*dz + offset, lvWorld,G4Colour(1.0, 0.0, 0.0));
  place(new G4Box("Plane",              dxy,dxy,0.5*dz), mat,        0, 0, 0.5*dz + offset, lvWorld,G4Colour(0.0, 1.0, 0.0));
  place(new G4Box("ForwardScatterPlane",dxy,dxy,0.5*dz), mat_vacuum, 0, 0, 1.5*dz + offset, lvWorld,G4Colour(0.0, 0.0, 1.0));

  return pvWorld;
}
