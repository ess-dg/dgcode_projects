#ifndef G4GeoLoki_BcsTube_hh
#define G4GeoLoki_BcsTube_hh

#include "G4Materials/NamedMaterialProvider.hh"

class BcsTube {
public:
    /// materials ///
  static G4Material* tubeWallMaterial;
  static G4Material* tubeInnerGas;
  static G4Material* strawWallMaterial;
  static G4Material* converterMaterial;
  static G4Material* countingGas;

    /// straw ///
  static double getStrawOuterRadius();
  static double getStrawInnerRadius();
  static double getStrawWallThickness();

  static double getFrontTubeConverterThickness();
  static double getBackTubeConverterThickness();

    /// tube ///
  static double getTubeOuterRadius();
  static double getTubeInnerRadius();

  static double getStrawPositionX(const int strawId);
  static double getStrawPositionY(const int strawId);

private:
  const static double strawOuterRadius;
  const static double strawWallThickness;

  const static double frontTubeConverterThickness;
  const static double backTubeConverterThickness;

  const static double tubeOuterRadius;
  const static double tubeWallThickness;

  const static double strawPositionsInTube[7][2];
  const static double strawStrawDistance;
};

#endif
