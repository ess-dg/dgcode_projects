#ifndef G4GeoLoki_CalibMasks_hh
#define G4GeoLoki_CalibMasks_hh

#include <iostream>
#include <vector>
#include <map>
#include <numeric>
#include "G4Materials/NamedMaterialProvider.hh"

class CalibMasks {
public:
  class CalibMasksBase {
    public:
      CalibMasksBase(std::string name, double thickness, double height, double leftTubeEndDistance, std::vector<double> pattern);
      CalibMasksBase(std::string name, double thickness, double height, double leftTubeEndDistance, double elevationFromMaskFront, std::vector<double> pattern);

      std::string getName() const;
      double getThickness() const;
      double getHeight() const;
      double getLeftTubeEndDistance() const;
      double getElevationFromMaskFront() const;
      std::vector<double> getPattern() const;
      double getWidth() const;
    private:
      const std::string m_name;
      const double m_thickness;
      const double m_height;
      const double m_leftTubeEndDistance;
      const double m_elevationFromMaskFront = 0.0;
      const std::vector<double> m_pattern;
  };
  static CalibMasksBase getCalibMask(std::string name);
  static G4Material* maskMaterial;
  static G4Material* maskBoxMaterial;

private:
  const static std::map<std::string, CalibMasks::CalibMasksBase> m_masks;
};

#endif
