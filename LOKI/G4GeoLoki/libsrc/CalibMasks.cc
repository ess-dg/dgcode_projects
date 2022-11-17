#include "G4GeoLoki/CalibMasks.hh"
#include "Core/Units.hh"
#include <cmath>
#include <cassert>

// NOTE: All values are expressed in [mm] unit


//Planned calibration slit masks for LoKI, with 94 mm wide parts and 6 mm wide slits
const CalibMasks::CalibMasksBase mask0("lokiStandard-0", 0.3, 1600.0, -47.0, 
  {94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.});
const CalibMasks::CalibMasksBase mask1("lokiStandard-1", 0.3, 500.0, -47.0, 
  {94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.});
const CalibMasks::CalibMasksBase mask2("lokiStandard-2", 0.3, 500.0, -47.0, 
  {94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.});
const CalibMasks::CalibMasksBase mask3("lokiStandard-3", 0.3, 500.0, -47.0,
  {94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.});
const CalibMasks::CalibMasksBase mask4("lokiStandard-4", 0.3, 500.0, -47.0, 
  {94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.});
const CalibMasks::CalibMasksBase mask5("lokiStandard-5", 0.3, 820.0, -47.0, 
  {94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.});
const CalibMasks::CalibMasksBase mask6("lokiStandard-6", 0.3, 820.0, -47.0, 
  {94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.});
const CalibMasks::CalibMasksBase mask7("lokiStandard-7", 0.3, 600.0, -47.0, 
  {94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.});
const CalibMasks::CalibMasksBase mask8("lokiStandard-8", 0.3, 820.0, -47.0, 
  {94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.,6., 94.});

// B4C sheet(cadmium in real life) with holes(slits) cut into it, used for calibration at the LoKI rear bank experiment at Larmor(ISIS)
// From the right 76 mm (gap) – 74 mm (Cd) - 3 mm (slit) - 103 (Cd) - 3 (slit) - 103 (Cd) - 3 (slit) - 103 (Cd) - 3 (slit) - 103 (Cd) - 3 (slit) - 100 (Cd) - 3 (slit) - 100 (Cd) - 3 (slit) - 100 (Cd) - 3 (slit) - 100 (Cd) - 3 (slit) - 63 (Cd)
// '76 mm (gap)' means the distance from the right end of the detectors tubes, which translates to -50 mm distance from the left end of the BCS tubes 
CalibMasks::CalibMasksBase maskLarmor("larmorCdCalibMask", 0.3, 800., -50., (75.-59.0188),
  {63., 3.,100.,3.,100.,3.,100.,3.,100.,3., 103.,3.,103.,3.,103.,3.,103.,3., 74.});

const std::map<std::string, CalibMasks::CalibMasksBase> CalibMasks::m_masks { 
  {maskLarmor.getName(), maskLarmor}, 
  {mask0.getName(), mask0},
  {mask1.getName(), mask1},
  {mask2.getName(), mask2},
  {mask3.getName(), mask3},
  {mask4.getName(), mask4},
  {mask5.getName(), mask5},
  {mask6.getName(), mask5},
  {mask7.getName(), mask7},
  {mask8.getName(), mask8},
};


CalibMasks::CalibMasksBase::CalibMasksBase(std::string name, double thickness, double height, double leftTubeEndDistance, std::vector<double> pattern)
  : m_name(name), 
    m_thickness(thickness), 
    m_height(height), 
    m_leftTubeEndDistance(leftTubeEndDistance), 
    m_pattern(pattern)
{
}
CalibMasks::CalibMasksBase::CalibMasksBase(std::string name, double thickness, double height, double leftTubeEndDistance, double elevationFromMaskFront, std::vector<double> pattern)
  : m_name(name), 
    m_thickness(thickness), 
    m_height(height), 
    m_leftTubeEndDistance(leftTubeEndDistance), 
    m_elevationFromMaskFront(elevationFromMaskFront),
    m_pattern(pattern)
{
}

CalibMasks::CalibMasksBase CalibMasks::getCalibMask(std::string name){
  return m_masks.find(name)->second;
}

std::string CalibMasks::CalibMasksBase::getName() const{
  return this->m_name;
}
double CalibMasks::CalibMasksBase::getThickness() const{
  return this->m_thickness;
}
double CalibMasks::CalibMasksBase::getHeight() const{
  return this->m_height;
}
double CalibMasks::CalibMasksBase::getLeftTubeEndDistance() const{
  return this->m_leftTubeEndDistance;
}
double CalibMasks::CalibMasksBase::getElevationFromMaskFront() const{
  return this->m_elevationFromMaskFront;
}
std::vector<double> CalibMasks::CalibMasksBase::getPattern() const{
  return this->m_pattern;
}
double CalibMasks::CalibMasksBase::getWidth() const{
  return std::accumulate(this->m_pattern.begin(), this->m_pattern.end(), decltype(this->m_pattern)::value_type(0));
}

G4Material* CalibMasks::maskMaterial = NamedMaterialProvider::getMaterial("ESS_B4C:b10_enrichment=0.95");
G4Material* CalibMasks::maskBoxMaterial = NamedMaterialProvider::getMaterial("G4_Vacuum");