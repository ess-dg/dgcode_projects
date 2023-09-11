#include "G4GeoLoki/BoronMasks.hh"
#include "G4Units/Units.hh"
#include <cmath>
#include <iostream>
#include <array>
#include <cassert>

G4Material* BoronMasks::maskMaterial = NamedMaterialProvider::getMaterial("ESS_B4C:b10_enrichment=0.95");

/// regular rectangular masks ///

const std::array<std::array<double, 7>, 6> BoronMasks::rearBoronMasks = {{
    // dx, dy, dz, x,y,z, rot  (x,y,z of its top,left,front corner (in front view))
    {132.0, 2 * 870.0, 5.0, 0.0, 0.0, 0.0, 0.0},                 // left
    {1031.0, 101.0, 5.0, 117.5, 0.0, 59.0, 0.0},                 // top face
    {132.0, 2 * 870.0, 5.0, 1132.0, 0.0, 0.0, 0.0},              // right face
    {1031.0, 42.0, 5.0, 117.5, 870.0 + 870.0 - 42.0, 59.0, 0.0}, // bottom face
    {1216.0, 5.0, 125.0, 24.5, 95.39, 66.05, 13.45},             // top
    {1216.0, 5.0, 125.0, 24.5, 1698.06, 64.47, 13.45}            // bottom
}};

const std::array<std::array<double, 7>, 8> BoronMasks::midVerticalBoronMasks = {{
    // dx, dy, dz, x,y,z,  (x,y,z of its top,left,front corner (in front view))
    {150.0, 377.0, 5.0, 0.0, 70.45, 0.0, 0.0},                            // left face
    {150.0, 377.0, 5.0, 150.0 + 465.0, 70.45, 0.0, 0.0},                  // right face
    {465.0, 90.0, 5.0, 150.0, 447.45 - 90.0, 0.0, 0.0},                   // bottom face
    {700.0, 5.0, 215.0, 32.5, 416.49, 27.13, 13.45},                      // bottom
    {720.0, 47.0, 5.0, 22.5, 8.92, 59.24, 13.45},                         // ortogonal close to beam
    {135.5, 5.0, 300.0 - 0.05, 14.5, 70.03, 6.04, 13.45},                 // top left part
    {135.5, 5.0, 300.0 - 0.05, 14.5 + 135.5 + 465.0, 70.03, 6.04, 13.45}, // top right part
    {465.0, 5.0, 245.0 - 0.05, 14.5 + 135.5, 57.23, 59.53, 13.45},        // top middle part //57.24=70.03-cos(a)*55 59.53=6.04+sin(a)
}};

const std::array<std::array<double, 7>, 8> BoronMasks::midHorizontalBoronMasks = {{
    // dx, dy, dz, x,y,z,  (x,y,z of its top,left,front corner (in front view))
    {175.0, 519.5, 5.0, 0.0, 70.68, 0.0, 0.0},                          // left face
    {175.0, 519.5, 5.0, 175.0 + 915.0, 70.68, 0.0, 0.0},                // right face
    {915.0, 129.5, 5.0, 175.0, 590.18 - 129.5, 0.0, 0.0},               // top face
    {1200.0, 5.0, 215.0, 20.0 + 12.5, 532.53, 27.13, 13.45},            // bottom
    {1220.0, 49.0, 5.0, 22.5, 6.98, 57.68, 13.45},                      // ortogonal close to beam
    {160.5, 5.0, 300.0 - 0.05, 14.5, 70.03, 4.94, 13.45},               // top left part
    {160.5, 5.0, 300.0 - 0.05, 14.5 + 160.5 + 915, 70.03, 4.94, 13.45}, // top right part
    {915.0, 5.0, 245.0 - 0.05, 14.5 + 160.5, 57.23, 58.43, 13.45}       // top middle part //57.24=70.03-cos(a)*55 58.43=4.94+sin(a)
}};

const std::array<std::array<double, 7>, 8> BoronMasks::frontTopBoronMasks = {{
    // dx, dy, dz, x,y,z,  (x,y,z of its top,left,front corner (in front view))
    {210.0, 884.0, 5.0, 0.0, 61.0, 0.0, 0.0}, // left face
    //{37.0, 864.0-162.0, 5.0,  173.0, 71.0, 0.0,  0.0}, //left face addition
    {210.0, 884.0, 5.0, 210.0 + 1045.0, 61.0, 0.0, 0.0}, // right face
    //{37.0, 864.0-162.0, 5.0,  173.0+1119.0-37.0, 71.0, 0.0,  0.0}, //right face addition
    {1045.0, 162.0, 5.0, 210.0, 945.0 - 162.0 - 10.0, 0.0, 0.0},          // top face
    {1300.0, 5.0, 215.0, 82.5, 872.62, 22.27, 13.45},                     // bottom
    {1465.0, 43.0, 5.0, 0.0, 12.84, 64.25, 13.45},                        // ortogonal close to beam
    {173.0, 5.0, 305.0 - 0.05, 0.0, 71.2, 5.35 - 0.05, 13.45},            // top left part
    {173.0, 5.0, 305.0 - 0.05, 173.0 + 1119.0, 71.2, 5.35 - 0.05, 13.45}, // top right part
    {1119, 5.0, 250.0 - 0.05, 173, 58.40, 58.84 - 0.05, 13.45}            // top middle part //58.40=71.20-cos(a)*55 58.84=5.35+sin(a)
}};

const std::array<std::array<double, 7>, 8> BoronMasks::frontBottomBoronMasks = {{
    // dx, dy, dz, x,y,z,  (x,y,z of its top,left,front corner (in front view))
    {210.0, 620, 5.0, 0.0, 61.42, 0.0, 0.0}, // left face
    //{37.0, 600-125.5, 5.0,  173.0, 71.43, 0.0,  0.0}, //left face addition
    {210.0, 620, 5.0, 210.0 + 1045.0, 61.42, 0.0, 0.0}, // right face
    //{37.0, 600-125.5, 5.0,  173.0+1119.0-37.0, 71.43, 0.0,  0.0}, //right face addition
    {1045.0, 125.5, 5.0, 210.0, 681.43 - 125.5 - 10.0, 0.0, 0.0},  // bottom face
    {1300.0, 5, 215.0, 82.5, 645.51, 21.78, 13.45},                // bottom face
    {1465.0, 43.0, 5.0, 0.0, 12.84, 64.15, 13.45},                 // ortogonal close to beam
    {173.0, 5.0, 305.0 - 0.05, 0.0, 71.20, 5.25, 13.45},           // top: left part
    {1119.0, 5.0, 250.0 - 0.05, 173.0, 58.40, 58.74, 13.45},       // top: middle part //58.40=71.20-cos(a)*55 58.74=5.25+sin(a)*55
    {173.0, 5.0, 305.0 - 0.05, 173.0 + 1119.0, 71.20, 5.25, 13.45} // top: right
}};

const std::array<std::array<double, 7>, 8> BoronMasks::frontVerticalBoronMasks = {{
    // dx, dy, dz, x,y,z,  (x,y,z of its top,left,front corner (in front view))
    {173.0, 944.0, 5.0, 0.0, 1016.11 - 944.0, 0.0, 0.0},            // left
    {173.0, 944.0, 5.0, 173.0 + 1119.0, 1016.11 - 944.0, 0.0, 0.0}, // right face
    {1119.0, 129.5, 5.0, 173.0, 1016.11 - 129.5, 0.0, 0.0},         // bottom face
    {1300.0, 5.0, 215.0, 83.0, 986.31, 21.78, 13.45},               // bottom
    {1465.0, 48.0, 5.0, 0.0, 7.98, 62.98, 13.45},                   // ortogonal close to beam
    {173.0, 5.0, 305.0 - 0.05, 0.0, 71.20, 5.25, 13.45},            // top: left part
    {1119.0, 5.0, 250.0 - 0.05, 173.0, 58.40, 58.74, 13.45},        // top: middle part //58.40=71.20-cos(a)*55 58.74=5.25+sin(a)*55
    {173.0, 5.0, 305.0 - 0.05, 173.0 + 1119.0, 71.20, 5.25, 13.45}  // top: right
}};

double BoronMasks::getBoronMaskParameter(const int bankId, const int maskId, const int parameterIndex) {
  assert(0 <= bankId && bankId <= 8);
  assert(0 <= parameterIndex && parameterIndex <= 7);
  assert(0 <= maskId && maskId <= 8);
  switch (bankId) {
  case 0:
    assert(0 <= maskId && maskId <= 6);
    return rearBoronMasks[maskId][parameterIndex] *Units::mm;
  case 1:
  case 3:
    return midHorizontalBoronMasks[maskId][parameterIndex] *Units::mm;
  case 2:
  case 4:
    return midVerticalBoronMasks[maskId][parameterIndex] *Units::mm;
  case 5:
    return frontTopBoronMasks[maskId][parameterIndex] *Units::mm;
  case 6:
  case 8:
    return frontVerticalBoronMasks[maskId][parameterIndex] *Units::mm;
  case 7:
    return frontBottomBoronMasks[maskId][parameterIndex] *Units::mm;
  default:
    return 0;
  }
}

int BoronMasks::getNumberOfBoronMasks(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  switch (bankId) {
  case 0:
    return rearBoronMasks.size();
  case 1:
  case 3:
    return midHorizontalBoronMasks.size();
  case 2:
  case 4:
    return midVerticalBoronMasks.size();
  case 5:
    return frontTopBoronMasks.size();
  case 6:
  case 8:
    return frontVerticalBoronMasks.size();
  case 7:
    return frontBottomBoronMasks.size();
  default:
    return 0;
  }
}

double BoronMasks::getSize(const int bankId, const int maskId, const int axisIndex) {
  assert(0 <= axisIndex && axisIndex <= 2);
  return getBoronMaskParameter(bankId, maskId, axisIndex);
}
double BoronMasks::getPosition(const int bankId, const int maskId, const int axisIndex) {
  assert(0 <= axisIndex && axisIndex <= 2);
  return getBoronMaskParameter(bankId, maskId, axisIndex + 3);
}
double BoronMasks::getRotation(const int bankId, const int maskId) {
  return getBoronMaskParameter(bankId, maskId, 6) *Units::degree;
}

/// triangular masks ///

const std::array<std::array<double, 10>, 4> BoronMasks::triangularBoronMasks = {{
    // dxHalf dyHalf dzHalf, xCutPoint, yCutPoint, xCutDir(up/down), yCutDir(up/down), bankId, xCentreDistanceFromBankCentre, yCentreDistanceFromBankCentre
    {200, 215, 2.5, 64, 65, 1, 1, 5, 892.5, -137.53 + 16.5},  // SI-7615-691
    {280, 240, 2.5, 69, 52, -1, 1, 5, -972.5, -22.53 + 16.5}, // SI-7615-692
    {225, 190, 2.5, 74, 73, 1, -1, 7, 917.5, 119.51 + 16.5},  // SI-7615-693
    {185, 210, 2.5, 60, 0, -1, -1, 7, -877.5, 69.52 + 16.5}   // SI-7615-694
}};

double BoronMasks::getTriangularBoronMaskParameter(const int maskId, const int parameterIndex) {
  assert(0 <= maskId && maskId <= 3);
  assert(0 <= parameterIndex && parameterIndex <= 9);
  return triangularBoronMasks[maskId][parameterIndex];
}
double BoronMasks::getHalfSizeOfTriangularMask(const int maskId, const int axisIndex) {
  assert(0 <= axisIndex && axisIndex <= 2);
  return getTriangularBoronMaskParameter(maskId, axisIndex) *Units::mm;
}
double BoronMasks::getCutPointOfTriangularMask(const int maskId, const int axisIndex) {
  assert(0 <= axisIndex && axisIndex <= 1);
  return getTriangularBoronMaskParameter(maskId, 3 + axisIndex) *Units::mm;
}
double BoronMasks::getCutDirOfTriangularMask(const int maskId, const int axisIndex) {
  assert(0 <= axisIndex && axisIndex <= 1);
  return getTriangularBoronMaskParameter(maskId, 5 + axisIndex);
}
double BoronMasks::getBankIdOfTriangularMask(const int maskId) {
  return getTriangularBoronMaskParameter(maskId, 7);
}
double BoronMasks::getPosInBankOfTriangularMask(const int maskId, const int axisIndex) { //CentreDistanceFromBankCentre
  assert(0 <= axisIndex && axisIndex <= 1);
  return getTriangularBoronMaskParameter(maskId, 8 + axisIndex) *Units::mm;
}
