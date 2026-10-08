#include "G4GeoLoki/BcsBanks.hh"
#include "Units/Units.hh"
#include <cmath>
#include <array>
#include <cassert>

// The banks: the tables (also those of the nominal placement, used in BcsBanksPlacement.cc), the dimensions and
// content of the banks, the tube numbering, the masks and the beamstop (see BcsBanks.hh).

const double BcsBanks::packHolderDistanceFromPackTop = 7.6 *Units::mm;
const double BcsBanks::packHolderDistanceFromPackFront = 7.5 *Units::mm;

/// bank ///
const double BcsBanks::strawLengthInBank[9] = { // all in mm
    1000.0, // 0 - rear
    1000.0,  // 1 - mid top
    500.0, // 2 - mid left
    1000.0,  // 3 - mid bottom
    500.0, // 4 - mid right
    1200.0, // 5 - front top
    1200.0, // 6 - front left
    1200.0, // 7 - front bottom
    1200.0, // 8 - front right
};

const int BcsBanks::numberOfPacksInBank[9] = {
    28, // 0 - rear
    8,  // 1 - mid top
    6, // 2 - mid left
    8,  // 3 - mid bottom
    6, // 4 - mid right
    14, // 5 - front top
    16, // 6 - front left
    10, // 7 - front bottom
    16, // 8 - front right
};

const double BcsBanks::bankPositionAngle[9] = {
    0, // 0 - rear
    9.7,  // 1 - mid top
    6.5, // 2 - mid left
    9.7,  // 3 - mid bottom
    6.5, // 4 - mid right
    32.20, // 5 - front top
    25.6, // 6 - front left
    27.5, // 7 - front bottom
    25.6, // 8 - front right
};

const double BcsBanks::bankTiltAngle[9] = {
    90.0, // 0 - rear
    94.0,  // 1 - mid top
    92.0, // 2 - mid left
    94.0,  // 3 - mid bottom
    92.0, // 4 - mid right
    104.7, // 5 - front top
    105.0, // 6 - front left
    100.0, // 7 - front bottom //180-80=100
    105.0, // 8 - front right
};

const double BcsBanks::bankSideDirection[9][3] = { // unit vector from the beam axis towards the bank
    {0.0, -1.0, 0.0}, // 0 - rear (on the beam axis; its frame is oriented like the bottom banks)
    {0.0, 1.0, 0.0},  // 1 - mid top
    {1.0, 0.0, 0.0},  // 2 - mid left (x is positive to the left, looking along the beam)
    {0.0, -1.0, 0.0}, // 3 - mid bottom
    {-1.0, 0.0, 0.0}, // 4 - mid right
    {0.0, 1.0, 0.0},  // 5 - front top
    {1.0, 0.0, 0.0},  // 6 - front left
    {0.0, -1.0, 0.0}, // 7 - front bottom
    {-1.0, 0.0, 0.0}, // 8 - front right
};

const double BcsBanks::bankDistance[9] = {
    0.0, // 0 - rear //TODO ssd
    2950.0,  // 1 - mid top
    3350.0, // 2 - mid left
    2950.0,  // 3 - mid bottom
    3350.0, // 4 - mid right
    1364.68, // 5 - front top
    1750.0, // 6 - front left
    1340.0, // 7 - front bottom
    1750.0, // 8 - front right
};

const double BcsBanks::bankPositionOffset[9][3] = {
    {0.0, 0.0, 0.0}, // 0 - rear
    {0.0, 0.0, 0.0},  // 1 - mid top
    {0.0, 0.0, 0.0}, // 2 - mid left
    {0.0, 0.0, 0.0},  // 3 - mid bottom
    {0.0, 0.0, 0.0}, // 4 - mid right
    {100.0, 0.0, 0.0}, // 5 - front top
    {0.0, -15.0, 0.0}, // 6 - front left
    {-100.0, 0.0, 0.0}, // 7 - front bottom
    {0.0, 50.0, 0.0}, // 8 - front right
};
const double BcsBanks::bankSize[9][3] = { // x (width), y (height), z (depth)  //NOTE: x-y swap for vertical banks
    {1265.0, 1806.0, 285.00+15}, //0 - rear  {1265, 870+870+66, 285} //+20mm to avoid volume overlap of packs and the bank
    {1265.0, 590.18, 297.82},  // 1 - mid top {175+915+175, }
    {765.0, 447.45, 298.92}, // 2 - mid left
    {1265.0, 590.18, 297.82},  // 3 - mid bottom
    {765.0, 447.45, 298.92}, // 4 - mid right
    {1465.0, 945.00, 303.09}, // 5 - front top
    {1465.0, 1016.11, 302.99}, // 6 - front left
    {1465.0, 681.43, 302.99}, // 7 - front bottom
    {1465.0, 1016.11, 302.99}, // 8 - front right // changed from 1016.08
};

const double BcsBanks::topmostPackHolderPositionInBankFromTopFront[9][2] = { // y (height from bank top(*)), z (depth from bank front) (*)rotation?
    {117.26, 35.65}, //0 - rear
    {77.12, 31.33},  // 1 - mid top (upside down)
    {77.38, 31.33}, // 2 - mid left (upside down)
    {77.12, 31.33},  // 3 - mid bottom
    {77.38, 31.33}, // 4 - mid right
    {77.94, 31.33}, // 5 - front top (upside down)
    {77.92, 31.33}, // 6 - front left (upside down)
    {77.92, 31.33}, // 7 - front bottom
    {77.92, 31.33}, // 8 - front right
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 double BcsBanks::getPackRotation() {
   return BcsPack::getTubeGridParallelogramAngle();
 }

double BcsBanks::getPackPackDistance() {
  return BcsPack::getTubeGridParallelogramSide() * 2;
}


/// bank ///
double BcsBanks::getStrawLengthByBankId(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  return strawLengthInBank[bankId] * Units::mm;
}
int BcsBanks::getNumberOfPacksByBankId(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  return numberOfPacksInBank[bankId];
}
int BcsBanks::getNumberOfTubes(const int bankId){
  assert(0 <= bankId && bankId <= 8);
  return numberOfPacksInBank[bankId] * 8;
}

std::array<double,3> BcsBanks::getBankHalfSizeInBankFrame(const int bankId) {
  return {0.5 * getBankSize(bankId, 2), 0.5 * getBankSize(bankId, 1), 0.5 * getBankSize(bankId, 0)};
}

std::array<double,3> BcsBanks::getPackPositionInBankFrame(const int bankId, const int packNumber) {
  return {getPackPositionInBank(bankId, packNumber, 2), getPackPositionInBank(bankId, packNumber, 1),
          getPackPositionInBank(bankId, packNumber, 0)};
}

std::array<double,3> BcsBanks::getBoronMaskPositionInBankFrame(const int bankId, const int maskId) {
  return {getBoronMaskPosition(bankId, maskId, 2), getBoronMaskPosition(bankId, maskId, 1),
          getBoronMaskPosition(bankId, maskId, 0)};
}

std::array<double,3> BcsBanks::getCalibMaskPositionInBankFrame(CalibMasks::CalibMasksBase calibMask, const int bankId) const {
  return {getCalibMaskPosition(calibMask, bankId, 2), getCalibMaskPosition(calibMask, bankId, 1),
          getCalibMaskPosition(calibMask, bankId, 0)};
}

BcsBanks::PlacementInBank BcsBanks::getTriangularBoronMaskPlacementInBankFrame(const int maskId) {
  const int bankId = BoronMasks::getBankIdOfTriangularMask(maskId); // 5 or 7
  const double cutDir = BoronMasks::getCutDirOfTriangularMask(maskId, 1); // +1 (bank 5) or -1 (bank 7)
  PlacementInBank placement;
  // in front of the bank volume (touching it), at the given distances from the detector system centre: across the
  // tubes, and along the tubes (positive to the left looking along the beam, i.e. along cutDir * bank z)
  placement.position = {-0.5 * getBankSize(bankId, 2) - BoronMasks::getHalfSizeOfTriangularMask(maskId, 2),
                        detectorSystemCentreOffsetInBank(bankId, 1) - BoronMasks::getPosInBankOfTriangularMask(maskId, 1),
                        cutDir * BoronMasks::getPosInBankOfTriangularMask(maskId, 0)};
  // the mask axes: x along the tubes (to the left, looking along the beam), y across the tubes (upwards), z = depth
  // (the thickness); bank 5 (top) and bank 7 (bottom) have opposite bank y and z axes (see bankSideDirection)
  placement.rotation = {{{0.0, 0.0, 1.0}, {0.0, -cutDir, 0.0}, {cutDir, 0.0, 0.0}}}; // rows; columns = the mask axes
  return placement;
}

double BcsBanks::getBankSize(const int bankId, const int axisIndex) {
  assert(0 <= bankId && bankId <= 8);
  assert(0 <= axisIndex && axisIndex <= 2);
  return bankSize[bankId][axisIndex]+0.05 *Units::mm;
}

double BcsBanks::packHolderToPackCentreCoordsInPack(const int axisIndex) {
  assert(1 <= axisIndex && axisIndex <= 2);
  if(axisIndex == 1) { //y
    return 0.5*BcsPack::getPackBoxHeight() - packHolderDistanceFromPackTop;
  }
  else {//z
    return 0.5*BcsPack::getPackBoxWidth() - packHolderDistanceFromPackFront;
  }
}

double BcsBanks::getTopmostPackPositionInBank(const int bankId, const int axisIndex) {
  assert(0 <= bankId && bankId <= 8);
  assert(0 <= axisIndex && axisIndex <= 2);
  if(axisIndex == 0){ //x direction
    return 0;
  }
  else if(axisIndex == 1) { //y direction
  const double packHolderPositionYFromBankTop = topmostPackHolderPositionInBankFromTopFront[bankId][0];
  const double packHolderPositionY = 0.5* getBankSize(bankId, 1)  - packHolderPositionYFromBankTop;

  const double positionOfPackCentreFromPackHolderY = packHolderToPackCentreCoordsInPack(2) *std::sin(getPackRotation()) - packHolderToPackCentreCoordsInPack(1) *std::cos(getPackRotation());

  return packHolderPositionY + positionOfPackCentreFromPackHolderY;
  }
  else{ //z direction
  const double packHolderPositionZFromBankFront = topmostPackHolderPositionInBankFromTopFront[bankId][1];
  const double packHolderPositionZ = - (0.5*getBankSize(bankId, 2) - packHolderPositionZFromBankFront);

  const double positionOfPackCentreFromPackHolderZ = packHolderToPackCentreCoordsInPack(2) *std::cos(getPackRotation()) + packHolderToPackCentreCoordsInPack(1) *std::sin(getPackRotation());
  return packHolderPositionZ + positionOfPackCentreFromPackHolderZ;
  }
}

double BcsBanks::getPackPositionInBank(const int bankId, const int packNumber, const int axisIndex) {
  assert(0 <= axisIndex && axisIndex <= 2);
  if(axisIndex == 0 || axisIndex == 2) {
    return getTopmostPackPositionInBank(bankId, axisIndex);
  }
  else{
    return getTopmostPackPositionInBank(bankId, axisIndex) - packNumber * getPackPackDistance();
  }
}


double BcsBanks::packHolderToFirstTubeCentreCoordsInPack(const int axisIndex) {
  assert(1 <= axisIndex && axisIndex <= 2);
  if(axisIndex == 1) { //y
    return 0.5*BcsPack::getPackBoxHeight() + 0.5*BcsPack::getVerticalTubeDistanceInPack() - packHolderDistanceFromPackTop;
  }
  else {//z
    return BcsPack::getTubeCentreDistanceFromPackFront() - packHolderDistanceFromPackFront;
  }
}

double BcsBanks::detectorSystemFrontDistanceFromBankFront(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  const double packHolderPositionZFromBankFront = topmostPackHolderPositionInBankFromTopFront[bankId][1];
  const double positionOfFirstTubeCentreFromPackHolderZ = packHolderToFirstTubeCentreCoordsInPack(2) *std::cos(getPackRotation()) + packHolderToFirstTubeCentreCoordsInPack(1) *std::sin(getPackRotation());

  return packHolderPositionZFromBankFront + positionOfFirstTubeCentreFromPackHolderZ - BcsTube::getTubeOuterRadius();
}

double BcsBanks::detectorSystemCentreDistanceFromBankTop(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  const double packHolderDistanceYFromBankTop = topmostPackHolderPositionInBankFromTopFront[bankId][0];

  const double distanceOfSecondRowTubeCentreFromPackHolderY = std::abs(packHolderToFirstTubeCentreCoordsInPack(2) * std::sin(getPackRotation()) - packHolderToFirstTubeCentreCoordsInPack(1) * std::cos(getPackRotation()));
  const double distanceOfFirstRowTubeCentreFromPackHolderY = distanceOfSecondRowTubeCentreFromPackHolderY - BcsPack::getTubeGridParallelogramSide();
  const double detectorSystemSizeY = (numberOfPacksInBank[bankId] * 2 - 1) * BcsPack::getTubeGridParallelogramSide(); // from top row centre to lowest row centre

  const double distanceOfDetectorSystemCentreFromPackHolderY = distanceOfFirstRowTubeCentreFromPackHolderY + 0.5 * detectorSystemSizeY;

  return packHolderDistanceYFromBankTop + distanceOfDetectorSystemCentreFromPackHolderY;
}

double BcsBanks::detectorSystemCentreOffsetInBank(const int bankId, const int axisIndex) {
  assert(1 <= axisIndex && axisIndex <= 2);
  if (axisIndex == 1) {
    const double distanceFromBankTop = detectorSystemCentreDistanceFromBankTop(bankId);
    return 0.5 * getBankSize(bankId, 1) - distanceFromBankTop;
  }
  else { //axisIndex = 2
    const double distanceFromBankFront = detectorSystemFrontDistanceFromBankFront(bankId);
    return 0.5 * getBankSize(bankId, 2) - distanceFromBankFront;
  }
}

bool BcsBanks::isVertical(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  return (bankId == 2 || bankId == 4 || bankId == 6 || bankId == 8);
}
bool BcsBanks::areTubesInverselyNumbered(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  return (bankId == 1 || bankId == 2 || bankId == 5 || bankId == 6);
}

namespace {
  int getNumberOfTubeRows(const int bankId) {
    return 2 * BcsBanks::getNumberOfPacksByBankId(bankId);
  }
  // the tube row of a tube, counted from pack 0 (the first row of pack 0 is row 0)
  int getTubeRowFromPack0(const int bankId, const int tubeId) {
    const int row = tubeId % getNumberOfTubeRows(bankId); // in the numbering order
    return BcsBanks::areTubesInverselyNumbered(bankId) ? (getNumberOfTubeRows(bankId) - 1) - row : row;
  }
}

int BcsBanks::getTubeIdInBank(const int bankId, const int packId, const int inPackTubeId) {
  assert(0 <= packId && packId < getNumberOfPacksByBankId(bankId));
  assert(0 <= inPackTubeId && inPackTubeId <= 7);
  const int rowFromPack0 = 2 * packId + inPackTubeId / 4;
  const int row = areTubesInverselyNumbered(bankId) ? (getNumberOfTubeRows(bankId) - 1) - rowFromPack0 : rowFromPack0;
  const int layer = inPackTubeId % 4;
  return layer * getNumberOfTubeRows(bankId) + row;
}

int BcsBanks::getPackId(const int bankId, const int tubeId) {
  assert(0 <= tubeId && tubeId < getNumberOfTubes(bankId));
  return getTubeRowFromPack0(bankId, tubeId) / 2;
}

int BcsBanks::getInPackTubeId(const int bankId, const int tubeId) {
  assert(0 <= tubeId && tubeId < getNumberOfTubes(bankId));
  return (getTubeRowFromPack0(bankId, tubeId) % 2) * 4 + getTubeLayerId(bankId, tubeId);
}

int BcsBanks::getTubeLayerId(const int bankId, const int tubeId) {
  return tubeId / getNumberOfTubeRows(bankId);
}

int BcsBanks::getNumberOfBanks() const {
  return m_numberOfBanks;
}

/// Borom Masks ///

double BcsBanks::getBoronMaskPosition(const int bankId, const int maskId, const int axisIndex) {
  assert(0 <= axisIndex && axisIndex <= 2);
  const double thickness =  BoronMasks::getSize(bankId, maskId, axisIndex);
  const double position =  BoronMasks::getPosition(bankId, maskId, axisIndex);
  const double bankSizeHalf = 0.5*getBankSize(bankId, axisIndex);
  const double rotation =  BoronMasks::getRotation(bankId, maskId);

  if(axisIndex == 0){//x direction
    return bankSizeHalf - (position + 0.5*thickness);
  }
  else if(axisIndex == 1){  //y direction
    const double rotationCorrection = (1 - std::cos(rotation)) * 0.5*BoronMasks::getSize(bankId, maskId, 1) - std::sin(rotation) * 0.5* (-BoronMasks::getSize(bankId, maskId, 2));
    return bankSizeHalf - (position + 0.5*thickness) + rotationCorrection;
  }
  else{ //z direction
    const double rotationCorrection = (1 - std::cos(rotation)) * 0.5*(-BoronMasks::getSize(bankId, maskId, 2)) + std::sin(rotation) * 0.5*BoronMasks::getSize(bankId, maskId, 1);
    return - bankSizeHalf + position + 0.5*thickness + rotationCorrection;
  }
}

double BcsBanks::getCalibMaskPosition(CalibMasks::CalibMasksBase calibMask, const int bankId, const int axisIndex) const {
  assert(0 <= axisIndex && axisIndex <= 2);
  if((axisIndex == 0 && !isVertical(bankId)) || (axisIndex == 1 && isVertical(bankId))) {
    return 0.5*getStrawLengthByBankId(bankId) - calibMask.getLeftTubeEndDistance() - 0.5*calibMask.getWidth();
  }
  else if((axisIndex == 0 && isVertical(bankId)) || (axisIndex == 1 && !isVertical(bankId))) {
    return 0.0;
  }
  else{ //z direction
    return - (detectorSystemCentreOffsetInBank(bankId, 2) + calibMask.getElevationFromTubeFront() + 0.5*calibMask.getThickness());
  }
 }

const double BcsBanks::beamstopSize[5][3] = { // width, height, thickness [mm] (from ESS-1178830 Table 4.6 Selected beamstop sizes)
    {30., 35., 1.},   // id = 1 ("includes transmission detector")
    {20., 25., 1.},   // id = 2
    {50., 60., 1.},   // id = 3
    {65., 75., 1.},   // id = 4
    {100., 105., 1.}, // id = 5
};

double BcsBanks::getBeamstopSize(const int beamstopId, const int axisIndex) { // 0 - x, 1 - y, 2 - z
  assert(1 <= beamstopId && beamstopId <= 5);
  assert(0 <= axisIndex && axisIndex <= 2);
  return beamstopSize[beamstopId-1][axisIndex]*Units::mm;
}
// default beamstop could depending on collen: (smallest sufficient bs)
// collen = 3 m
//   rear det dist = 5 -> bs 5
// collen = 5 m
//   rear det dist = 5 -> bs 3
//   rear det dist = 10 -> bs 5
// collen = 8 m
//   rear det dist = 5 -> bs 3
//   rear det dist = 10 -> bs 4
