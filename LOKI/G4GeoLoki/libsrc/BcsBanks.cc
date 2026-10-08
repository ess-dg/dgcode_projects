
#include "G4GeoLoki/BcsBanks.hh"
#include "Units/Units.hh"
#include <cmath>
#include <iostream>
#include <array>
#include <cassert>
#include <algorithm>

BcsBanks::BcsBanks(double rearBankDistance, int numberOfBanks, const std::string& bankCalibration)
  : m_rearBankDistance(rearBankDistance),
    m_numberOfBanks(numberOfBanks),
    m_bankCalibration(BankCalibration::load(bankCalibration))
{
  for (int bankId = 0; bankId <= 8; bankId++) {
    if (!m_bankCalibration.hasBank(bankId))
      continue;
    BankCalibration::Bank bank = m_bankCalibration.getBank(bankId);
    if (bankId == 0) // the rear bank moves along the beam: its distance is the rear detector distance (as nominal)
      bank.frontFaceCentre[2] = m_rearBankDistance;
    m_calibratedTransforms[bankId] = transformFromPlacement(bankId, bank);
  }
}

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

double BcsBanks::calcBankRotation(const int bankId){ //27.5+ (80-90)
  assert(0 <= bankId && bankId <= 8);
  return (90 - (bankTiltAngle[bankId] - bankPositionAngle[bankId])) *Units::deg;
}

const double BcsBanks::bankRotation[9][3] = { // all in mm
    {0.0,      0.5*M_PI, calcBankRotation(0)}, // 0 - rear
    {M_PI,     0.5*M_PI, calcBankRotation(1)},  // 1 - mid top
    {1.5*M_PI, 0.5*M_PI, calcBankRotation(2)}, // 2 - mid left
    {0.0,      0.5*M_PI, calcBankRotation(3)},  // 3 - mid bottom
    {0.5*M_PI, 0.5*M_PI, calcBankRotation(4)}, // 4 - mid right
    {M_PI,     0.5*M_PI, calcBankRotation(5)}, // 5 - front top
    {1.5*M_PI, 0.5*M_PI, calcBankRotation(6)}, // 6 - front left
    {0.0,      0.5*M_PI, calcBankRotation(7)}, // 7 - front bottom
    {0.5*M_PI, 0.5*M_PI, calcBankRotation(8)}, // 8 - front right
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

double BcsBanks::calcBankPositionZ(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  double intendedPosition = bankDistance[bankId] * std::cos(bankPositionAngle[bankId]*Units::deg);
  double bankCentreOffsetZ = detectorSystemCentreOffsetInBank(bankId, 2) * std::cos(calcBankRotation(bankId)) - detectorSystemCentreOffsetInBank(bankId, 1) * std::sin(calcBankRotation(bankId));
  return intendedPosition + bankCentreOffsetZ;
}
double BcsBanks::calcBankPositionXY(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  double intendedPosition = bankDistance[bankId] * std::sin(bankPositionAngle[bankId]*Units::deg);
  double bankCentreOffsetXY = detectorSystemCentreOffsetInBank(bankId, 2) * std::sin(calcBankRotation(bankId)) + detectorSystemCentreOffsetInBank(bankId, 1) * std::cos(calcBankRotation(bankId));

  return intendedPosition + bankCentreOffsetXY;
}

const int BcsBanks::bankPosDir[9] = { 1, 1, 1, -1, -1, 1, 1, -1, -1};

const double BcsBanks::bankPosition[9][3] = {
    {0, 0, 0}, // 0 - rear !calculated in getBankPosition function!
    {0, calcBankPositionXY(1), calcBankPositionZ(1)},  // 1 - mid top
    {calcBankPositionXY(2), 0, calcBankPositionZ(2)}, // 2 - mid left
    {0, calcBankPositionXY(3)*bankPosDir[3], calcBankPositionZ(3)},  // 3 - mid bottom
    {calcBankPositionXY(4)*bankPosDir[3], 0, calcBankPositionZ(4)}, // 4 - mid right
    {0, calcBankPositionXY(5), calcBankPositionZ(5)}, // 5 - front top
    {calcBankPositionXY(6), 0, calcBankPositionZ(6)}, // 6 - front left
    {0, calcBankPositionXY(7)*bankPosDir[3], calcBankPositionZ(7)}, // 7 - front bottom
    {calcBankPositionXY(8)*bankPosDir[3], 0, calcBankPositionZ(8)}, // 8 - front right
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

double BcsBanks::getBankRotation(const int bankId, const int axisIndex) {
  assert(0 <= bankId && bankId <= 8);
  assert(0 <= axisIndex && axisIndex <= 2);
  return bankRotation[bankId][axisIndex];
}

double BcsBanks::getBankPosition(const int bankId, const int axisIndex) const {
  assert(0 <= bankId && bankId <= 8);
  assert(0 <= axisIndex && axisIndex <= 2);
  if(bankId == 0){
    double rearBankPosition[] = {0, -detectorSystemCentreOffsetInBank(bankId, 1), this->m_rearBankDistance + detectorSystemCentreOffsetInBank(bankId, 2)};
    return rearBankPosition[axisIndex];
  }
  else{
    return bankPosition[bankId][axisIndex] + bankPositionOffset[bankId][axisIndex];
  }
}
std::array<double,3> BankTransform::toGlobal(const std::array<double,3>& local) const {
  std::array<double,3> global;
  for (int i = 0; i < 3; i++)
    global[i] = rotation[i][0] * local[0] + rotation[i][1] * local[1] + rotation[i][2] * local[2] + translation[i];
  return global;
}

std::array<double,3> BankTransform::toLocal(const std::array<double,3>& global) const {
  std::array<double,3> local;
  for (int j = 0; j < 3; j++)
    local[j] = rotation[0][j] * (global[0] - translation[0])
             + rotation[1][j] * (global[1] - translation[1])
             + rotation[2][j] * (global[2] - translation[2]);
  return local;
}

BankTransform BcsBanks::getBankTransform(const int bankId) const {
  assert(0 <= bankId && bankId <= 8);
  if (m_bankCalibration.hasBank(bankId))
    return m_calibratedTransforms[bankId];
  return getNominalBankTransform(bankId);
}

std::array<double,3> BcsBanks::getFrontFaceCentreInBank(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  // bank volume: local x = depth (towards the back), y = across the tubes, z = along the tubes
  return {-detectorSystemCentreOffsetInBank(bankId, 2), detectorSystemCentreOffsetInBank(bankId, 1), 0.0};
}

BankTransform BcsBanks::getNominalBankTransform(const int bankId) const {
  assert(0 <= bankId && bankId <= 8);
  return transformFromPlacement(bankId, nominalBankPlacement(bankId));
}

namespace {
  std::array<double,3> cross(const std::array<double,3>& a, const std::array<double,3>& b) {
    return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
  }
}

BankCalibration::Bank BcsBanks::nominalBankPlacement(const int bankId) const {
  assert(0 <= bankId && bankId <= 8);
  const double distance = bankId == 0 ? m_rearBankDistance : bankDistance[bankId] *Units::mm;
  const double positionAngle = bankPositionAngle[bankId] *Units::deg;
  const double normalAngle = calcBankRotation(bankId); // between the layer normal and the beam axis
  const auto& side = bankSideDirection[bankId];
  const std::array<double,3> beam = {0., 0., 1.};

  BankCalibration::Bank bank;
  std::array<double,3> towardsBeamAxis; // in the section plane, perpendicular to n: the local y axis of the bank
  for (int i = 0; i < 3; i++) {
    bank.frontFaceCentre[i] = distance * (std::sin(positionAngle) * side[i] + std::cos(positionAngle) * beam[i])
                              + bankPositionOffset[bankId][i] *Units::mm;
    bank.layerNormal[i] = std::sin(normalAngle) * side[i] + std::cos(normalAngle) * beam[i];
    towardsBeamAxis[i] = -std::cos(normalAngle) * side[i] + std::sin(normalAngle) * beam[i];
  }
  // the bank frame is [n, towardsBeamAxis, n x towardsBeamAxis] (see transformFromPlacement)
  const double s = areTubesInverselyNumbered(bankId) ? -1.0 : 1.0;
  const auto alongTubesAxis = cross(bank.layerNormal, towardsBeamAxis);
  for (int i = 0; i < 3; i++) {
    bank.acrossTubes[i] = s * towardsBeamAxis[i];
    bank.alongTubes[i] = -s * alongTubesAxis[i];
  }
  return bank;
}

BankTransform BcsBanks::transformFromPlacement(const int bankId, const BankCalibration::Bank& placement) {
  assert(0 <= bankId && bankId <= 8);
  const double s = areTubesInverselyNumbered(bankId) ? -1.0 : 1.0; // upside down banks: rotated by 180 deg about n
  BankTransform transform;
  for (int i = 0; i < 3; i++) {
    transform.rotation[i][0] = placement.layerNormal[i];
    transform.rotation[i][1] = s * placement.acrossTubes[i];
    transform.rotation[i][2] = -s * placement.alongTubes[i];
  }
  const auto frontFaceCentreInBank = getFrontFaceCentreInBank(bankId);
  for (int i = 0; i < 3; i++) {
    transform.translation[i] = placement.frontFaceCentre[i];
    for (int j = 0; j < 3; j++)
      transform.translation[i] -= transform.rotation[i][j] * frontFaceCentreInBank[j];
  }
  return transform;
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

/// bank notch ///
const double BcsBanks::nominalNotchCentre[3] = {-145.0 *Units::mm, -42.0 *Units::mm, -425.0 *Units::mm};
const double BcsBanks::nominalNotchHalfSize[3] = {15.0 *Units::mm, 70.0 *Units::mm, 60.0 *Units::mm};
const double BcsBanks::notchMargin = 1.0 *Units::mm;

bool BcsBanks::hasBankNotch(const int bankId) {
  return bankId == 6 || bankId == 8;
}

int BcsBanks::getBankNotchNeighbour(const int bankId) {
  assert(hasBankNotch(bankId));
  return bankId - 1;
}

namespace {
  struct Box { // oriented box: p_world = rotation * p_box + centre, |p_box[k]| <= halfSize[k]
    BankTransform transform;
    std::array<double,3> halfSize;
  };

  // The part of the segment p0-p1 inside the box (Liang-Barsky clipping); false if it is outside.
  bool clipSegment(const std::array<double,3>& p0, const std::array<double,3>& p1, const Box& box,
                   std::array<double,3>& q0, std::array<double,3>& q1) {
    const auto a = box.transform.toLocal(p0);
    const auto b = box.transform.toLocal(p1);
    double tMin = 0.0, tMax = 1.0;
    for (int k = 0; k < 3; k++) {
      const double d = b[k] - a[k];
      if (std::abs(d) < 1e-12) {
        if (std::abs(a[k]) > box.halfSize[k])
          return false;
        continue;
      }
      double t1 = (-box.halfSize[k] - a[k]) / d;
      double t2 = (box.halfSize[k] - a[k]) / d;
      if (t1 > t2)
        std::swap(t1, t2);
      tMin = std::max(tMin, t1);
      tMax = std::min(tMax, t2);
      if (tMin > tMax)
        return false;
    }
    for (int i = 0; i < 3; i++) {
      q0[i] = p0[i] + tMin * (p1[i] - p0[i]);
      q1[i] = p0[i] + tMax * (p1[i] - p0[i]);
    }
    return true;
  }

  // The parts of the 12 edges of box a that are inside box b
  void addEdgePointsInside(const Box& a, const Box& b, std::vector<std::array<double,3>>& points) {
    for (int corner = 0; corner < 8; corner++) {
      for (int axis = 0; axis < 3; axis++) {
        if (corner & (1 << axis))
          continue; // each edge once: from the corner with the lower coordinate along the axis
        std::array<double,3> c0, c1;
        for (int k = 0; k < 3; k++) {
          c0[k] = ((corner >> k) & 1) ? a.halfSize[k] : -a.halfSize[k];
          c1[k] = c0[k];
        }
        c1[axis] = a.halfSize[axis];
        std::array<double,3> q0, q1;
        if (clipSegment(a.transform.toGlobal(c0), a.transform.toGlobal(c1), b, q0, q1)) {
          points.push_back(q0);
          points.push_back(q1);
        }
      }
    }
  }
}

std::vector<std::array<double,3>> BcsBanks::getBankBoxOverlapPoints(const int bankId, const int otherBankId) const {
  // The overlap region of two boxes is a convex polyhedron; its vertices are on the edges of the two boxes, so
  // it is spanned by the parts of the edges of each box that are inside the other box.
  auto bankBox = [this](const int id) {
    return Box{getBankTransform(id), {0.5 * getBankSize(id, 2), 0.5 * getBankSize(id, 1), 0.5 * getBankSize(id, 0)}};
  };
  const Box a = bankBox(bankId);
  const Box b = bankBox(otherBankId);
  std::vector<std::array<double,3>> points;
  addEdgePointsInside(a, b, points);
  addEdgePointsInside(b, a, points);
  for (auto& point : points)
    point = a.transform.toLocal(point);
  return points;
}

BcsBanks::BankNotch BcsBanks::getBankNotch(const int bankId) const {
  assert(hasBankNotch(bankId));
  BankNotch notch;
  notch.isNominal = true;
  for (int k = 0; k < 3; k++) {
    notch.centre[k] = nominalNotchCentre[k];
    notch.halfSize[k] = nominalNotchHalfSize[k];
  }
  const int neighbour = getBankNotchNeighbour(bankId);
  if (!isBankCalibrated(bankId) && !isBankCalibrated(neighbour))
    return notch; // nominal placements
  const auto points = getBankBoxOverlapPoints(bankId, neighbour);
  bool insideNominalNotch = true;
  for (const auto& point : points)
    for (int k = 0; k < 3; k++)
      if (std::abs(point[k] - nominalNotchCentre[k]) > nominalNotchHalfSize[k])
        insideNominalNotch = false;
  if (insideNominalNotch)
    return notch; // (also if the boxes don't overlap)
  for (int k = 0; k < 3; k++) {
    double low = nominalNotchCentre[k] - nominalNotchHalfSize[k];
    double high = nominalNotchCentre[k] + nominalNotchHalfSize[k];
    for (const auto& point : points) {
      low = std::min(low, point[k] - notchMargin);
      high = std::max(high, point[k] + notchMargin);
    }
    notch.centre[k] = 0.5 * (low + high);
    notch.halfSize[k] = 0.5 * (high - low);
  }
  notch.isNominal = false;
  return notch;
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
