#include "G4GeoLoki/BcsBanks.hh"
#include "Units/Units.hh"
#include <cmath>
#include <array>
#include <cassert>

// The placement of the banks in the world: BankTransform, the nominal placement from the tables and the placement
// from a bank calibration (see BcsBanks.hh). (The tables are in BcsBanks.cc.)

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

double BcsBanks::calcBankRotation(const int bankId){ //27.5+ (80-90)
  assert(0 <= bankId && bankId <= 8);
  return (90 - (bankTiltAngle[bankId] - bankPositionAngle[bankId])) *Units::deg;
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
