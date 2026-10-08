#include "Larmor/Larmor2022Bank.hh"
#include "Units/Units.hh"
#include <cassert>
#include <cmath>

Larmor2022Bank::Larmor2022Bank(double rearBankDistance, int numberOfPixelsInStraw)
  : m_rearBankDistance(rearBankDistance),
    m_numberOfPixelsInStraw(numberOfPixelsInStraw)
{
  assert(numberOfPixelsInStraw > 0);
}

int Larmor2022Bank::getNumberOfPacks() {
  return BcsBanks::getNumberOfPacksByBankId(bankId);
}

int Larmor2022Bank::getNumberOfTubes() {
  return getNumberOfPacks() * 8;
}

int Larmor2022Bank::getTotalNumberOfPixels() const {
  return getNumberOfTubes() * 7 * m_numberOfPixelsInStraw;
}

double Larmor2022Bank::getBankPositionY() {
  // beam centre at 1155 mm above the floor; the bank stood on a 4 mm electrical isolation layer, 33 mm above the
  // Larmor floor (due to the wheels)
  return 0.5 * BcsBanks::getBankSize(bankId, 1) - 1155 *Units::mm + (4+33) *Units::mm;
}

BankTransform Larmor2022Bank::getBankTransform() const {
  // the bank frame axes in the world: depth (x) along the beam, across the tubes (y) up, along the tubes (z) to
  // the right (-x)
  BankTransform transform;
  transform.rotation = {{{0.0, 0.0, -1.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}}};
  // the front face centre of the detector system on the beam axis at the rear bank distance (in z), the bank at
  // the beam height of the experiment (in y)
  const auto frontFaceCentreInBank = BcsBanks::getFrontFaceCentreInBank(bankId); // (depth, across, 0)
  transform.translation = {0.0, getBankPositionY(), m_rearBankDistance - frontFaceCentreInBank[0]};
  return transform;
}

int Larmor2022Bank::getTubeVolumeNumber(const int packNumber, const int inPackTubeId, const bool oldTubeNumbering) {
  assert(0 <= inPackTubeId && inPackTubeId <= 7);
  const double rowNumber = packNumber + 0.5 * ((int)inPackTubeId/4); // +0.5 for the second row (inPackTubeId > 3)
  const int layerNr = inPackTubeId % 4; // [0-3]
  if (oldTubeNumbering)
    return rowNumber * 8 + layerNr;
  return rowNumber * 2 + layerNr * getNumberOfPacks() * 2;
}

int Larmor2022Bank::getPackId(const int tubeId, const bool oldTubeNumbering) {
  return oldTubeNumbering ? (int) tubeId / 8 : (int) (tubeId % (getNumberOfPacks() * 2)) / 2;
}

int Larmor2022Bank::getInPackTubeId(const int tubeId, const bool oldTubeNumbering) {
  if (oldTubeNumbering)
    return tubeId % 8;
  const int newTubeIdConvertedToOldId = ((tubeId % 2) * 4) + ((int) tubeId / (getNumberOfPacks() * 2));
  return newTubeIdConvertedToOldId % 8;
}

int Larmor2022Bank::getTubeLayerId(const int tubeId, const bool oldTubeNumbering) {
  const int tubePerLayer = getNumberOfTubes() / 4;
  return oldTubeNumbering ? (tubeId % 4) : (int) tubeId / tubePerLayer;
}

int Larmor2022Bank::getTubeId(const int pixelId) const {
  return (int) pixelId / (m_numberOfPixelsInStraw * 7);
}

int Larmor2022Bank::getStrawId(const int pixelId) const {
  const int pixelIdInTube = pixelId - getTubeId(pixelId) * 7 * m_numberOfPixelsInStraw;
  return (int) pixelIdInTube / m_numberOfPixelsInStraw;
}

namespace {
  void coordinateRotation(double &x, double &y, const double angle) {
    double tempX = std::cos(angle) * x - std::sin(angle) * y;
    double tempY = std::sin(angle) * x + std::cos(angle) * y;
    x = tempX;
    y = tempY;
  }
}

std::array<double,3> Larmor2022Bank::getPixelCentre(const int pixelId, const bool oldTubeNumbering) const {
  assert(0 <= pixelId && pixelId < getTotalNumberOfPixels());
  const int tubeId = getTubeId(pixelId);
  const int inPackTubeId = getInPackTubeId(tubeId, oldTubeNumbering);
  const int packId = getPackId(tubeId, oldTubeNumbering);
  const int strawId = getStrawId(pixelId);

  // pixel in straw (the bank z axis is along the straws, centred on them)
  const int pixelIdInStraw = pixelId % m_numberOfPixelsInStraw;
  const double strawLength = BcsBanks::getStrawLengthByBankId(bankId);
  const double pixelLength = strawLength / m_numberOfPixelsInStraw;
  double positionZ = - 0.5 * strawLength + (pixelIdInStraw + 0.5) * pixelLength;

  // straw in tube
  double positionX = BcsTube::getStrawPositionX(strawId);
  double positionY = BcsTube::getStrawPositionY(strawId);

  // tube in pack
  coordinateRotation(positionX, positionY, BcsPack::getTubeRotation());
  positionX += BcsPack::getHorizontalTubeOffset(inPackTubeId);
  positionY += BcsPack::getVerticalTubeOffset(inPackTubeId);

  // pack in bank
  coordinateRotation(positionX, positionY, BcsBanks::getPackRotation());
  const auto packPosition = BcsBanks::getPackPositionInBankFrame(bankId, packId);
  positionX += packPosition[0];
  positionY += packPosition[1];
  positionZ += packPosition[2];

  // bank in world
  return getBankTransform().toGlobal({positionX, positionY, positionZ});
}

int Larmor2022Bank::getPixelId(const int tubeId, const int strawId, const double x, const double y, const double z) const {
  const double strawLength = BcsBanks::getStrawLengthByBankId(bankId);
  const double pixelLength = strawLength / m_numberOfPixelsInStraw;
  const double localZ = getBankTransform().toLocal({x, y, z})[2];
  const int pixelIdInStraw = std::floor((localZ + 0.5 * strawLength) / pixelLength);
  return (tubeId * 7 + strawId) * m_numberOfPixelsInStraw + pixelIdInStraw;
}

CalibMasks::CalibMasksBase Larmor2022Bank::getCalibMask() {
  // B4C sheet (cadmium in real life) with holes (slits) cut into it.
  // From the right 76 mm (gap) – 74 mm (Cd) - 3 mm (slit) - 103 (Cd) - 3 (slit) - 103 (Cd) - 3 (slit) - 103 (Cd) - 3 (slit) - 103 (Cd) - 3 (slit) - 100 (Cd) - 3 (slit) - 100 (Cd) - 3 (slit) - 100 (Cd) - 3 (slit) - 100 (Cd) - 3 (slit) - 63 (Cd)
  // '76 mm (gap)' means the distance from the right end of the detectors tubes, which translates to -50 mm distance from the left end of the BCS tubes
  // Placed 75 mm from the front of the detectors.
  return CalibMasks::CalibMasksBase("larmorCdCalibMask", 0.3, 800., -50., 75.,
    {63., 3.,100.,3.,100.,3.,100.,3.,100.,3., 103.,3.,103.,3.,103.,3.,103.,3., 74.});
}

std::array<double,3> Larmor2022Bank::getCalibMaskPosition() const {
  const auto calibMask = getCalibMask();
  const auto frontFaceCentreInBank = BcsBanks::getFrontFaceCentreInBank(bankId); // (depth, across, 0)
  const double offsetAlongTubes = 0.5 * BcsBanks::getStrawLengthByBankId(bankId) - calibMask.getLeftTubeEndDistance() - 0.5 * calibMask.getWidth();
  const double elevationFromBankCentre = -frontFaceCentreInBank[0] + calibMask.getElevationFromTubeFront() + 0.5 * calibMask.getThickness();
  // NOTE: the y is that of the LOKI rear bank centre (at 5 m), not the bank height of the experiment (getBankPositionY)
  return {offsetAlongTubes,
          -frontFaceCentreInBank[1],
          (m_rearBankDistance - frontFaceCentreInBank[0]) - elevationFromBankCentre};
}
