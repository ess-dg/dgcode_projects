#include "G4GeoLoki/AimHelper.hh"

#include <cmath>
#include <iostream>
#include <array>
#include <cassert>

//////// Utilities for getting the centre coordinates of a pixel ////////
std::tuple<double,double,double> AimHelper::getPixelCentreCoordinates(const int pixelId) const {
  const int bankId = getBankId(pixelId);
  const int tubeId = getTubeId(pixelId, bankId);
  const int inPackTubeId = getInPackTubeId(bankId, tubeId);
  const int packId = getPackId(bankId, tubeId);
  const int strawId = getStrawId(pixelId, bankId, tubeId);

  ///////// pixel in straw /////////
  double positionZ = getPixelPositionInStraw(pixelId, bankId);

  ///////// straw in tube /////////
  double positionX = BcsTube::getStrawPositionX(strawId);
  double positionY = BcsTube::getStrawPositionY(strawId);

  ///////// tube in pack /////////
  // apply tube rotation
  coordinateRotation(positionX, positionY, BcsPack::getTubeRotation());
  // place tube in pack
  positionX += BcsPack::getHorizontalTubeOffset(inPackTubeId);
  positionY += BcsPack::getVerticalTubeOffset(inPackTubeId);

  ///////// pack in bank /////////
  // apply pack rotation
  coordinateRotation(positionX, positionY, getPackRotation());
  // place pack in bank
  const auto packPosition = getPackPositionInBankFrame(bankId, packId);
  positionX += packPosition[0];
  positionY += packPosition[1];
  positionZ += packPosition[2];

  ///////// bank in world /////////
  const auto global = getBankTransform(bankId).toGlobal({positionX, positionY, positionZ});

  return { global[0], global[1], global[2] };
}

void AimHelper::coordinateRotation(double &x, double &y, const double angle) {
  double tempX = std::cos(angle) * x - std::sin(angle) * y;
  double tempY = std::sin(angle) * x + std::cos(angle) * y;
  x = tempX;
  y = tempY;
}

int AimHelper::getBankId(const int pixelId) const {
  for (int bankId = 0; bankId < getNumberOfBanks(); bankId++){
    if(pixelId < getBankPixelOffset(bankId+1)){
      return bankId;
    }
  }
  throw std::runtime_error("Pixel id is out of the range for the banks in the geometry");
}

int AimHelper::getPackId(const int bankId, const int tubeId) {
  const int numberOfPacks = getNumberOfPacksByBankId(bankId);
  const int normalPackId = (int) (tubeId % (numberOfPacks * 2)) / 2;
  return !areTubesInverselyNumbered(bankId) ? normalPackId : ((numberOfPacks - 1) - normalPackId);
}

int AimHelper::getInPackTubeId(const int bankId, const int tubeId) {
  const int numberOfPacks = getNumberOfPacksByBankId(bankId);
  // the tube in the pack (0-3: the 4 layers of the first tube row of the pack, 4-7: of the second row) from the tube
  // id (layer by layer, 2 tube rows per pack, see GeoBCSBanks)
  const int inPackTubeId = ((tubeId % 2) * 4) + ((int) tubeId / (numberOfPacks * 2));
  return areTubesInverselyNumbered(bankId) ? (inPackTubeId + 4) % 8 : inPackTubeId % 8;
}

int AimHelper::getTubeId(const int pixelId, const int bankId) const {
  const int pixelIdInBank = pixelId - getBankPixelOffset(bankId);
  const int numberOfPixelsInATube = getNumberOfPixelsInStraw(bankId) * 7;
  return (int) pixelIdInBank / numberOfPixelsInATube;
}

int AimHelper::getStrawId(const int pixelId, const int bankId, const int tubeId) const {
  const int pixelIdInBank = pixelId - getBankPixelOffset(bankId);
  const int pixelIdInTube = pixelIdInBank - tubeId * 7 * getNumberOfPixelsInStraw(bankId);
  return (int) pixelIdInTube / getNumberOfPixelsInStraw(bankId);
}

double AimHelper::getPixelPositionInStraw(const int pixelId, const int bankId) const {
  const int locPixelId = pixelId % getNumberOfPixelsInStraw(bankId);
  const double pixelLength = getStrawLengthByBankId(bankId) / getNumberOfPixelsInStraw(bankId);
  const double position = - 0.5* getStrawLengthByBankId(bankId) + (locPixelId + 0.5) * pixelLength;
  return areTubesInverselyNumbered(bankId) ? - position : position;
}
