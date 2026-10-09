#include "G4GeoLoki/PixelatedBanks.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <array>
#include <cassert>
#include <stdexcept>

PixelatedBanks::PixelatedBanks(double rearBankDistance)
  : BcsBanks(rearBankDistance)
{
  setNumberOfPixelsInStraw(defaultNumberOfPixelsInStraw);
}

void PixelatedBanks::setNumberOfPixelsInStraw(const int strawPixelNumber) {
  if (strawPixelNumber < 1)
    throw std::invalid_argument("PixelatedBanks: the number of pixels per straw must be positive, not " + std::to_string(strawPixelNumber));
  m_numberOfPixelsInStraw.fill(strawPixelNumber);
}
PixelatedBanks::PixelatedBanks(double rearBankDistance, int strawPixelNumber)
  : BcsBanks(rearBankDistance)
{
  setNumberOfPixelsInStraw(strawPixelNumber);
}
PixelatedBanks::PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks)
  : BcsBanks(rearBankDistance, numberOfBanks)
{
  setNumberOfPixelsInStraw(strawPixelNumber);
}
PixelatedBanks::PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks, const std::string& bankCalibration)
  : BcsBanks(rearBankDistance, numberOfBanks, bankCalibration)
{
  setNumberOfPixelsInStraw(strawPixelNumber);
}
PixelatedBanks::PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks, const BankCalibration& bankCalibration)
  : BcsBanks(rearBankDistance, numberOfBanks, bankCalibration)
{
  setNumberOfPixelsInStraw(strawPixelNumber);
}

int PixelatedBanks::getNumberOfPixelsInStraw(const int bankId) const {
  checkBankId(bankId);
  return m_numberOfPixelsInStraw[bankId];
}

int PixelatedBanks::getNumberOfPixels(const int bankId) const {
  const int numberOfStrawsInBank = getNumberOfTubes(bankId) * 7;
  return numberOfStrawsInBank * getNumberOfPixelsInStraw(bankId);
}

int PixelatedBanks::getTotalNumberOfPixels() const {
  return getBankPixelOffset(getNumberOfBanks());
}


int PixelatedBanks::getBankPixelOffset(const int bankId) const {
  if (bankId < 0 || bankId > getNumberOfBanks())
    throw std::out_of_range("PixelatedBanks::getBankPixelOffset: bank " + std::to_string(bankId) + " (0-"
                            + std::to_string(getNumberOfBanks()) + ")");
  int offset = 0;
  for (int bankIndex = 0; bankIndex < bankId; bankIndex++) {
    const int numberOfStrawsInBank = getNumberOfTubes(bankIndex) * 7;
    offset += numberOfStrawsInBank * getNumberOfPixelsInStraw(bankIndex);
  }
  return offset;
}

int PixelatedBanks::getLocalPositionPixelId(const int bankId, const double positionX, const double positionY, const double positionZ) const{
  const double strawLength = getStrawLengthByBankId(bankId);
  const double pixelLength = strawLength / getNumberOfPixelsInStraw(bankId);
  // bank-local z is the straw axis, centred on the straw (see AimHelper::getPixelPositionInStraw)
  const double localZ = getBankTransform(bankId).toLocal({positionX, positionY, positionZ})[2];
  const double distanceFromFirstPixelEnd = areTubesInverselyNumbered(bankId) ? 0.5 * strawLength - localZ
                                                                             : localZ + 0.5 * strawLength;
  // a position at or beyond a straw end (e.g. exactly at the end) gives the end pixel of this straw, never a pixel
  // of the neighbouring straw (the counting gas ends inside the straw, so real hits are not affected)
  const int pixel = std::floor(distanceFromFirstPixelEnd / pixelLength);
  return std::min(std::max(pixel, 0), getNumberOfPixelsInStraw(bankId) - 1);
}

int PixelatedBanks::getPixelId(const int bankId, const int tubeId, const int strawId, const double positionX, const double positionY, const double positionZ) const{
  if (bankId < 0 || bankId >= getNumberOfBanks() || tubeId < 0 || tubeId >= getNumberOfTubes(bankId) || strawId < 0 || strawId > 6)
    throw std::out_of_range("PixelatedBanks::getPixelId: no bank " + std::to_string(bankId) + " tube " + std::to_string(tubeId)
                            + " straw " + std::to_string(strawId) + " in the geometry");
  const int bankPixelOffset = getBankPixelOffset(bankId);
  const int strawPixelOffset = (tubeId * 7 + strawId) * getNumberOfPixelsInStraw(bankId);
  return bankPixelOffset + strawPixelOffset + getLocalPositionPixelId(bankId, positionX, positionY, positionZ);
}

void PixelatedBanks::dumpInfo() const {
  int totalTumberOfPacks = 0;
  int totalTumberOfTubes = 0;
  int totalTumberOfStraws = 0;
  int totalTumberOfPixels = 0;
  for (int bankIndex = 0; bankIndex < 9; bankIndex++) {
    const int nPacksInBank = getNumberOfPacksByBankId(bankIndex);
    totalTumberOfPacks+=nPacksInBank;
    const int nTubesInBank = getNumberOfTubes(bankIndex);
    totalTumberOfTubes+=nTubesInBank;
    const int nStrawsInBank = nTubesInBank * 7;
    totalTumberOfStraws+=nStrawsInBank;
    const int pixelPerStraw = getNumberOfPixelsInStraw(bankIndex);
    const int nPixelsInBank = nStrawsInBank * pixelPerStraw;
    totalTumberOfPixels+=nPixelsInBank;

    auto indent = "    ";
    std::cout<<"Bank "<<bankIndex<<"\n";
    std::cout<<indent<<"Detector length: "<<getStrawLengthByBankId(bankIndex)<<" mm"<<"\n";
    std::cout<<indent<<"Pixels per straw: "<<pixelPerStraw<<"\n";
    std::cout<<indent<<"Number of packs: "<<nPacksInBank<<"\n";
    std::cout<<indent<<"Number of tubes: "<<nTubesInBank<<"\n";
    std::cout<<indent<<"Number of straws: "<<nStrawsInBank<<"\n";
    std::cout<<indent<<"Number of pixels: "<<nPixelsInBank<<" (starting from: "<<getBankPixelOffset(bankIndex)<<")"<<"\n";
  }
  std::cout<<"Total number of packs: "<<totalTumberOfPacks<<"\n";
  std::cout<<"Total number of tubes: "<<totalTumberOfTubes<<"\n";
  std::cout<<"Total number of straws: "<<totalTumberOfStraws<<"\n";
  std::cout<<"Total number of pixels: "<<totalTumberOfPixels<<"\n";
  
  std::cout<<"\nBeamstop options:\n";
  for (int beamstopId = 1; beamstopId <= 5; beamstopId++) {
    std::cout<<" id: "<< beamstopId; 
    std::cout<<", width: " <<getBeamstopSize(beamstopId, 0) << " mm";
    std::cout<<", height: " <<getBeamstopSize(beamstopId, 1) <<" mm " << "\n";
  }
}