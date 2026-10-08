#include "G4GeoLoki/PixelatedBanks.hh"

#include <cmath>
#include <iostream>
#include <array>
#include <cassert>
#include <stdexcept>

PixelatedBanks::PixelatedBanks(double rearBankDistance)
  : BcsBanks(rearBankDistance)
{
  m_numberOfPixelsInStraw.fill(defaultNumberOfPixelsInStraw);
}
PixelatedBanks::PixelatedBanks(double rearBankDistance, int strawPixelNumber)
  : BcsBanks(rearBankDistance)
{
  m_numberOfPixelsInStraw.fill(defaultNumberOfPixelsInStraw);
  for(int i=0; i<getNumberOfBanks(); i++) {
    m_numberOfPixelsInStraw[i] = strawPixelNumber;
  }
}
PixelatedBanks::PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks)
  : BcsBanks(rearBankDistance, numberOfBanks)
{
  m_numberOfPixelsInStraw.fill(defaultNumberOfPixelsInStraw);
  for(int i=0; i<getNumberOfBanks(); i++) {
    m_numberOfPixelsInStraw[i] = strawPixelNumber;
  }
}
PixelatedBanks::PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks, const std::string& bankCalibration)
  : BcsBanks(rearBankDistance, numberOfBanks, bankCalibration)
{
  m_numberOfPixelsInStraw.fill(defaultNumberOfPixelsInStraw);
  for(int i=0; i<getNumberOfBanks(); i++) {
    m_numberOfPixelsInStraw[i] = strawPixelNumber;
  }
}

int PixelatedBanks::getNumberOfPixelsInStraw(const int bankId) const {
  assert(0 <= bankId && bankId <= 8);
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
  assert(0 <= bankId && bankId <= 9);
  int offset = 0;
  for (int bankIndex = 0; bankIndex < bankId; bankIndex++) {
    const int numberOfStrawsInBank = getNumberOfTubes(bankIndex) * 7;
    offset += numberOfStrawsInBank * getNumberOfPixelsInStraw(bankIndex);
  }
  return offset;
}

int PixelatedBanks::getPositionPixelId(const int bankId, const double positionX, const double positionY) const{
  if (isBankCalibrated(bankId)) // the (x, y) version assumes the nominal, axis-aligned bank placement
    throw std::logic_error("PixelatedBanks::getPixelId(bank, tube, straw, x, y) can not be used with a bank calibration, use getPixelId(bank, tube, straw, x, y, z)");
  const double pixelLength = getStrawLengthByBankId(bankId) / getNumberOfPixelsInStraw(bankId);

  if (isVertical(bankId)) { //vertical straw
    const double strawBegin = getBankPosition(bankId, 1) - 0.5* getStrawLengthByBankId(bankId);
    return std::floor((positionY - strawBegin) / pixelLength);
  }
  else { //horizontal straw
    const double strawBegin = getBankPosition(bankId, 0) - 0.5* getStrawLengthByBankId(bankId);
    const int invertedPixelId = std::floor((positionX - strawBegin) / pixelLength);
    return (getNumberOfPixelsInStraw(bankId) - 1) - invertedPixelId; //pixels are numbered in minus x direction
  }
}

int PixelatedBanks::getLocalPositionPixelId(const int bankId, const double positionX, const double positionY, const double positionZ) const{
  const double strawLength = getStrawLengthByBankId(bankId);
  const double pixelLength = strawLength / getNumberOfPixelsInStraw(bankId);
  // bank-local z is the straw axis, centred on the straw (see AimHelper::getPixelPositionInStraw)
  const double localZ = getBankTransform(bankId).toLocal({positionX, positionY, positionZ})[2];
  const double distanceFromFirstPixelEnd = areTubesInverselyNumbered(bankId) ? 0.5 * strawLength - localZ
                                                                             : localZ + 0.5 * strawLength;
  return std::floor(distanceFromFirstPixelEnd / pixelLength);
}

int PixelatedBanks::getPixelId(const int bankId, const int tubeId, const int strawId, const double positionX, const double positionY, const double positionZ) const{
  const int bankPixelOffset = getBankPixelOffset(bankId);
  const int strawPixelOffset = (tubeId * 7 + strawId) * getNumberOfPixelsInStraw(bankId);
  return bankPixelOffset + strawPixelOffset + getLocalPositionPixelId(bankId, positionX, positionY, positionZ);
}

int PixelatedBanks::getPixelId(const int bankId, const int tubeId, const int strawId, const double positionX, const double positionY) const{
  const int bankPixelOffset = getBankPixelOffset(bankId);
  const int strawPixelOffset = (tubeId * 7 + strawId) * getNumberOfPixelsInStraw(bankId);
  const int positionPixelId = getPositionPixelId(bankId, positionX, positionY);
  return bankPixelOffset + strawPixelOffset + positionPixelId;
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