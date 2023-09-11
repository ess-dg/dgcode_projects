#include "G4GeoLoki/PixelatedBanks.hh"
#include "G4Units/Units.hh"
#include <cmath>
#include <iostream>
#include <array>
#include <cassert>

PixelatedBanks::PixelatedBanks(double rearBankDistance)
  : BcsBanks(rearBankDistance)
{
}
PixelatedBanks::PixelatedBanks(double rearBankDistance, int strawPixelNumber)
  : BcsBanks(rearBankDistance)
{
  for(int i=0; i<getNumberOfBanks(); i++) {
    numberOfPixelsInStraw[i] = strawPixelNumber;
  }
}
PixelatedBanks::PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks)
  : BcsBanks(rearBankDistance, numberOfBanks)
{
  for(int i=0; i<getNumberOfBanks(); i++) {
    numberOfPixelsInStraw[i] = strawPixelNumber;
  }
}

int PixelatedBanks::numberOfPixelsInStraw[9] = { // number of pixels along the straws
    256, // 0 - rear
    256,  // 1 - mid top
    256, // 2 - mid left
    256,  // 3 - mid bottom
    256, // 4 - mid right
    256, // 5 - front top
    256, // 6 - front left
    256, // 7 - front bottom
    256, // 8 - front right
};
int PixelatedBanks::getNumberOfPixelsInStraw(const int bankId) {
  assert(0 <= bankId && bankId <= 8);
  return numberOfPixelsInStraw[bankId];
}

int PixelatedBanks::getNumberOfPixels(const int bankId) {
  const int numberOfStrawsInBank = getNumberOfTubes(bankId) * 7;
  return numberOfStrawsInBank * getNumberOfPixelsInStraw(bankId);
}

int PixelatedBanks::getTotalNumberOfPixels() {
  return getBankPixelOffset(getNumberOfBanks());
}

int PixelatedBanks::getTubeLayerId(const int bankId, const int tubeId, const bool oldTubeNumbering) {
  const int tubePerLayer = getNumberOfTubes(bankId) / 4;
  return oldTubeNumbering ? (tubeId % 4) : (int) tubeId / tubePerLayer;
}


int PixelatedBanks::getBankPixelOffset(const int bankId) {
  assert(0 <= bankId && bankId <= 9);
  int offset = 0;
  for (int bankIndex = 0; bankIndex < bankId; bankIndex++) {
    const int numberOfStrawsInBank = getNumberOfTubes(bankIndex) * 7;
    offset += numberOfStrawsInBank * getNumberOfPixelsInStraw(bankIndex);
  }
  return offset;
}

int PixelatedBanks::getPositionPixelId(const int bankId, const double positionX, const double positionY) const{
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

int PixelatedBanks::getPixelId(const int bankId, const int tubeId, const int strawId, const double positionX, const double positionY) const{
  const int bankPixelOffset = getBankPixelOffset(bankId);
  const int strawPixelOffset = (tubeId * 7 + strawId) * getNumberOfPixelsInStraw(bankId);
  const int positionPixelId = getPositionPixelId(bankId, positionX, positionY);
  return bankPixelOffset + strawPixelOffset + positionPixelId;
}

void PixelatedBanks::dumpInfo(){
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