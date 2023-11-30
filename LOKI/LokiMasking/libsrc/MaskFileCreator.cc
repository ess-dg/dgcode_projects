#include "LokiMasking/MaskFileCreator.hh"
#include <fstream>
#include <algorithm>

MaskFileCreator::MaskFileCreator(const char* fileName, const int indexOffset, const std::vector<int>& bankPixelLimits, const int aimingBankId):
  m_fileName(fileName),
  m_numberOfPixels(bankPixelLimits[bankPixelLimits.size()-1]),
  m_indexOffset(indexOffset),
  m_enteredPixels{new bool[m_numberOfPixels]},
  m_enteredPixelsAimingCheck{new bool[m_numberOfPixels]},
  m_bankPixelLimits(bankPixelLimits),
  m_aimingBankId(aimingBankId){
    for(int i = 0; i < m_numberOfPixels; i++){
        m_enteredPixels[i] = false;
        m_enteredPixelsAimingCheck[i] = false;
    }
}

bool MaskFileCreator::isPixelEntered(const int pixelNumber) const {
  return m_enteredPixels[pixelNumber];
}

void MaskFileCreator::setPixelEntered(const int pixelNumber) {
  if(pixelNumber < m_numberOfPixels) {
    m_enteredPixels[pixelNumber] = true;
  }
  //TODO handle error
}

bool MaskFileCreator::isPixelEnteredAimingCheck(const int pixelNumber) const {
  return m_enteredPixelsAimingCheck[pixelNumber];
}

void MaskFileCreator::setPixelEnteredAimingCheck(const int pixelNumber) {
  if(pixelNumber < m_numberOfPixels) {
    m_enteredPixelsAimingCheck[pixelNumber] = true;
  }
  //TODO handle error
}

void MaskFileCreator::createMaskFile() const {
  const int numberOfBanks = m_bankPixelLimits.size()-1;
  for (int bankId = 0; bankId<numberOfBanks; bankId++){
    if(m_aimingBankId>=0 && bankId != m_aimingBankId) {
      continue;
    }
    std::ofstream maskFile;
    std::string filename = m_fileName;
    if(numberOfBanks>1){
      filename = "bank" + std::to_string(bankId) + "_" + filename;
    }
    maskFile.open(filename);
    maskFile << "<?xml version=\"1.0\"?>\n";
    maskFile << "<detector-masking>\n";
    maskFile << "\t<group>\n";
    maskFile << "\t\t<detids> ";

    for (int i = m_bankPixelLimits[bankId]; i < m_bankPixelLimits[bankId+1]; i++) {
      if (m_enteredPixels[i] == false) {
        maskFile << i + m_indexOffset << ", ";
      }
    }
    maskFile.seekp(-2, std::ios_base::cur); //Go back with the write pointer to override the last coma and space ", "
    maskFile << " </detids>\n";
    maskFile << "\t</group>\n";
    maskFile << "</detector-masking>";
    maskFile.close();

    std::cout << "Created masking file: " << filename << std::endl;
    checkAimingPixelCoverage(bankId);
  }
}

void MaskFileCreator::checkAimingPixelCoverage(int bankId) const {
  int nonEnteredPixels = 0;
  for (int i = m_bankPixelLimits[bankId]; i < m_bankPixelLimits[bankId+1]; i++) {
    if(!m_enteredPixelsAimingCheck[i]) {
      nonEnteredPixels++;
    }
  }
  if (nonEnteredPixels) {
    const int pixelNumberInBank =  m_bankPixelLimits[bankId+1] - m_bankPixelLimits[bankId];
    std::cout<< "WARNING: " << nonEnteredPixels << "/" << pixelNumberInBank << " pixels were not hit by any of the geantinos!\n";
  }
}
