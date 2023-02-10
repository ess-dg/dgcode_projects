#include "LOKI/MaskFileCreator.hh"
#include <fstream>
#include<algorithm>

MaskFileCreator::MaskFileCreator(const char* fileName, const int numberOfPixels, const int indexOffset):
  m_fileName(fileName),
  m_numberOfPixels(numberOfPixels),
  m_indexOffset(indexOffset),
  m_enteredPixels{new bool[numberOfPixels]},
  m_enteredPixelsAimingCheck{new bool[numberOfPixels]} {
    for(int i = 0; i < numberOfPixels; i++){
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
  std::ofstream maskFile;
  maskFile.open(m_fileName);
  maskFile << "<?xml version=\"1.0\"?>\n";
  maskFile << "<detector-masking>\n";
  maskFile << "\t<group>\n";
  maskFile << "\t\t<detids> ";

  for (int i = 0; i < m_numberOfPixels; i++) {
    if (m_enteredPixels[i] == false) {
      maskFile << i + m_indexOffset << ", ";
    }
  }
  maskFile.seekp(-2, std::ios_base::cur); //Go back with the write pointer to override the last coma and space ", "
  maskFile << " </detids>\n";
  maskFile << "\t</group>\n";
  maskFile << "</detector-masking>";
  maskFile.close();

  std::cout << "Created masking file: " << this->m_fileName << std::endl;
  checkAimingPixelCoverage();
}

void MaskFileCreator::checkAimingPixelCoverage() const {
  int nonEnteredPixels = 0;
  for (int i=0; i < m_numberOfPixels; i++) {
    if(!m_enteredPixelsAimingCheck[i]) {
      nonEnteredPixels++;
    }
  }
  if (nonEnteredPixels) {
    std::cout<< "WARNING: " << nonEnteredPixels << "/" << m_numberOfPixels << " pixels were not hit by any of the geantinos!\n";
  }
}
