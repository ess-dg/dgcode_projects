#ifndef LOKI_MaskFileCreator_hh
#define LOKI_MaskFileCreator_hh

#include <iostream>
#include <memory>

class MaskFileCreator {
public:
  MaskFileCreator(const char* fileName, const int numberOfPixels, const int indexOffset);

  bool isPixelEntered(const int pixelNumber) const;
  void setPixelEntered(const int pixelNumber);
  bool isPixelEnteredAimingCheck(const int pixelNumber) const;
  void setPixelEnteredAimingCheck(const int pixelNumber);
  void createMaskFile() const;
private:
  const std::string m_fileName;
  const int m_numberOfPixels;
  const int m_indexOffset;
  std::unique_ptr<bool[]> m_enteredPixels;
  std::unique_ptr<bool[]> m_enteredPixelsAimingCheck;

  void checkAimingPixelCoverage() const;
};

#endif
