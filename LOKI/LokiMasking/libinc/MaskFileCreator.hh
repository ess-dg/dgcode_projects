#ifndef LokiMasking_MaskFileCreator_hh
#define LokiMasking_MaskFileCreator_hh

#include <iostream>
#include <memory>
#include <vector>

class MaskFileCreator {
public:
  MaskFileCreator(const char* fileName, const int indexOffset, const std::vector<int>& bankPixelLimits, const int aimingBankId=0);

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
  std::vector<int> m_bankPixelLimits;
  const int m_aimingBankId;

  void checkAimingPixelCoverage(const int bankId) const;
};

#endif
