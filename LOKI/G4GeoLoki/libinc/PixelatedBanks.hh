#ifndef G4GeoLoki_PixelatedBanks_hh
#define G4GeoLoki_PixelatedBanks_hh

#include "G4GeoLoki/BcsBanks.hh"
#include <array>

class PixelatedBanks : public BcsBanks{
public:
  PixelatedBanks(double rearBankDistance);
  PixelatedBanks(double rearBankDistance, int strawPixelNumber);
  PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks);
  /// bankCalibration: see BcsBanks / BankCalibration.hh
  PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks, const std::string& bankCalibration);
  PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks, const BankCalibration& bankCalibration);

  /// The number of pixels per straw is a property of each object (constructor), 256 if not given.
  static constexpr int defaultNumberOfPixelsInStraw = 256;
  int getTotalNumberOfPixels() const;
  /// Pixel id from the global position of a hit (x, y, z): the pixel along the straw is taken from the
  /// position along the tubes in the bank frame (getBankTransform), so it also holds for rotated banks. A position
  /// at or beyond a straw end gives the end pixel of that straw. Throws std::out_of_range for a bank, tube or straw
  /// that is not in the geometry.
  int getPixelId(const int bankId, const int tubeId, const int strawId, const double positionX, const double positionY, const double positionZ) const;
  int getNumberOfPixels(const int bankId) const;
  int getNumberOfPixelsInStraw(const int bankId) const;

  int getBankPixelOffset(const int bankId) const;
  void dumpInfo() const;

private:
  std::array<int,9> m_numberOfPixelsInStraw; // number of pixels along the straws, per bank
  void setNumberOfPixelsInStraw(const int strawPixelNumber);
  int getLocalPositionPixelId(const int bankId, const double positionX, const double positionY, const double positionZ) const;
};

#endif
