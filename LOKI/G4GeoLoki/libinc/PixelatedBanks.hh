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

  /// The number of pixels per straw is a property of each object (constructor), 256 if not given.
  static constexpr int defaultNumberOfPixelsInStraw = 256;
  int getTotalNumberOfPixels() const;
  /// Pixel id from the global (x, y) of a hit, assuming the nominal bank placement (throws for calibrated banks).
  int getPixelId(const int bankId, const int tubeId, const int strawId, const double positionX, const double positionY) const;
  /// Pixel id from the global position of a hit (x, y, z): the pixel along the straw is taken from the
  /// position along the tubes in the bank frame (getBankTransform), so it also holds for rotated banks.
  int getPixelId(const int bankId, const int tubeId, const int strawId, const double positionX, const double positionY, const double positionZ) const;
  int getNumberOfPixels(const int bankId) const;
  int getNumberOfPixelsInStraw(const int bankId) const;

  static int getTubeLayerId(const int bankId, const int tubeId, const bool oldTubeNumbering);

  int getBankPixelOffset(const int bankId) const;
  void dumpInfo() const;

private:
  std::array<int,9> m_numberOfPixelsInStraw; // number of pixels along the straws, per bank
  int getPositionPixelId(const int bankId, const double positionX, const double positionY) const;
  int getLocalPositionPixelId(const int bankId, const double positionX, const double positionY, const double positionZ) const;
};

#endif
