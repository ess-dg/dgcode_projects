#ifndef G4GeoLoki_PixelatedBanks_hh
#define G4GeoLoki_PixelatedBanks_hh

#include "G4GeoLoki/BcsBanks.hh"

class PixelatedBanks : public BcsBanks{
public:
  PixelatedBanks(double rearBankDistance);
  PixelatedBanks(double rearBankDistance, int strawPixelNumber);
  PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks);
  /// bankCalibration: see BcsBanks / BankCalibration.hh
  PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks, const std::string& bankCalibration);

  int getTotalNumberOfPixels();
  /// Pixel id from the global (x, y) of a hit, assuming the nominal bank placement (throws for calibrated banks).
  int getPixelId(const int bankId, const int tubeId, const int strawId, const double positionX, const double positionY) const;
  /// Pixel id from the global position of a hit (x, y, z): the pixel along the straw is taken from the
  /// position along the tubes in the bank frame (getBankTransform), so it also holds for rotated banks.
  int getPixelId(const int bankId, const int tubeId, const int strawId, const double positionX, const double positionY, const double positionZ) const;
  static int getNumberOfPixels(const int bankId);
  static int getNumberOfPixelsInStraw(const int bankId);

  static int getTubeLayerId(const int bankId, const int tubeId, const bool oldTubeNumbering);

  static int getBankPixelOffset(const int bankId);
  static void dumpInfo();

private:
  static int numberOfPixelsInStraw[9];
  int getPositionPixelId(const int bankId, const double positionX, const double positionY) const;
  int getLocalPositionPixelId(const int bankId, const double positionX, const double positionY, const double positionZ) const;
};

#endif
