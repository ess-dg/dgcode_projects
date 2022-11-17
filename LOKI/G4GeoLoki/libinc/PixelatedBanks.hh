#ifndef G4GeoLoki_PixelatedBanks_hh
#define G4GeoLoki_PixelatedBanks_hh

#include "G4GeoLoki/BcsBanks.hh"

class PixelatedBanks : public BcsBanks{
public:
  PixelatedBanks(double rearBankDistance);
  PixelatedBanks(double rearBankDistance, int strawPixelNumber);
  PixelatedBanks(double rearBankDistance, int strawPixelNumber, int numberOfBanks);

  int getTotalNumberOfPixels();
  int getPixelId(const int bankId, const int tubeId, const int strawId, const double positionX, const double positionY) const;
  static int getNumberOfPixels(const int bankId);
  static int getNumberOfPixelsInStraw(const int bankId);

  static int getTubeLayerId(const int bankId, const int tubeId, const bool oldTubeNumbering);

protected:  
  static int getBankPixelOffset(const int bankId); 

private:
  static int numberOfPixelsInStraw[9];
  int getPositionPixelId(const int bankId, const double positionX, const double positionY) const;
};

#endif
