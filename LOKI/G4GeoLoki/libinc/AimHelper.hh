#ifndef G4GeoLoki_AimHelper_hh
#define G4GeoLoki_AimHelper_hh

#include "G4GeoLoki/PixelatedBanks.hh"
#include <tuple>

class AimHelper : public PixelatedBanks {
public:
  using PixelatedBanks::PixelatedBanks; //inherit constructors of PixelatedBanks

  std::tuple<double,double,double> getPixelCentreCoordinates(const int pixelId) const;

  int getBankId(const int pixelId) const;
  static int getPackId(const int bankId, const int tubeId);
  int getTubeId(const int pixelId, const int bankId) const;
  int getStrawId(const int pixelId, const int bankId, const int tubeId) const;

private:
  static void coordinateRotation(double &x, double &y, const double angle);
  static int getInPackTubeId(const int bankId, const int tubeId);
  double getPixelPositionInStraw(const int pixelId, const int bankId) const;
};

#endif
