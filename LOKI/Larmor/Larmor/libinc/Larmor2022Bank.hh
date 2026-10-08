#ifndef Larmor_Larmor2022Bank_hh
#define Larmor_Larmor2022Bank_hh

#include <array>
#include "G4GeoLoki/BcsBanks.hh"

/// The LOKI rear bank (bank 0 of G4GeoLoki, 28 packs) in the test at the Larmor instrument (ISIS) in 2022.
///
/// A frozen copy of the bank placement, of the tube numbering (the current one and the old one of ICD v1, used by
/// the experiment) and of the pixel positions and ids, so that changes of the LOKI geometry code do not change the
/// Larmor results. The detector hardware (tubes, packs, the bank box) comes from G4GeoLoki.
///
/// Frames: the world frame of G4GeoLoki (x to the left looking along the beam, y up, z along the beam, the sample
/// at the origin) and the bank frame (x = depth, y = across the tubes, z = along the tubes, see BcsBanks.hh).
class Larmor2022Bank {
public:
  static constexpr int bankId = 0; // the LOKI rear bank in the tables of G4GeoLoki
  static constexpr int defaultNumberOfPixelsInStraw = 512;

  /// rearBankDistance: the distance of the detector front from the sample along the beam (4099 mm in the experiment)
  Larmor2022Bank(double rearBankDistance, int numberOfPixelsInStraw = defaultNumberOfPixelsInStraw);

  double getRearBankDistance() const { return m_rearBankDistance; }
  int getNumberOfPixelsInStraw() const { return m_numberOfPixelsInStraw; }
  int getTotalNumberOfPixels() const;
  static int getNumberOfPacks();
  static int getNumberOfTubes();

  /// The bank placement: the rear bank (oriented as in LOKI) at the beam height of the experiment.
  /// p_world = rotation * p_bank + translation.
  BankTransform getBankTransform() const;
  /// The y of the bank volume centre: the beam was 1155 mm above the floor.
  static double getBankPositionY();

  /// Tube numbering. tubeId: the copy number of the TubeWall volume. Current numbering (ICD v2): layer by layer
  /// (the front layer first), top to bottom in each layer; old numbering (ICD v1): pack by pack, top to bottom, the
  /// 4 layers of a tube row after each other.
  static int getTubeVolumeNumber(const int packNumber, const int inPackTubeId, const bool oldTubeNumbering);
  static int getPackId(const int tubeId, const bool oldTubeNumbering);
  static int getInPackTubeId(const int tubeId, const bool oldTubeNumbering);
  static int getTubeLayerId(const int tubeId, const bool oldTubeNumbering);

  /// Pixel ids: (tubeId * 7 + strawId) * pixels per straw + the pixel along the straw (from the left end, looking
  /// along the beam).
  int getTubeId(const int pixelId) const;
  int getStrawId(const int pixelId) const;
  std::array<double,3> getPixelCentre(const int pixelId, const bool oldTubeNumbering) const;
  /// The pixel id of a hit at the global position (x, y, z) in a straw.
  int getPixelId(const int tubeId, const int strawId, const double x, const double y, const double z) const;

  /// The Cd calibration slit mask in front of the bank (outside the bank volume): its description and the position of
  /// its centre in the world (it is placed with the frame rotation rotateY(90 deg), like the bank in 2020).
  static CalibMasks::CalibMasksBase getCalibMask();
  std::array<double,3> getCalibMaskPosition() const;

private:
  const double m_rearBankDistance;
  const int m_numberOfPixelsInStraw;
};

#endif
