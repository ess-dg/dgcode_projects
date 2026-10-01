#ifndef G4GeoLoki_BcsBanks_hh
#define G4GeoLoki_BcsBanks_hh

#include <array>
#include "G4GeoLoki/BcsTube.hh"
#include "G4GeoLoki/BcsPack.hh"
#include "G4GeoLoki/BoronMasks.hh"
#include "G4GeoLoki/CalibMasks.hh"

/// Placement of a bank in the world: p_world = rotation * p_local + translation.
/// rotation[i][j] is the world component i of the bank-local axis j (the columns are
/// the images of the local x (depth), y (across the tubes) and z (along the tubes) axes);
/// translation is the centre of the bank box in the world [Geant4 length units].
struct BankTransform {
  std::array<std::array<double,3>,3> rotation;
  std::array<double,3> translation;
  std::array<double,3> toGlobal(const std::array<double,3>& local) const;
  std::array<double,3> toLocal(const std::array<double,3>& global) const;
};

class BcsBanks {
public:
  BcsBanks(double rearBankDistance, int numberOfBanks = 9);
  /// banks ///
  static double getPackRotation();
  static double getPackPackDistance();// TODO better name?
  static double getPackPositionInBank(const int bankId, const int packNumber, const int axisIndex);// 0 - x, 1 - y, 2 - z

  static double getStrawLengthByBankId(const int bankId);
  static int getNumberOfPacksByBankId(const int bankId);
  static int getNumberOfTubes(const int bankId);

  static double getBankRotation(const int bankId, const int axisIndex); // 0 - x, 1 - y, 2 - z
  double getBankPosition(const int bankId, const int axisIndex) const; // 0 - x, 1 - y, 2 - z
  /// The placement of a bank (single source for the Geant4 geometry, AimHelper and PixelatedBanks).
  /// Nominal: built from getBankRotation / getBankPosition (and the Larmor 2022 bank height).
  BankTransform getBankTransform(const int bankId, const bool isLarmor2022Experiment = false) const;
  static double getBankSize(const int bankId, const int axisIndex); // 0 - x, 1 - y, 2 - z

  static double detectorSystemFrontDistanceFromBankFront(const int bankId);

  static double getLarmor2022ExperimentBankPositionY();

  static bool isVertical(const int bankId);
  static bool areTubesInverselyNumbered(const int bankId);

  int getNumberOfBanks() const;
  /// boron masks ///
  static double getBoronMaskPosition(const int bankId, const int maskId, const int axisIndex);
  static double getTriangularBoronMaskPosition(const int maskId, const int axisIndex);

  /// calibration masks ///
  double getCalibMaskPosition(CalibMasks::CalibMasksBase calibMask,const int bankId, const int axisIndex) const;
  double getCalibMaskPositionOutsideBank(CalibMasks::CalibMasksBase calibMask,const int bankId, const int axisIndex) const;

  /// beamstop ///
  static double getBeamstopSize(const int beamstopId, const int axisIndex); // 0 - x, 1 - y, 2 - z

private:
  const double m_rearBankDistance;
  const int m_numberOfBanks;

  const static double packHolderDistanceFromPackTop;
  const static double packHolderDistanceFromPackFront;

  /// bank //
  const static double strawLengthInBank[9];
  const static int numberOfPacksInBank[9];

  const static double bankRotation[9][3];
  const static double bankPositionAngle[9];
  const static double bankTiltAngle[9];
  static double calcBankRotation(const int bankId);

  const static int bankPosDir[9]; //indicate direction along respective (X or Y) axis
  const static double bankPosition[9][3];
  const static double bankPositionOffset[9][3];
  const static double bankSize[9][3];
  const static double topmostPackHolderPositionInBankFromTopFront[9][2];
  static double getTopmostPackPositionInBank(const int bankId, const int axisIndex);
  static double packHolderToPackCentreCoordsInPack(const int axisIndex);

  const static double bankDistance[9];
  static double calcBankPositionZ(const int bankId);
  static double calcBankPositionXY(const int bankId);

  static double packHolderToFirstTubeCentreCoordsInPack(const int axisIndex);
  static double detectorSystemCentreDistanceFromBankTop(const int bankId);
  static double detectorSystemCentreOffsetInBank(const int bankId, const int axisIndex);

  const static double beamstopSize[5][3];
};

#endif
