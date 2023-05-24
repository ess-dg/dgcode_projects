#ifndef G4GeoLoki_BcsBanks_hh
#define G4GeoLoki_BcsBanks_hh

#include <array>
#include "G4GeoLoki/BcsTube.hh"
#include "G4GeoLoki/BcsPack.hh"
#include "G4GeoLoki/BoronMasks.hh"
#include "G4GeoLoki/CalibMasks.hh"

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
};

#endif
