#ifndef G4GeoLoki_BcsBanks_hh
#define G4GeoLoki_BcsBanks_hh

#include <array>
#include <string>
#include <vector>
#include "G4GeoLoki/BcsTube.hh"
#include "G4GeoLoki/BcsPack.hh"
#include "G4GeoLoki/BoronMasks.hh"
#include "G4GeoLoki/CalibMasks.hh"
#include "G4GeoLoki/BankCalibration.hh"

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
  /// bankCalibration: the name of the bank calibration (see BankCalibration.hh), the default is the
  /// nominal geometry.
  BcsBanks(double rearBankDistance, int numberOfBanks = 9, const std::string& bankCalibration = BankCalibration::nominalName);
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
  /// Calibrated banks (see getBankCalibration): R = [n | s*w | -s*u] (s = -1 for the banks mounted
  /// upside down, see areTubesInverselyNumbered), translation = F - R * (front face centre in the bank).
  /// For the rear bank (0) the z of F is replaced by the rear detector distance.
  BankTransform getBankTransform(const int bankId, const bool isLarmor2022Experiment = false) const;
  /// The nominal placement of a bank (the same as getBankTransform without a calibration).
  BankTransform getNominalBankTransform(const int bankId, const bool isLarmor2022Experiment = false) const;
  /// The bank calibration in use (BankCalibration::nominalName for the nominal geometry).
  const BankCalibration& getBankCalibration() const { return m_bankCalibration; }
  bool isBankCalibrated(const int bankId) const { return m_bankCalibration.hasBank(bankId); }
  /// The centre of the front face of the detector system (front layer) in the bank frame.
  static std::array<double,3> getFrontFaceCentreInBank(const int bankId);
  static double getBankSize(const int bankId, const int axisIndex); // 0 - x, 1 - y, 2 - z

  /// The notch: an empty part of the bank box of the front left / right banks (6, 8) cut out where the
  /// bank box of the front top / bottom bank (5, 7) would overlap with it. Box in the bank frame (BankTransform
  /// local coordinates). The nominal notch is used if the overlap region of the two bank boxes (with the
  /// placements in use) is inside it; otherwise the notch is enlarged to cover the overlap region plus notchMargin.
  /// (Whether an enlarged notch stays clear of the packs is checked by the Geant4 overlap check of the packs.)
  struct BankNotch {
    std::array<double,3> centre;
    std::array<double,3> halfSize;
    bool isNominal;
  };
  static bool hasBankNotch(const int bankId);
  static int getBankNotchNeighbour(const int bankId); // 6 -> 5, 8 -> 7
  BankNotch getBankNotch(const int bankId) const;
  static const double notchMargin;
  /// The points spanning the overlap region of the bank boxes of two banks (empty if they don't overlap), in
  /// the frame of the first bank.
  std::vector<std::array<double,3>> getBankBoxOverlapPoints(const int bankId, const int otherBankId) const;

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
  BankCalibration m_bankCalibration;
  std::array<BankTransform,9> m_calibratedTransforms; // valid for the calibrated banks

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
  const static double nominalNotchCentre[3];
  const static double nominalNotchHalfSize[3];
};

#endif
