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

/// The LOKI detector banks: their dimensions and content (BcsBanks.cc), their placement in the world
/// (BcsBanksPlacement.cc) and the notch of the banks 6 and 8 (BcsBanksNotch.cc). The tubes and packs are described
/// by BcsTube and BcsPack, the masks by BoronMasks and CalibMasks.
class BcsBanks {
public:
  /// bankCalibration: the name of the bank calibration (see BankCalibration.hh), the default is the
  /// nominal geometry.
  BcsBanks(double rearBankDistance, int numberOfBanks = 9, const std::string& bankCalibration = BankCalibration::nominalName);
  /// With a bank calibration object (e.g. BankCalibration::fromText of the calibration recorded in a simulation output).
  BcsBanks(double rearBankDistance, int numberOfBanks, const BankCalibration& bankCalibration);
  int getNumberOfBanks() const;
  /// Throws std::out_of_range if bankId is not 0-8 (the functions of the bank placement and of the pixel ids check
  /// their input with it; the table accessors only assert, which release builds do not check).
  static void checkBankId(const int bankId);

  ///////// Banks: dimensions, packs, tubes (BcsBanks.cc) /////////
  static double getStrawLengthByBankId(const int bankId);
  static int getNumberOfPacksByBankId(const int bankId);
  static int getNumberOfTubes(const int bankId);
  static bool isVertical(const int bankId);
  static bool areTubesInverselyNumbered(const int bankId);
  static double getPackRotation();
  static double getPackPackDistance();// TODO better name?
  /// The tables of this class follow the drawings: axisIndex 0 = x = the width (along the tubes), 1 = y = the
  /// height (across the tubes), 2 = z = the depth (see the bank frame below).
  static double getBankSize(const int bankId, const int axisIndex); // 0 - x, 1 - y, 2 - z
  static double getPackPositionInBank(const int bankId, const int packNumber, const int axisIndex);// 0 - x, 1 - y, 2 - z
  static double detectorSystemFrontDistanceFromBankFront(const int bankId);

  /// The bank frame: the frame of a bank volume, of its contents and of BankTransform. x = depth (from the front,
  /// facing the sample, to the back), y = across the tubes, z = along the tubes (the axis of the tubes and straws,
  /// as of G4Tubs). The *InBankFrame functions give the values of the tables in this order.
  static std::array<double,3> getBankHalfSizeInBankFrame(const int bankId);
  static std::array<double,3> getPackPositionInBankFrame(const int bankId, const int packNumber);
  /// The centre of the front face of the detector system (front layer) in the bank frame.
  static std::array<double,3> getFrontFaceCentreInBank(const int bankId);

  /// Tube numbering: the tube id is the copy number of the TubeWall volume (and the tube part of the pixel id). A
  /// pack has 2 tube rows of 4 tubes, one tube in each layer (inPackTubeId 0-3: the first row, 4-7: the second row,
  /// layer = inPackTubeId % 4, see BcsPack); the tube rows of a bank are numbered from pack 0, the first row first
  /// (reversed for the banks mounted upside down, see areTubesInverselyNumbered). The tubes are numbered layer by
  /// layer (from the front), along the tube rows in each layer: tubeId = layer * (number of rows) + row.
  static int getTubeIdInBank(const int bankId, const int packId, const int inPackTubeId);
  /// The inverse of getTubeIdInBank: the pack, the tube in the pack and the layer (0: front ... 3: back) of a tube.
  static int getPackId(const int bankId, const int tubeId);
  static int getInPackTubeId(const int bankId, const int tubeId);
  static int getTubeLayerId(const int bankId, const int tubeId);

  ///////// Masks and beamstop (BcsBanks.cc) /////////
  static double getBoronMaskPosition(const int bankId, const int maskId, const int axisIndex);
  static std::array<double,3> getBoronMaskPositionInBankFrame(const int bankId, const int maskId);
  double getCalibMaskPosition(CalibMasks::CalibMasksBase calibMask,const int bankId, const int axisIndex) const;
  std::array<double,3> getCalibMaskPositionInBankFrame(CalibMasks::CalibMasksBase calibMask, const int bankId) const;
  /// A placement in the bank frame: position, and rotation (its columns: the axes of the placed volume).
  struct PlacementInBank {
    std::array<double,3> position;
    std::array<std::array<double,3>,3> rotation;
  };
  /// The triangular boron masks of banks 5 and 7 are mounted on the front of the bank (they are placed in the
  /// world, as they are outside the bank volume, with the transform of their bank).
  static PlacementInBank getTriangularBoronMaskPlacementInBankFrame(const int maskId);
  static double getBeamstopSize(const int beamstopId, const int axisIndex); // 0 - x, 1 - y, 2 - z

  ///////// Bank placement in the world (BcsBanksPlacement.cc) /////////
  /// The placement of a bank in the world (the single source for the Geant4 geometry, AimHelper and
  /// PixelatedBanks): from the front face centre F and the bank axes u, w, n (see BankCalibration.hh),
  /// R = [n | s*w | -s*u] (s = -1 for the banks mounted upside down, see areTubesInverselyNumbered) and
  /// translation = F - R * getFrontFaceCentreInBank. F, u, w, n come from the bank calibration in use, or, for
  /// the nominal geometry, from the tables (see nominalBankPlacement). For the rear bank (0) the z of F is the rear
  /// detector distance.
  BankTransform getBankTransform(const int bankId) const;
  /// The nominal placement of a bank (the same as getBankTransform without a calibration).
  BankTransform getNominalBankTransform(const int bankId) const;
  /// The nominal front face centre F and axes u, w, n of a bank, from the tables (in the terms of the drawing):
  /// the bank is on one side of the beam (bankSideDirection), in the section plane spanned by the beam axis and that
  /// side; F is at the bank distance along the position angle in that plane, plus the panel offset; the layer
  /// normal n is tilted from the beam axis by calcBankRotation (= 90 deg - (face angle - position angle)) towards
  /// the side.
  BankCalibration::Bank nominalBankPlacement(const int bankId) const;
  /// The bank calibration in use (BankCalibration::nominalName for the nominal geometry).
  const BankCalibration& getBankCalibration() const { return m_bankCalibration; }
  bool isBankCalibrated(const int bankId) const { return m_bankCalibration.hasBank(bankId); }

  ///////// The notch of the banks 6 and 8 (BcsBanksNotch.cc) /////////
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

private:
  const double m_rearBankDistance;
  const int m_numberOfBanks;
  BankCalibration m_bankCalibration;
  std::array<BankTransform,9> m_calibratedTransforms; // valid for the calibrated banks

  ///////// Tables and helpers of the banks (BcsBanks.cc) /////////
  const static double packHolderDistanceFromPackTop;
  const static double packHolderDistanceFromPackFront;
  const static double strawLengthInBank[9];
  const static int numberOfPacksInBank[9];
  const static double bankSize[9][3];
  const static double topmostPackHolderPositionInBankFromTopFront[9][2];
  static double getTopmostPackPositionInBank(const int bankId, const int axisIndex);
  static double packHolderToPackCentreCoordsInPack(const int axisIndex);
  static double packHolderToFirstTubeCentreCoordsInPack(const int axisIndex);
  static double detectorSystemCentreDistanceFromBankTop(const int bankId);
  static double detectorSystemCentreOffsetInBank(const int bankId, const int axisIndex);
  const static double beamstopSize[5][3];

  ///////// Tables of the nominal placement (BcsBanks.cc) and placement helpers (BcsBanksPlacement.cc) /////////
  const static double bankSideDirection[9][3]; // from the beam axis towards the bank (in the section plane)
  const static double bankDistance[9];
  const static double bankPositionAngle[9];
  const static double bankTiltAngle[9];
  const static double bankPositionOffset[9][3];
  static double calcBankRotation(const int bankId);
  /// The transform of a bank placed by its front face centre and axes (see getBankTransform).
  static BankTransform transformFromPlacement(const int bankId, const BankCalibration::Bank& placement);

  ///////// The nominal notch (BcsBanksNotch.cc) /////////
  const static double nominalNotchCentre[3];
  const static double nominalNotchHalfSize[3];
};

#endif
