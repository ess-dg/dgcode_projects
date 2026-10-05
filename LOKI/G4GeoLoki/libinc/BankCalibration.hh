#ifndef G4GeoLoki_BankCalibration_hh
#define G4GeoLoki_BankCalibration_hh

#include <array>
#include <string>

/// Named calibration of the LOKI detector bank placements (e.g. from a survey).
///
/// The name "nominal-geometry" means the built-in nominal bank placements (no file).
/// Any other name selects the data file G4GeoLoki/data/bank_calibration_<name>.txt; a value
/// containing a '/' is used as the path of a calibration file directly.
///
/// The file is written by the loki-geometry package (loki-calibrate --geant4-calibration-file).
/// Lines starting with '#' are comments. The first other line is the format line
/// "loki_bank_calibration 1", then "name <name>", then one line per calibrated bank:
///   bank <id> Fx Fy Fz ux uy uz wx wy wz nx ny nz
/// in the NeXus / Geant4 frame (x horizontal, positive to the left looking along the beam,
/// y up, z along the beam, origin at the sample position) in mm:
///   F: centre of the front face of the bank (the plane touching the front layer tubes,
///      in the middle of the front layer and of the straw length),
///   u: unit vector along the tubes, from the last pixel end to the first pixel end,
///   w: unit vector in the tube plane, perpendicular to the tubes, towards tube 0,
///   n: u x w, the layer normal, pointing away from the sample.
/// All banks (0-8) can be calibrated; banks without a line keep their nominal placement. The rear
/// bank (0) moves along the beam: its axes and the x and y of F are used, the z of F is replaced by
/// the rear detector distance (see BcsBanks).
class BankCalibration {
public:
  struct Bank {
    std::array<double,3> frontFaceCentre; // F [Geant4 length units]
    std::array<double,3> alongTubes;      // u
    std::array<double,3> acrossTubes;     // w
    std::array<double,3> layerNormal;     // n
  };

  static const std::string nominalName; // "nominal-geometry"
  /// The calibration used by default (the default of the bank_calibration geometry parameter and of the scripts).
  static const std::string defaultName;

  /// The nominal calibration (no calibrated bank).
  BankCalibration();
  /// Load a calibration by name (see above); throws std::runtime_error on any problem.
  static BankCalibration load(const std::string& name);

  const std::string& name() const { return m_name; }
  bool isNominal() const;
  bool hasBank(const int bankId) const;
  const Bank& getBank(const int bankId) const;
  /// The data file used ("" for the nominal geometry).
  const std::string& fileName() const { return m_fileName; }

private:
  std::string m_name;
  std::string m_fileName;
  std::array<bool,9> m_hasBank;
  std::array<Bank,9> m_banks;
  static BankCalibration loadFile(const std::string& fileName, const std::string& expectedName);
};

#endif
