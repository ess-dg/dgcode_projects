#ifndef G4GeoLoki_BankCalibration_hh
#define G4GeoLoki_BankCalibration_hh

#include <array>
#include <string>

/// Named calibration of the LOKI detector bank placements (e.g. from a survey).
///
/// The name "nominal-geant4-geometry" means the built-in nominal bank placements (no file).
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
/// The file must contain at least one bank line; |F| must be plausible (0.5-20 m, e.g. not in m instead of mm).
/// All banks (0-8) can be calibrated; banks without a line keep their nominal placement. The rear
/// bank (0) moves along the beam: its axes and the x and y of F are used, the z of F is replaced by
/// the rear detector distance (see BcsBanks).
///
/// Adding a new calibration (e.g. from a new survey), with the loki-geometry package
/// (https://github.com/MilanKlausz/loki-geometry, see its README "A new survey"):
///   1. add the survey file to loki-geometry (its entry in survey_files gives the calibration name),
///   2. loki-calibrate --survey-file NEW.txt --geant4-calibration-file <dgcode_projects>/LOKI/G4GeoLoki/data/
///      writes G4GeoLoki/data/bank_calibration_<name>.txt (all 9 banks; banks not in the survey from the
///      intended CAD geometry, marked with a comment),
///   3. rebuild (sb) to install the file; it is then selected with bank_calibration=<name>,
///   4. check the geometry for overlaps (add the name to the configurations of
///      G4GeoLokiTests/scripts/test_geo_dense_overlaps), and
///      verify the pixel JSON of loki-pixel-json -c <name> with sb_g4geoloki_verifyjson.
/// To make it the default, change defaultName (BankCalibration.cc); the tests that pin the default
/// (test_pixel_position_default and the SAM-606 sections of the tests) then need new logs.
/// Without a survey entry, any file of this format can be used by its path (bank_calibration=/path/to/file.txt).
class BankCalibration {
public:
  struct Bank {
    std::array<double,3> frontFaceCentre; // F [Geant4 length units]
    std::array<double,3> alongTubes;      // u
    std::array<double,3> acrossTubes;     // w
    std::array<double,3> layerNormal;     // n
  };

  static const std::string nominalName; // "nominal-geant4-geometry"
  /// The calibration used by default (the default of the bank_calibration geometry parameter and of the scripts).
  static const std::string defaultName;

  /// The nominal calibration (no calibrated bank).
  BankCalibration();
  /// Load a calibration by name (see above); throws std::runtime_error on any problem.
  static BankCalibration load(const std::string& name);
  /// A calibration from the text of its file (e.g. recorded in a simulation output, so that its analysis uses the
  /// placements of the simulation even if the file has changed or is not there): name is the bank_calibration value
  /// it was selected with; an empty text is the nominal geometry (only for the nominal name).
  static BankCalibration fromText(const std::string& name, const std::string& text);

  const std::string& name() const { return m_name; }
  bool isNominal() const;
  bool hasBank(const int bankId) const;
  const Bank& getBank(const int bankId) const;
  /// The data file used ("" for the nominal geometry, or a calibration from a text).
  const std::string& fileName() const { return m_fileName; }
  /// The text of the calibration file ("" for the nominal geometry).
  const std::string& text() const { return m_text; }

private:
  std::string m_name;
  std::string m_fileName;
  std::string m_text;
  std::array<bool,9> m_hasBank;
  std::array<Bank,9> m_banks;
  static BankCalibration loadFile(const std::string& fileName, const std::string& expectedName);
  /// Parse the text of a calibration file (source: the file name or another label for the error messages).
  static BankCalibration parse(const std::string& text, const std::string& source, const std::string& expectedName);
};

#endif
