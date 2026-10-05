#include "G4GeoLoki/BankCalibration.hh"
#include "Core/FindData.hh"
#include "Units/Units.hh"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

const std::string BankCalibration::nominalName = "nominal-geometry";
const std::string BankCalibration::defaultName = "2026-September-SAM-606";

std::string BankCalibration::effectiveName(const std::string& name, const bool isLarmor2022Experiment) {
  return isLarmor2022Experiment ? nominalName : name;
}

BankCalibration::BankCalibration()
  : m_name(nominalName), m_fileName("")
{
  m_hasBank.fill(false);
}

bool BankCalibration::isNominal() const {
  for (bool has : m_hasBank)
    if (has)
      return false;
  return true;
}

bool BankCalibration::hasBank(const int bankId) const {
  return 0 <= bankId && bankId <= 8 && m_hasBank[bankId];
}

const BankCalibration::Bank& BankCalibration::getBank(const int bankId) const {
  if (!hasBank(bankId))
    throw std::runtime_error("BankCalibration: bank " + std::to_string(bankId) + " is not calibrated in '" + m_name + "'");
  return m_banks[bankId];
}

BankCalibration BankCalibration::load(const std::string& name) {
  if (name == nominalName)
    return BankCalibration();
  if (name.empty())
    throw std::runtime_error("BankCalibration: empty calibration name (use '" + nominalName + "' for the nominal geometry)");
  BankCalibration calibration;
  if (name.find('/') != std::string::npos) {
    calibration = loadFile(name, "");
  }
  else {
    const std::string dataFileName = "bank_calibration_" + name + ".txt";
    const std::string path = Core::findData("G4GeoLoki", dataFileName);
    if (path.empty())
      throw std::runtime_error("BankCalibration: unknown bank calibration '" + name + "' (no data file G4GeoLoki/" + dataFileName
                               + "; use '" + nominalName + "' or the name of a file in G4GeoLoki/data)");
    calibration = loadFile(path, name);
  }
  calibration.m_name = name;
  return calibration;
}

namespace {
  double dot(const std::array<double,3>& a, const std::array<double,3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }
  std::array<double,3> cross(const std::array<double,3>& a, const std::array<double,3>& b) {
    return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
  }
}

BankCalibration BankCalibration::loadFile(const std::string& fileName, const std::string& expectedName) {
  std::ifstream file(fileName);
  if (!file)
    throw std::runtime_error("BankCalibration: can not open the bank calibration file " + fileName);
  auto fail = [&fileName](int lineNumber, const std::string& message) {
    throw std::runtime_error("BankCalibration: " + fileName + ":" + std::to_string(lineNumber) + ": " + message);
  };

  BankCalibration calibration;
  calibration.m_fileName = fileName;
  bool formatSeen = false;
  bool nameSeen = false;
  std::string line;
  int lineNumber = 0;
  while (std::getline(file, line)) {
    lineNumber++;
    std::istringstream words(line);
    std::string key;
    if (!(words >> key) || key[0] == '#')
      continue; // empty or comment line
    if (!formatSeen) {
      int version = 0;
      if (key != "loki_bank_calibration" || !(words >> version))
        fail(lineNumber, "the first line must be 'loki_bank_calibration <version>'");
      if (version != 1)
        fail(lineNumber, "unsupported format version " + std::to_string(version) + " (supported: 1)");
      formatSeen = true;
    }
    else if (key == "name") {
      std::string fileCalibrationName;
      if (nameSeen || !(words >> fileCalibrationName))
        fail(lineNumber, "invalid or repeated 'name' line");
      if (!expectedName.empty() && fileCalibrationName != expectedName)
        fail(lineNumber, "the calibration name '" + fileCalibrationName + "' does not match the requested '" + expectedName + "'");
      nameSeen = true;
    }
    else if (key == "bank") {
      if (!nameSeen)
        fail(lineNumber, "the 'name' line must come before the bank lines");
      int bankId = -1;
      std::array<double,12> values;
      if (!(words >> bankId))
        fail(lineNumber, "invalid bank id");
      for (double& value : values)
        if (!(words >> value) || !std::isfinite(value))
          fail(lineNumber, "a bank line needs the bank id and 12 finite numbers (F, u, w, n)");
      std::string extra;
      if (words >> extra)
        fail(lineNumber, "unexpected value '" + extra + "' after the 12 numbers of the bank line");
      if (bankId < 0 || bankId > 8)
        fail(lineNumber, "invalid bank id " + std::to_string(bankId) + " (0-8)");
      if (calibration.m_hasBank[bankId])
        fail(lineNumber, "bank " + std::to_string(bankId) + " appears twice");

      Bank bank;
      for (int i = 0; i < 3; i++) {
        bank.frontFaceCentre[i] = values[i] * Units::mm;
        bank.alongTubes[i] = values[3 + i];
        bank.acrossTubes[i] = values[6 + i];
        bank.layerNormal[i] = values[9 + i];
      }
      // (u, w, n) must be a right-handed orthonormal basis, n pointing away from the sample
      const double tolerance = 1e-9;
      const auto& u = bank.alongTubes;
      const auto& w = bank.acrossTubes;
      const auto& n = bank.layerNormal;
      if (std::abs(dot(u, u) - 1) > tolerance || std::abs(dot(w, w) - 1) > tolerance || std::abs(dot(n, n) - 1) > tolerance)
        fail(lineNumber, "u, w and n must be unit vectors");
      if (std::abs(dot(u, w)) > tolerance || std::abs(dot(u, n)) > tolerance || std::abs(dot(w, n)) > tolerance)
        fail(lineNumber, "u, w and n must be perpendicular to each other");
      if (dot(cross(u, w), n) < 0)
        fail(lineNumber, "n must be u x w (right-handed basis), got -(u x w)");
      if (dot(n, bank.frontFaceCentre) <= 0)
        fail(lineNumber, "the layer normal n must point away from the sample");

      calibration.m_banks[bankId] = bank;
      calibration.m_hasBank[bankId] = true;
    }
    else {
      fail(lineNumber, "unknown line '" + key + "'");
    }
  }
  if (!formatSeen || !nameSeen)
    fail(lineNumber, "incomplete file (the 'loki_bank_calibration' and 'name' lines are required)");
  return calibration;
}
