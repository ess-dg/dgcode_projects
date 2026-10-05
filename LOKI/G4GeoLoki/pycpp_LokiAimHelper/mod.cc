#include "Core/Python.hh"
#include "G4GeoLoki/AimHelper.hh"

namespace {
  py::tuple pyAimHelper_getPixelCentreCoordinates( const AimHelper& this_,
                                                    const int pixelId,
                                                    const bool isOldPixelNumbering,
                                                    const bool isLarmor2022Experiment)
  {
    auto xyz = this_.getPixelCentreCoordinates( pixelId,
                                                isOldPixelNumbering,
                                                isLarmor2022Experiment );
    return py::make_tuple( std::get<0>( xyz ),
                           std::get<1>( xyz ),
                           std::get<2>( xyz ) );
  }
  // ((rotation rows), (translation)): p_world = rotation * p_bank + translation [mm]
  py::tuple pyAimHelper_getBankTransform( const AimHelper& this_, const int bankId, const bool isLarmor2022Experiment )
  {
    const auto t = this_.getBankTransform( bankId, isLarmor2022Experiment );
    const auto& R = t.rotation;
    return py::make_tuple( py::make_tuple( py::make_tuple( R[0][0], R[0][1], R[0][2] ),
                                           py::make_tuple( R[1][0], R[1][1], R[1][2] ),
                                           py::make_tuple( R[2][0], R[2][1], R[2][2] ) ),
                           py::make_tuple( t.translation[0], t.translation[1], t.translation[2] ) );
  }
  std::string pyAimHelper_getBankCalibrationName( const AimHelper& this_ )
  {
    return this_.getBankCalibration().name();
  }
}

PYTHON_MODULE( mod )
{
  // the bank calibration names (see G4GeoLoki/BankCalibration.hh)
  mod.attr("NOMINAL_BANK_CALIBRATION") = BankCalibration::nominalName;
  mod.attr("DEFAULT_BANK_CALIBRATION") = BankCalibration::defaultName;

  py::class_<AimHelper>(mod, "AimHelper")
    .def(py::init<double>())
    .def(py::init<double, int>())
    .def(py::init<double, int, int>())
    .def(py::init<double, int, int, std::string>())
    .def("getPixelCentreCoordinates",&pyAimHelper_getPixelCentreCoordinates)
    .def("getTotalNumberOfPixels",&AimHelper::getTotalNumberOfPixels)
    .def_static("getNumberOfPixels",&AimHelper::getNumberOfPixels)
    .def_static("getBankPixelOffset",&AimHelper::getBankPixelOffset)
    .def_static("getNumberOfPixelsInStraw",&AimHelper::getNumberOfPixelsInStraw)
    .def("getBankPosition",&AimHelper::getBankPosition)
    .def("getBankTransform",&pyAimHelper_getBankTransform, py::arg("bankId"), py::arg("isLarmor2022Experiment") = false)
    .def("getBankCalibrationName",&pyAimHelper_getBankCalibrationName)
    .def("isBankCalibrated",&AimHelper::isBankCalibrated)
    .def_static("dumpInfo",&AimHelper::dumpInfo)
    .def("getBankId",&AimHelper::getBankId)
    .def_static("getPackId",&AimHelper::getPackId)
    .def_static("getTubeId",&AimHelper::getTubeId)
    .def_static("getStrawId",&AimHelper::getStrawId)
    ;
}
