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



}

PYTHON_MODULE
{
#ifdef DGCODE_USEPYBIND11
  py::class_<AimHelper>(m, "AimHelper")
    .def(py::init<double>())
    .def(py::init<double, int>())
    .def(py::init<double, int, int>())
    .def("getPixelCentreCoordinates",&pyAimHelper_getPixelCentreCoordinates)
    .def("getTotalNumberOfPixels",&AimHelper::getTotalNumberOfPixels)
    .def_static("getNumberOfPixels",&AimHelper::getNumberOfPixels)
    .def_static("getBankPixelOffset",&AimHelper::getBankPixelOffset)
    .def("getBankPosition",&AimHelper::getBankPosition)
    .def_static("dumpInfo",&AimHelper::dumpInfo)
    .def("getBankId",&AimHelper::getBankId)
    .def_static("getPackId",&AimHelper::getPackId)
    .def_static("getTubeId",&AimHelper::getTubeId)
    .def_static("getStrawId",&AimHelper::getStrawId)
    ;
#else
  py::class_<AimHelper>("AimHelper", py::init<double>())
    .def(py::init<double, int>())
    .def(py::init<double, int, int>())
    .def("getPixelCentreCoordinates",&pyAimHelper_getPixelCentreCoordinates)
    .def("getTotalNumberOfPixels",&AimHelper::getTotalNumberOfPixels)
    .def("getNumberOfPixels",&AimHelper::getNumberOfPixels).staticmethod("getNumberOfPixels")
    .def("getBankPixelOffset",&AimHelper::getBankPixelOffset).staticmethod("getBankPixelOffset")
    .def("getBankPosition",&AimHelper::getBankPosition)
    .def("dumpInfo",&AimHelper::dumpInfo).staticmethod("dumpInfo")
    .def("getBankId",&AimHelper::getBankId)
    .def("getPackId",&AimHelper::getPackId).staticmethod("getPackId")
    .def("getTubeId",&AimHelper::getTubeId).staticmethod("getTubeId")
    .def("getStrawId",&AimHelper::getStrawId).staticmethod("getStrawId")
    ;
#endif
}
