#include "Core/Python.hh"
#include "G4GeoLoki/AimHelper.hh"

PYTHON_MODULE
{
  py::class_<AimHelper>("AimHelper", py::init<double>())
    .def(py::init<double, int>())
    .def(py::init<double, int, int>())
    .def("getPixelCentreCoordinates",&AimHelper::getPixelCentreCoordinates)
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
}
