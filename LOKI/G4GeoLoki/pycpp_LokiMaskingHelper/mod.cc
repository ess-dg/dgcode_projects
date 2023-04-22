#include "Core/Python.hh"
#include "G4GeoLoki/MaskingHelper.hh"

PYTHON_MODULE
{
  py::class_<MaskingHelper>("MaskingHelper", py::init<double>())
    .def(py::init<double, int>())
    .def(py::init<double, int, int>())
    .def("getPixelCentrePositionsForMasking",&MaskingHelper::getPixelCentrePositionsForMasking)
    .def("getTotalNumberOfPixels",&MaskingHelper::getTotalNumberOfPixels)
    .def("getNumberOfPixels",&MaskingHelper::getNumberOfPixels).staticmethod("getNumberOfPixels")
    .def("getBankPixelOffset",&MaskingHelper::getBankPixelOffset).staticmethod("getBankPixelOffset")
    .def("getBankPosition",&MaskingHelper::getBankPosition)
    .def("dumpInfo",&MaskingHelper::dumpInfo).staticmethod("dumpInfo")
    .def("getBankId",&MaskingHelper::getBankId)
    .def("getPackId",&MaskingHelper::getPackId).staticmethod("getPackId")
    .def("getTubeId",&MaskingHelper::getTubeId).staticmethod("getTubeId")
    .def("getStrawId",&MaskingHelper::getStrawId).staticmethod("getStrawId")
    ;
}
