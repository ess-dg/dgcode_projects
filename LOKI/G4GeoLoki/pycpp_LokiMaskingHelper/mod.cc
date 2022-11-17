#include "Core/Python.hh"
#include "G4GeoLoki/MaskingHelper.hh"

PYTHON_MODULE
{
  py::class_<MaskingHelper>("MaskingHelper", py::init<double>())
    .def(py::init<double, int>())
    .def(py::init<double, int, int>())
    .def("getPixelCentrePositionsForMasking",&MaskingHelper::getPixelCentrePositionsForMasking)
    .def("getTotalNumberOfPixels",&MaskingHelper::getTotalNumberOfPixels)
    // .def("getTotalNumberOfPixels",&MaskingHelper::getTotalNumberOfPixels).staticmethod("getTotalNumberOfPixels")
    ;
}
