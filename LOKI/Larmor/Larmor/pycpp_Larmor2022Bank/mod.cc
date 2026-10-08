#include "Core/Python.hh"
#include "Larmor/Larmor2022Bank.hh"

namespace {
  py::tuple toTuple( const std::array<double,3>& v )
  {
    return py::make_tuple( v[0], v[1], v[2] );
  }
  py::tuple pyLarmor2022Bank_getPixelCentre( const Larmor2022Bank& this_, const int pixelId, const bool oldTubeNumbering )
  {
    return toTuple( this_.getPixelCentre( pixelId, oldTubeNumbering ) );
  }
  // ((rotation rows), (translation)): p_world = rotation * p_bank + translation [mm]
  py::tuple pyLarmor2022Bank_getBankTransform( const Larmor2022Bank& this_ )
  {
    const auto t = this_.getBankTransform();
    const auto& R = t.rotation;
    return py::make_tuple( py::make_tuple( toTuple( R[0] ), toTuple( R[1] ), toTuple( R[2] ) ), toTuple( t.translation ) );
  }
  py::tuple pyLarmor2022Bank_getCalibMaskPosition( const Larmor2022Bank& this_ )
  {
    return toTuple( this_.getCalibMaskPosition() );
  }
}

PYTHON_MODULE( mod )
{
  mod.attr("DEFAULT_NUMBER_OF_PIXELS_IN_STRAW") = Larmor2022Bank::defaultNumberOfPixelsInStraw;

  py::class_<Larmor2022Bank>(mod, "Larmor2022Bank")
    .def(py::init<double>(), py::arg("rearBankDistance"))
    .def(py::init<double, int>(), py::arg("rearBankDistance"), py::arg("numberOfPixelsInStraw"))
    .def("getRearBankDistance", &Larmor2022Bank::getRearBankDistance)
    .def("getNumberOfPixelsInStraw", &Larmor2022Bank::getNumberOfPixelsInStraw)
    .def("getTotalNumberOfPixels", &Larmor2022Bank::getTotalNumberOfPixels)
    .def_static("getNumberOfPacks", &Larmor2022Bank::getNumberOfPacks)
    .def_static("getNumberOfTubes", &Larmor2022Bank::getNumberOfTubes)
    .def("getBankTransform", &pyLarmor2022Bank_getBankTransform)
    .def_static("getBankPositionY", &Larmor2022Bank::getBankPositionY)
    .def_static("getTubeVolumeNumber", &Larmor2022Bank::getTubeVolumeNumber)
    .def_static("getPackId", &Larmor2022Bank::getPackId)
    .def_static("getInPackTubeId", &Larmor2022Bank::getInPackTubeId)
    .def_static("getTubeLayerId", &Larmor2022Bank::getTubeLayerId)
    .def("getTubeId", &Larmor2022Bank::getTubeId)
    .def("getStrawId", &Larmor2022Bank::getStrawId)
    .def("getPixelCentre", &pyLarmor2022Bank_getPixelCentre, py::arg("pixelId"), py::arg("oldTubeNumbering") = false)
    .def("getPixelId", &Larmor2022Bank::getPixelId)
    .def("getCalibMaskPosition", &pyLarmor2022Bank_getCalibMaskPosition)
    ;
}
