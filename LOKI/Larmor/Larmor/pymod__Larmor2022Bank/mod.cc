#include "sbldpy/sbldpy.hh"
#include "Larmor/Larmor2022Bank.hh"

//Private module used by the Larmor2022Bank Python module (which provides the
//default value of the numberOfPixelsInStraw constructor argument).

SBLDPY_MODULE( _Larmor2022Bank )
{
  m.def( "default_number_of_pixels_in_straw", []() { return std::int64_t( Larmor2022Bank::defaultNumberOfPixelsInStraw ); } );

  using B = Larmor2022Bank;
  using I = std::int64_t;
  m.add_class<B>( "Larmor2022Bank" )
    .def_init<double, I>( { "rearBankDistance", "numberOfPixelsInStraw" } )
    .def( "getRearBankDistance", &B::getRearBankDistance )
    .def( "getNumberOfPixelsInStraw", &B::getNumberOfPixelsInStraw )
    .def( "getTotalNumberOfPixels", &B::getTotalNumberOfPixels )
    .def_static( "getNumberOfPacks", &B::getNumberOfPacks )
    .def_static( "getNumberOfTubes", &B::getNumberOfTubes )
    // ((rotation rows), (translation)): p_world = rotation * p_bank + translation [mm]
    .def( "getBankTransform", []( const B& b )
    {
      const auto t = b.getBankTransform();
      return std::make_pair( t.rotation, t.translation );
    } )
    .def_static( "getBankPositionY", &B::getBankPositionY )
    .def_static( "getTubeVolumeNumber", []( I packNumber, I inPackTubeId, bool oldTubeNumbering )
                 { return B::getTubeVolumeNumber( packNumber, inPackTubeId, oldTubeNumbering ); },
                 { "packNumber", "inPackTubeId", "oldTubeNumbering" } )
    .def_static( "getPackId", []( I tubeId, bool oldTubeNumbering ) { return B::getPackId( tubeId, oldTubeNumbering ); },
                 { "tubeId", "oldTubeNumbering" } )
    .def_static( "getInPackTubeId", []( I tubeId, bool oldTubeNumbering ) { return B::getInPackTubeId( tubeId, oldTubeNumbering ); },
                 { "tubeId", "oldTubeNumbering" } )
    .def_static( "getTubeLayerId", []( I tubeId, bool oldTubeNumbering ) { return B::getTubeLayerId( tubeId, oldTubeNumbering ); },
                 { "tubeId", "oldTubeNumbering" } )
    .def( "getTubeId", []( const B& b, I pixelId ) { return b.getTubeId( pixelId ); }, { "pixelId" } )
    .def( "getStrawId", []( const B& b, I pixelId ) { return b.getStrawId( pixelId ); }, { "pixelId" } )
    .def( "getPixelCentre", []( const B& b, I pixelId, bool oldTubeNumbering ) { return b.getPixelCentre( pixelId, oldTubeNumbering ); },
          { "pixelId", "oldTubeNumbering" } )
    .def( "getPixelId", []( const B& b, I tubeId, I strawId, double x, double y, double z )
          { return b.getPixelId( tubeId, strawId, x, y, z ); }, { "tubeId", "strawId", "x", "y", "z" } )
    .def( "getCalibMaskPosition", &B::getCalibMaskPosition );
}
