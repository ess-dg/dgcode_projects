#include "sbldpy/sbldpy.hh"
#include "G4GeoLoki/AimHelper.hh"

//Private module used by the LokiAimHelper Python module (which provides the
//default values of the constructor arguments).

namespace {
  using T3 = std::array<double,3>;
}

SBLDPY_MODULE( _LokiAimHelper )
{
  // the bank calibration names (see G4GeoLoki/BankCalibration.hh)
  m.def( "nominal_bank_calibration", []() { return BankCalibration::nominalName; } );
  m.def( "default_bank_calibration", []() { return BankCalibration::defaultName; } );
  m.def( "default_number_of_pixels_in_straw", []() { return std::int64_t( PixelatedBanks::defaultNumberOfPixelsInStraw ); } );
  // the text of a bank calibration file ("" for the nominal geometry); throws for an unknown name
  m.def( "bankCalibrationText", []( const std::string& name ) { return BankCalibration::load( name ).text(); },
         { "name" } );

  m.add_class<AimHelper>( "AimHelper" )
    .def_init<double, std::int64_t, std::int64_t, std::string>( { "rearBankDistance", "strawPixelNumber",
                                                                  "numberOfBanks", "bankCalibration" } )
    .def( "getPixelCentreCoordinates", []( const AimHelper& a, std::int64_t pixelId )
    {
      auto xyz = a.getPixelCentreCoordinates( pixelId );
      return T3{ std::get<0>( xyz ), std::get<1>( xyz ), std::get<2>( xyz ) };
    }, { "pixelId" } )
    .def( "getTotalNumberOfPixels", []( const AimHelper& a ) { return a.getTotalNumberOfPixels(); } )
    .def( "getNumberOfPixels", []( const AimHelper& a, std::int64_t bankId ) { return a.getNumberOfPixels( bankId ); },
          { "bankId" } )
    .def( "getBankPixelOffset", []( const AimHelper& a, std::int64_t bankId ) { return a.getBankPixelOffset( bankId ); },
          { "bankId" } )
    .def( "getNumberOfPixelsInStraw", []( const AimHelper& a, std::int64_t bankId )
          { return a.getNumberOfPixelsInStraw( bankId ); }, { "bankId" } )
    // ((rotation rows), (translation)): p_world = rotation * p_bank + translation [mm]
    .def( "getBankTransform", []( const AimHelper& a, std::int64_t bankId )
    {
      const auto t = a.getBankTransform( bankId );
      return std::make_pair( t.rotation, t.translation );
    }, { "bankId" } )
    .def( "getBankCalibrationName", []( const AimHelper& a ) { return a.getBankCalibration().name(); } )
    .def( "isBankCalibrated", []( const AimHelper& a, std::int64_t bankId ) { return a.isBankCalibrated( bankId ); },
          { "bankId" } )
    .def( "dumpInfo", []( const AimHelper& a ) { a.dumpInfo(); } )
    .def( "getBankId", []( const AimHelper& a, std::int64_t pixelId ) { return a.getBankId( pixelId ); }, { "pixelId" } )
    .def_static( "getPackId", []( std::int64_t bankId, std::int64_t tubeId ) { return AimHelper::getPackId( bankId, tubeId ); },
                 { "bankId", "tubeId" } )
    .def( "getTubeId", []( const AimHelper& a, std::int64_t pixelId, std::int64_t bankId )
          { return a.getTubeId( pixelId, bankId ); }, { "pixelId", "bankId" } )
    .def( "getStrawId", []( const AimHelper& a, std::int64_t pixelId, std::int64_t bankId, std::int64_t tubeId )
          { return a.getStrawId( pixelId, bankId, tubeId ); }, { "pixelId", "bankId", "tubeId" } );
}
