"""Pixel positions and ids of the LoKI detector banks (see G4GeoLoki/AimHelper.hh)"""

from . import _LokiAimHelper as _impl

# the bank calibration names (see G4GeoLoki/BankCalibration.hh)
NOMINAL_BANK_CALIBRATION = _impl.nominal_bank_calibration()
DEFAULT_BANK_CALIBRATION = _impl.default_bank_calibration()
DEFAULT_NUMBER_OF_PIXELS_IN_STRAW = _impl.default_number_of_pixels_in_straw()

# the text of a bank calibration file ("" for the nominal geometry); throws for an unknown name
bankCalibrationText = _impl.bankCalibrationText

def AimHelper( rearBankDistance, strawPixelNumber = DEFAULT_NUMBER_OF_PIXELS_IN_STRAW,
               numberOfBanks = 9, bankCalibration = NOMINAL_BANK_CALIBRATION ):
    """Create an AimHelper object (the defaults are those of the C++ class)"""
    return _impl.AimHelper( float( rearBankDistance ), int( strawPixelNumber ),
                            int( numberOfBanks ), bankCalibration )
