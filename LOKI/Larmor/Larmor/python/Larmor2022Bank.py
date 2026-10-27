"""The Larmor 2022 detector bank (see Larmor/Larmor2022Bank.hh)"""

from . import _Larmor2022Bank as _impl

DEFAULT_NUMBER_OF_PIXELS_IN_STRAW = _impl.default_number_of_pixels_in_straw()

class Larmor2022Bank:
    """Larmor 2022 bank (wrapping the C++ class, to provide the default values of
    arguments)"""

    def __init__( self, rearBankDistance, numberOfPixelsInStraw = DEFAULT_NUMBER_OF_PIXELS_IN_STRAW ):
        b = _impl.Larmor2022Bank( float( rearBankDistance ), int( numberOfPixelsInStraw ) )
        self._b = b
        #Methods of the C++ object:
        for n in ( 'getRearBankDistance', 'getNumberOfPixelsInStraw', 'getTotalNumberOfPixels',
                   'getBankTransform', 'getTubeId', 'getStrawId', 'getPixelId', 'getCalibMaskPosition' ):
            setattr( self, n, getattr( b, n ) )

    def getPixelCentre( self, pixelId, oldTubeNumbering = False ):
        return self._b.getPixelCentre( pixelId, bool( oldTubeNumbering ) )

    getNumberOfPacks = staticmethod( _impl.Larmor2022Bank.getNumberOfPacks )
    getNumberOfTubes = staticmethod( _impl.Larmor2022Bank.getNumberOfTubes )
    getBankPositionY = staticmethod( _impl.Larmor2022Bank.getBankPositionY )
    getTubeVolumeNumber = staticmethod( _impl.Larmor2022Bank.getTubeVolumeNumber )
    getPackId = staticmethod( _impl.Larmor2022Bank.getPackId )
    getInPackTubeId = staticmethod( _impl.Larmor2022Bank.getInPackTubeId )
    getTubeLayerId = staticmethod( _impl.Larmor2022Bank.getTubeLayerId )
