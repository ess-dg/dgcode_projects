"""Test-only helpers exposing internals of the G4GeoLoki library and a minimal
Geant4 navigation utility (see the comments in pymod__TestUtils/mod.cc)"""

import numpy as _np
from . import _TestUtils as _impl
from G4GeoLoki.LokiAimHelper import NOMINAL_BANK_CALIBRATION as _NOMINAL

def _collect( fct, *args ):
    #Calls fct(*args,add) and returns the list of tuples passed to add:
    out = []
    fct( *args, lambda *a : out.append( a ) )
    return out

def construct_world( geo ):
    """Construct the world of a geometry module instance (e.g. created with
    G4GeoLoki.GeoBCSBanks.create()), and returns its number of daughters"""
    return _impl.construct_world( geo._g4handle )

def locate( x, y, z ):
    hist = []
    loc, ax, org = _impl.locate( x, y, z, lambda *a : hist.append( a ) )
    return hist, loc, ax, org

def world_daughters():
    return _collect( _impl.world_daughters )

def world_tree():
    return _collect( _impl.world_tree )

def dense_world_overlaps( npoints, tolerance ):
    return _collect( _impl.dense_world_overlaps, npoints, tolerance )

def dense_bank_containment( npoints, tolerance ):
    return _collect( _impl.dense_bank_containment, npoints, tolerance )

bank_daughters_distance_to_box = _impl.bank_daughters_distance_to_box

def MaskFileCreator( fileName, indexOffset, bankPixelLimits, aimingBankId = 0 ):
    """The mask file writer of the masking analysis (LokiMasking)"""
    return _impl.MaskFileCreator( fileName, indexOffset, _np.asarray( bankPixelLimits, dtype = _np.int64 ),
                                  aimingBankId )

class PixelCalc:
    """Thin wrapper around PixelatedBanks (the class used by the analysis programs
    to convert hits to pixel ids)"""

    def __init__( self, rear, n, nbanks, bankCalibration = _NOMINAL ):
        o = _impl.PixelCalc( float( rear ), int( n ), int( nbanks ), bankCalibration )
        self._o = o
        for name in ( 'getPixelId', 'getBankTransform', 'getNominalBankTransform', 'getBankNotch',
                      'isBankCalibrated', 'getBankCalibrationName', 'getBankCalibrationFile',
                      'getTotalNumberOfPixels', 'getNumberOfBanks' ):
            setattr( self, name, getattr( o, name ) )

    def getBankBoxOverlapPoints( self, bank, otherBank ):
        """The points spanning the overlap region of two bank boxes, in the frame of the first bank"""
        return [ tuple( p ) for p in self._o.getBankBoxOverlapPoints( bank, otherBank ).reshape( -1, 3 ).tolist() ]

#static getters (PixelatedBanks / BcsBanks / BcsPack / BcsTube):
getTubeLayerId = _impl.getTubeLayerId
getTubeIdInBank = _impl.getTubeIdInBank
getPackId = _impl.getPackId
getInPackTubeId = _impl.getInPackTubeId
getBankSize = _impl.getBankSize
getStrawLengthByBankId = _impl.getStrawLengthByBankId
getNumberOfPacksByBankId = _impl.getNumberOfPacksByBankId
getNumberOfTubes = _impl.getNumberOfTubes
getPackPositionInBank = _impl.getPackPositionInBank
getPackRotation = _impl.getPackRotation
getPackPackDistance = _impl.getPackPackDistance
detectorSystemFrontDistanceFromBankFront = _impl.detectorSystemFrontDistanceFromBankFront
isVertical = _impl.isVertical
areTubesInverselyNumbered = _impl.areTubesInverselyNumbered
getBeamstopSize = _impl.getBeamstopSize
getFrontFaceCentreInBank = _impl.getFrontFaceCentreInBank
getTubeRotation = _impl.getTubeRotation
getHorizontalTubeOffset = _impl.getHorizontalTubeOffset
getVerticalTubeOffset = _impl.getVerticalTubeOffset
getStrawPositionX = _impl.getStrawPositionX
getStrawPositionY = _impl.getStrawPositionY
getTubeOuterRadius = _impl.getTubeOuterRadius
getStrawInnerRadius = _impl.getStrawInnerRadius
