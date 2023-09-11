#include "G4GeoLoki/BcsPack.hh"
#include "G4Units/Units.hh"
#include <cmath>
#include <cassert>

const double BcsPack::tubeGridParallelogramBase = 27.00 *Units::mm;
const double BcsPack::tubeGridParallelogramSide = 28.40 *Units::mm;
const double BcsPack::tubeGridParallelogramAngle = 13.45 *Units::degree; // 90-76.55=13.45

const double BcsPack::tubeRotationAngle = (8.45-30) *Units::degree; //8.45 degree should be the final roataion after the pack is rotated, pack rotation must be subtracted // -30 degree needed because of different default positioning;

const double BcsPack::packBoxWidth = 265.7 *Units::mm; //excluding the handle
const double BcsPack::packBoxHeight = 55.20 *Units::mm;
const double BcsPack::packBoxIdleLengthOnOneEnd = 111.50 *Units::mm;
const double BcsPack::tubeCentreDistanceFromPackFront = 36.50 *Units::mm; //36.50=23.8+12.7 (packFrontToTubeFront+tubeOuterRadius)


double BcsPack::getHorizontalTubeDistanceInPack() {
  return tubeGridParallelogramBase;
}

double BcsPack::getVerticalTubeDistanceInPack() {  //28.40*cos(13.45 degree)=27.62
  return tubeGridParallelogramSide * std::cos(tubeGridParallelogramAngle);
}

double BcsPack::getTubeGridParallelogramAngle() {
  return tubeGridParallelogramAngle;
}
double BcsPack::getTubeRotation() { //compensated for the rotation of the packs
  return tubeRotationAngle - tubeGridParallelogramAngle;
}


double BcsPack::getTopRowOffsetInPack() {
  return tubeGridParallelogramSide * std::sin(tubeGridParallelogramAngle);
}

double BcsPack::getHorizontalTubeCentreOffsetInPack() {
  return -0.5*packBoxWidth + tubeCentreDistanceFromPackFront;
}

double BcsPack::getTubeGridParallelogramSide() {
  return tubeGridParallelogramSide;
}
double BcsPack::getTubeCentreDistanceFromPackFront() {
  return tubeCentreDistanceFromPackFront;
}


G4Material* BcsPack::packBoxFillMaterial = NamedMaterialProvider::getMaterial("G4_Vacuum");
double BcsPack::getPackBoxWidth() {
  return packBoxWidth;
}
double BcsPack::getPackBoxHeight() {
  return packBoxHeight;
}
double BcsPack::getPackBoxIdleLengthOnOneEnd() {
  return packBoxIdleLengthOnOneEnd;
}

/// Tube positioning in pack ///

double BcsPack::getHorizontalTubeOffset(const int inPackTubeId) {
  assert(0 <= inPackTubeId && inPackTubeId <= 7);
  const int inRowTubeId = inPackTubeId %4;
  return getHorizontalTubeCentreOffsetInPack() +
         inRowTubeId * getHorizontalTubeDistanceInPack() +
         (inPackTubeId < 4 ? getTopRowOffsetInPack() : 0.0);
}
double BcsPack::getVerticalTubeOffset(const int inPackTubeId) {
  assert(0 <= inPackTubeId && inPackTubeId <= 7);
  return 0.5 * getVerticalTubeDistanceInPack() * (inPackTubeId < 4 ? 1 : -1);
}


/// B4C panel parts ///
G4Material* BcsPack::B4CPanelMaterial = NamedMaterialProvider::getMaterial("ESS_B4C:b10_enrichment=0.95"); //TODO B4C S-DOUGH IS AN 18% EPOXY RESIN, 82% BORON CARBIDE MIX.IT IS MANUFACTURED BY STFC's ADVANCED MATERIALS GROUP.
const double BcsPack::B4CLengthOverStrawOnOneEnd = 12.0 *Units::mm;
const double BcsPack::B4CDistanceFromLastTubeCentre = 23.0 *Units::mm;

const double BcsPack::B4CPanelPartThickness[3] = { // all in mm
    5.0, // 0 - main part
    10.80,  // 1 - middle part
    3.00, // 2 - bottom part
};
const double BcsPack::B4CPanelPartHeight[3] = { // all in mm
    51.20, // 0 - main part
    5.00,  // 1 - middle part
    3.20, // 2 - bottom part
};

double BcsPack::getB4CLengthOverStrawOnOneEnd() {
  return B4CLengthOverStrawOnOneEnd;
}

double BcsPack::getB4CPartThickness(const int partId){
  assert(0 <= partId && partId <= 2);
  return B4CPanelPartThickness[partId] *Units::mm;
}

double BcsPack::getB4CPartHeight(const int partId) {
  assert(0 <= partId && partId <= 2);
  return B4CPanelPartHeight[partId] *Units::mm;
}

double BcsPack::getB4CPartHorizontalOffset(const int partId) {
  assert(0 <= partId && partId <= 2);
  if(partId == 0) {
    return getHorizontalTubeCentreOffsetInPack() + 3.0* tubeGridParallelogramBase + B4CDistanceFromLastTubeCentre + 0.5* getB4CPartThickness(0);
  }
  else if (partId == 1) {
    return getB4CPartHorizontalOffset(0) - 0.5* getB4CPartThickness(0) - 0.5* getB4CPartThickness(1);
  }
  else {
    return getB4CPartHorizontalOffset(1) - 0.5* getB4CPartThickness(1) + 0.5* getB4CPartThickness(2);
  }
}

double BcsPack::getB4CPartVerticalOffset(const int partId) {
  assert(0 <= partId && partId <= 2);
  if(partId == 0) {
    return 0.5* getB4CPartHeight(2);
  }
  else if (partId == 1) {
    return getB4CPartVerticalOffset(0) - 0.5* getB4CPartHeight(0) + 0.5* getB4CPartHeight(1);
  }
  else {
    return getB4CPartVerticalOffset(1) - 0.5* getB4CPartHeight(1) - 0.5* getB4CPartHeight(2);
  }
}

/// Al panel parts ///
G4Material* BcsPack::AlPanelMaterial = NamedMaterialProvider::getMaterial("NCrystal:cfg=Al_sg225.ncmat");
double BcsPack::getAlPartThickness(const int partId){
  assert(0 <= partId && partId <= 1);
  if(partId==0) {
    return getB4CPartHeight(2);
  }
  else {
    return getB4CPartThickness(1) - getB4CPartThickness(2) + getB4CPartThickness(0);
  }
}

double BcsPack::getAlPartHeight(const int partId){
  assert(0 <= partId && partId <= 1);
  if(partId==0) {
    return getB4CPartHeight(0) + getB4CPartHeight(2);
  }
  else {
    return getB4CPartHeight(2);
  }
}

double BcsPack::getAlPartHorizontalOffset(const int partId){
  assert(0 <= partId && partId <= 1);
  if(partId==0) {
    return getB4CPartHorizontalOffset(0) + 0.5*getB4CPartThickness(0)+ 0.5*getAlPartThickness(0);
  }
  else {
    return getAlPartHorizontalOffset(0) - 0.5*getAlPartThickness(0) - 0.5*getAlPartThickness(1);
  }
}

double BcsPack::getAlPartVerticalOffset(const int partId){
  assert(0 <= partId && partId <= 1);
  if(partId==0) {
    return 0.0;
  }
  else {
    return getB4CPartVerticalOffset(2);
  }
}
