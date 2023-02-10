#include "G4GeoLoki/BcsTube.hh"
#include "Core/Units.hh"
#include <cmath>
#include <cassert>

G4Material* BcsTube::tubeWallMaterial = NamedMaterialProvider::getMaterial("NCrystal:cfg=Al_sg225.ncmat");
G4Material* BcsTube::tubeInnerGas = NamedMaterialProvider::getMaterial("IdealGas:formula=0.8*Ar+0.2*CO2:pressure_atm=0.7");
G4Material* BcsTube::strawWallMaterial = NamedMaterialProvider::getMaterial("NCrystal:cfg=Cu_sg225.ncmat");
G4Material* BcsTube::converterMaterial = NamedMaterialProvider::getMaterial("ESS_B4C:b10_enrichment=0.95");
G4Material* BcsTube::countingGas = NamedMaterialProvider::getMaterial("IdealGas:formula=0.8*Ar+0.2*CO2:pressure_atm=0.7");

const double BcsTube::strawOuterRadius = 3.75 *Units::mm;
const double BcsTube::strawWallThickness = 0.0254 *Units::mm; //0.001 inch

double BcsTube::getStrawOuterRadius(){
  return strawOuterRadius;
}
double BcsTube::getStrawInnerRadius(){
  return strawOuterRadius - strawWallThickness;
}
double BcsTube::getStrawWallThickness(){
  return strawWallThickness;
}


const double BcsTube::frontTubeConverterThickness = 0.65 *Units::um;
const double BcsTube::backTubeConverterThickness = 1.00 *Units::um;

double BcsTube::getFrontTubeConverterThickness(){
  return frontTubeConverterThickness;
}
double BcsTube::getBackTubeConverterThickness(){
  return backTubeConverterThickness;
}


const double BcsTube::tubeOuterRadius = 12.7 *Units::mm; // 1 inch diameter
const double BcsTube::tubeWallThickness = 0.889 *Units::mm; //0.035 inch  //0.94

double BcsTube::getTubeOuterRadius() {
  return tubeOuterRadius;
}
double BcsTube::getTubeInnerRadius() {
  return tubeOuterRadius - tubeWallThickness;
}


const double BcsTube::strawStrawDistance = 7.75 *Units::mm;
const double BcsTube::strawPositionsInTube[7][2] = { //TODO reorder to change volume numbering
    { -0.5*strawStrawDistance * tan(M_PI/3.), -0.5*strawStrawDistance }, //20
    { 0, -strawStrawDistance }, //0
    { 0.5*strawStrawDistance * tan(M_PI/3.), -0.5*strawStrawDistance }, //10
    { 0, 0},  //30
    { -0.5*strawStrawDistance * tan(M_PI/3.), 0.5*strawStrawDistance }, //50
    { 0, strawStrawDistance}, //60
    { 0.5*strawStrawDistance * tan(M_PI/3.), 0.5*strawStrawDistance }, //40
  };

double BcsTube::getStrawPositionX(const int strawId) {
  assert(0 <= strawId && strawId <= 6);
  return strawPositionsInTube[strawId][0];
}
double BcsTube::getStrawPositionY(const int strawId) {
  assert(0 <= strawId && strawId <= 6);
  return strawPositionsInTube[strawId][1];
}
