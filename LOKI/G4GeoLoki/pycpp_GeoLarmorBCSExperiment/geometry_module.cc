/////////////////////////////////////////
// Declaration of our geometry module: //
/////////////////////////////////////////

#include "G4Interfaces/GeoConstructPyExport.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Transform3D.hh"
#include "G4Vector3D.hh"
#include "G4SubtractionSolid.hh"
#include <cmath>
#include <string>

#include "G4GeoLoki/BcsBanks.hh"

class GeoBCS : public G4Interfaces::GeoConstructBase
{
public:
  GeoBCS();
  virtual ~GeoBCS(){}
  virtual G4VPhysicalVolume* Construct();

  BcsBanks * banks;
protected:
  virtual bool validateParameters();
private:
  //Functions
  G4LogicalVolume * createTubeLV(double converter_thickness, double straw_length);
  G4LogicalVolume * createPackBoxLV(double strawLength, int packNumber, int numberOfPacksForInvertedNumbering);
  G4LogicalVolume * createBankLV(int bankId);

  int getTubeVolumeNumber(int packNumber, int localIndex, int numberOfPacksForInvertedNumbering);
};

// this line is necessary to be able to declare the geometry in the python simulation script
PYTHON_MODULE( mod ) { GeoConstructPyExport::exportGeo<GeoBCS>(mod, "GeoLarmorBCSExperiment"); }

////////////////////////////////////////////
// Implementation of our geometry module: //
////////////////////////////////////////////

GeoBCS::GeoBCS()
  : GeoConstructBase("G4GeoLoki/GeoLarmorBCSExperiment"){
  // declare all parameters that can be used from the command line,
  // define their type, pick a self-explanatory name and _unit
  // give the default value, the last 2 ones are optional (min, max)

  addParameterDouble("rear_detector_distance_m", 4.420, 3.0, 10.0);
  addParameterBoolean("withBeamstop", false);
  addParameterDouble("generator_detector_distance_cm", 400, 0, 1000); //

  addParameterString("world_material","G4_Vacuum");
  addParameterString("B4C_panel_material","MAT_B4C:b10_enrichment=0.95"); //TODO  B4C S-DOUGH IS AN 18% EPOXY RESIN, 82% BORON CARBIDE MIX.IT IS MANUFACTURED BY STFC's ADVANCED MATERIALS GROUP.
}

G4LogicalVolume * GeoBCS::createTubeLV(double converterThickness, double strawLength){
  const double effectiveStrawLength = strawLength - BcsTube::getStrawWallThickness(); //This is only epsilon difference...

  auto lv_tube = new G4LogicalVolume(new G4Tubs("TubeWall",0,  BcsTube::getTubeOuterRadius(), 0.5*strawLength, 0., 2*M_PI),
                                     BcsTube::tubeWallMaterial, "TubeWall");

  auto lv_empty_tube = place(new G4Tubs("EmptyTube", 0., BcsTube::getTubeInnerRadius(), 0.5*strawLength, 0., 2*M_PI),
                             BcsTube::tubeInnerGas, 0,0,0, lv_tube, G4Colour(0,1,1),-2,0,0).logvol;

  for (int cpNo = 0; cpNo <= 6; cpNo++){
    auto lv_straw_wall = place(new G4Tubs("StrawWall", 0, BcsTube::getStrawOuterRadius(), 0.5*strawLength, 0., 2 * M_PI),
                               BcsTube::strawWallMaterial, BcsTube::getStrawPositionX(cpNo), BcsTube::getStrawPositionY(cpNo), 0, lv_empty_tube, ORANGE, cpNo, 0, 0).logvol;

    auto lv_converter = place(new G4Tubs("Converter", 0., BcsTube::getStrawInnerRadius(), 0.5*effectiveStrawLength, 0., 2 * M_PI),
                              BcsTube::converterMaterial, 0, 0, 0, lv_straw_wall, G4Colour(0, 1, 1), cpNo + 100, 0, 0).logvol;

    place(new G4Tubs("CountingGas", 0., BcsTube::getStrawInnerRadius() - converterThickness, 0.5*effectiveStrawLength, 0., 2 * M_PI),
          BcsTube::countingGas, 0, 0, 0, lv_converter, G4Colour(0, 0, 1), 0, 0, 0);
  }
  return lv_tube;
}

int GeoBCS::getTubeVolumeNumber(int packNumber, int inPackTubeId, int numberOfPacksForInvertedNumbering){
  assert(0 <= inPackTubeId && inPackTubeId <= 7);
  const double rowNumber = packNumber + 0.5* ((int)inPackTubeId/4); // +0.5 for second row (inPackTubeId > 3)
  inPackTubeId %= 4;
  if(numberOfPacksForInvertedNumbering == 0){
    return rowNumber * 8 + inPackTubeId;
  }
  else{
    return (numberOfPacksForInvertedNumbering*8 - 4) - rowNumber*8 + inPackTubeId;
  }
}

///////////  CREATE PACK BOX LOGICAL VOLUME  //////////////////////////
G4LogicalVolume *GeoBCS::createPackBoxLV(double strawLength, int packNumber, int numberOfPacksForInvertedNumbering){
  auto lv_pack_box = new G4LogicalVolume(new G4Box("EmptyPackBox", 0.5*BcsPack::getPackBoxWidth(), 0.5*BcsPack::getPackBoxHeight(), 0.5 * strawLength + BcsPack::getPackBoxIdleLengthOnOneEnd()), BcsPack::packBoxFillMaterial, "EmptyPackBox");
  /// Add 8 bcs detector tubes ///
  double frontB4CThickness = 0.65 *Units::um;
  double backB4CThickness = packNumber != 3 ? 1.00 *Units::um : 0.65 *Units::um;
  auto lv_front_tube = createTubeLV(frontB4CThickness, strawLength);
  auto lv_back_tube = createTubeLV(backB4CThickness, strawLength);
  G4RotationMatrix* tubeRotationMatrix = new G4RotationMatrix(0, 0, BcsPack::getTubeRotation());

  for (int inPackTubeId = 0; inPackTubeId < 8; inPackTubeId++) {
    place( (inPackTubeId % 4 < 2) ? lv_front_tube : lv_back_tube,
          BcsPack::getHorizontalTubeOffset(inPackTubeId), BcsPack::getVerticalTubeOffset(inPackTubeId), 0,
          lv_pack_box, SILVER, getTubeVolumeNumber(packNumber, inPackTubeId, numberOfPacksForInvertedNumbering), 0, tubeRotationMatrix);
  }
 /// Add B4C panel in 3 parts ///
  const double B4CLengthHalf = 0.5*strawLength + BcsPack::getB4CLengthOverStrawOnOneEnd();

  for (int i = 0; i < 3; i++){
    place(new G4Box("B4CPanel", 0.5*BcsPack::getB4CPartThickness(i), 0.5*BcsPack::getB4CPartHeight(i), B4CLengthHalf),
          BcsPack::B4CPanelMaterial,
          BcsPack::getB4CPartHorizontalOffset(i), BcsPack::getB4CPartVerticalOffset(i), 0,
          lv_pack_box, G4Colour(0, 1, 0), -2, 0, new G4RotationMatrix());
  }
  /// Add Al behing the B4C panel in 2 parts ///
  for (int i = 0; i < 2; i++){
    place(new G4Box("AlPanel", 0.5*BcsPack::getAlPartThickness(i), 0.5*BcsPack::getAlPartHeight(i), B4CLengthHalf),
          BcsPack::AlPanelMaterial,
          BcsPack::getAlPartHorizontalOffset(i),  BcsPack::getAlPartVerticalOffset(i), 0,
          lv_pack_box, SILVER, -2, 0, new G4RotationMatrix());
  }
  return lv_pack_box;
}

///////////  CREATE DETECTOR BANK LOGICAL VOLUME  //////////////////////////
G4LogicalVolume *GeoBCS::createBankLV(int bankId){
  const double pack_pack_distance = banks->getPackPackDistance();
  const int numberOfPacks = 4;

  const double packRotation = banks->getPackRotation();

  const double bankSizeXHalf = 0.5* banks->getBankSize(bankId, 0) * 1.5;
  const double bankSizeYHalf = 0.5* banks->getBankSize(bankId, 1) * 0.3;
  const double bankSizeZHalf = 0.5* banks->getBankSize(bankId, 2);

  auto lv_bank = new G4LogicalVolume(new G4Box("EmptyPanelBox", bankSizeZHalf, bankSizeYHalf, bankSizeXHalf),
                                     BcsPack::packBoxFillMaterial, "Bank");

  const int numberOfPacksForInvertedNumbering = 0;

  const double topmostBankPositionY = 2 * pack_pack_distance + 45.00 *Units::mm;
  for (int packNumber = 0; packNumber < numberOfPacks; ++packNumber){
    double strawLength = packNumber < 2 ? 1500.0 * Units::mm : 1200.0 * Units::mm;
    auto lv_pack_box = createPackBoxLV(strawLength, packNumber, numberOfPacksForInvertedNumbering);
    place(lv_pack_box,
    banks->getPackPositionInBank(bankId, packNumber, 2), topmostBankPositionY - packNumber * pack_pack_distance, 0,
    lv_bank, G4Colour(0, 1, 1), -2, 0, new G4RotationMatrix(0, 0, packRotation));
  }

  // Add BeamStop to the Rear Bank volume
  const bool withBeamstop = getParameterBoolean("withBeamstop");
  if ( withBeamstop) {
    const std::string maskName = "BoronMask-Beamstop";
    const double detBankFrontDistance = banks->detectorSystemFrontDistanceFromBankFront(bankId);
    auto emptyRotation = new G4RotationMatrix();

    place(new G4Box(maskName, 0.5 * 1*Units::mm, 0.5 * 5*Units::cm, 0.5 * 5*Units::cm),
          BoronMasks::maskMaterial,
          -bankSizeZHalf + detBankFrontDistance - 5*Units::cm, -banks->getBankPosition(bankId, 1), 40*Units::mm,
          lv_bank, BLACK, -5, 0, emptyRotation);
  }

  return lv_bank;
 }




G4VPhysicalVolume* GeoBCS::Construct(){
  // this is where we put the entire geometry together, the private functions creating the logical volumes are meant to facilitate the code below

  auto world_material = getParameterMaterial("world_material");
  const double generator_detector_distance = getParameterDouble("generator_detector_distance_cm")*Units::cm;
  const double sdd = getParameterDouble("rear_detector_distance_m")*Units::m;

  banks = new BcsBanks(sdd);

  // calculate a value that is big enough to fit your world volume, the "super mother"
  double big_dimension = 1.1*( 1 *Units::m + generator_detector_distance + sdd);

  //World volume:
  auto worldvols = place(new G4Box("World", big_dimension, big_dimension, big_dimension), world_material, 0, 0, 0, 0, INVISIBLE);
  auto lvWorld = worldvols.logvol;
  auto pvWorld = worldvols.physvol;

  //for (int bankId = 0; bankId < 9; bankId++){
  const int bankId = 0;
  auto lv_bank = createBankLV(bankId);

  auto rotation = new G4RotationMatrix();
  rotation->rotateY(banks->getBankRotation(bankId, 1));
  rotation->rotateX(banks->getBankRotation(bankId, 0));
  rotation->rotateZ(banks->getBankRotation(bankId, 2));

  place(lv_bank, 0.0 /*banks->getBankPosition(bankId, 0)*/, banks->getBankPosition(bankId, 1), banks->getBankPosition(bankId, 2), lvWorld, ORANGE, bankId, 0, rotation);
  //}

  delete banks;
  return pvWorld;
}


bool GeoBCS::validateParameters() {
// you can apply conditions to control the sanity of the geometry parameters and warn the user of possible mistakes
  // a nice example: Projects/SingleCell/G4GeoSingleCell/libsrc/GeoB10SingleCell.cc
    return true;
}
