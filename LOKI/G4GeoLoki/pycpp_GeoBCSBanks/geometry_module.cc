/////////////////////////////////////////
// Declaration of our geometry module: //
/////////////////////////////////////////

#include "G4Interfaces/GeoConstructPyExport.hh"
#include "G4Para.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Transform3D.hh"
#include "G4Vector3D.hh"
#include "G4SubtractionSolid.hh"
#include <cmath>
#include <string>
#include <cassert>

#include "G4GeoLoki/BcsBanks.hh"

class GeoBCS : public G4Interfaces::GeoConstructBase
{
public:
  GeoBCS();
  virtual ~GeoBCS(){}
  virtual G4VPhysicalVolume* Construct();

protected:
  virtual bool validateParameters();
private:
  BcsBanks* banks;
  //Functions
  G4LogicalVolume * createTubeLV(double converter_thickness, double straw_length);
  G4LogicalVolume * createPackBoxLV(double strawLength, int packNumber, int numberOfPacksForInvertedNumbering, int numberOfPacks);
  G4LogicalVolume * createBankLV(int bankId);
  G4LogicalVolume * createTriangularMaskLV(int maskId);
  G4LogicalVolume * createCalibrationMaskLV(CalibMasks::CalibMasksBase calibMask);

  int getTubeVolumeNumber(int packNumber, int inPackTubeId, int numberOfPacksForInvertedNumbering, int numberOfPacks);
};

// this line is necessary to be able to declare the geometry in the python simulation script
PYTHON_MODULE { GeoConstructPyExport::exportGeo<GeoBCS>("GeoBCSBanks"); }

////////////////////////////////////////////
// Implementation of our geometry module: //
////////////////////////////////////////////

GeoBCS::GeoBCS()
  : GeoConstructBase("G4GeoLoki/GeoBCSBanks"){
  // declare all parameters that can be used from the command line,
  addParameterDouble("rear_detector_distance_m", 5.0, 4.0, 10.0); // default, min, max
  addParameterBoolean("with_beamstop", false);
  addParameterBoolean("larmor_2022_experiment", false);
  addParameterBoolean("with_calibration_slits", false);

  addParameterBoolean("old_tube_numbering", false);

  addParameterString("world_material","G4_Vacuum");
  addParameterString("B4C_panel_material","ESS_B4C:b10_enrichment=0.95");
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

int GeoBCS::getTubeVolumeNumber(int packNumber, int inPackTubeId, int numberOfPacksForInvertedNumbering, int numberOfPacks){
  assert(0 <= inPackTubeId && inPackTubeId <= 7);
  const double rowNumber = packNumber + 0.5* ((int)inPackTubeId/4); // +0.5 for second row (inPackTubeId > 3)
  inPackTubeId %= 4;

  const bool oldTubeNumbering = getParameterBoolean("old_tube_numbering");
  if(oldTubeNumbering){
    if(numberOfPacksForInvertedNumbering == 0){
      return rowNumber * 8 + inPackTubeId;
    }
    else{
      return (numberOfPacksForInvertedNumbering*8 - 4) - rowNumber*8 + inPackTubeId;
    }
  }
  else{ // new numbering
    const int layerNr = inPackTubeId; //[0-3]
    if(numberOfPacksForInvertedNumbering == 0){
      return rowNumber * 2 + layerNr * numberOfPacks * 2;
    }
    else{
      return (numberOfPacks*2 - 1) - rowNumber*2 + layerNr * numberOfPacks * 2;
    }
  }
}

///////////  CREATE PACK BOX LOGICAL VOLUME  //////////////////////////
G4LogicalVolume *GeoBCS::createPackBoxLV(double strawLength, int packNumber, int numberOfPacksForInvertedNumbering, int numberOfPacks){
  const double packRotation = banks->getPackRotation();
  // Instead of a rectangular box, a detector pack is encapsulated in parallelepiped, to avoid collision of the corners with the calibraion slits after applying the pack rotation.
  // The PackBoxWidth corresponds to the size of the volume encapsulating the electronics on the sides as well, not just the detectors, but that would cause collision with the calibration slits, so a multiplication factor of 0.799 is applied, to get a volume just large enough to fit in the detectors in the front.
  auto lv_pack_box = new G4LogicalVolume(
    new G4Para("EmptyPackBox", 0.799*0.5*BcsPack::getPackBoxWidth(), 0.5*BcsPack::getPackBoxHeight(), 0.5 * strawLength + BcsPack::getPackBoxIdleLengthOnOneEnd(), packRotation, 0, 0),
    BcsPack::packBoxFillMaterial, "EmptyPackBox");

  /// Add 8 BCS detector tubes ///
  auto lv_front_tube = createTubeLV(BcsTube::getFrontTubeConverterThickness(), strawLength);
  auto lv_back_tube = createTubeLV(BcsTube::getBackTubeConverterThickness(), strawLength);
  G4RotationMatrix* tubeRotationMatrix = new G4RotationMatrix(0, 0, BcsPack::getTubeRotation());

  for (int inPackTubeId = 0; inPackTubeId < 8; inPackTubeId++) {
    place((inPackTubeId % 4 < 2) ? lv_front_tube : lv_back_tube,
          BcsPack::getHorizontalTubeOffset(inPackTubeId), BcsPack::getVerticalTubeOffset(inPackTubeId), 0,
          lv_pack_box, SILVER, getTubeVolumeNumber(packNumber, inPackTubeId, numberOfPacksForInvertedNumbering, numberOfPacks), 0, tubeRotationMatrix);
  }
  /// Add B4C panel behind detectors in 3 parts ///
  const double B4CLengthHalf = 0.5*strawLength + BcsPack::getB4CLengthOverStrawOnOneEnd();

  for (int partId = 0; partId < 3; partId++){
    place(new G4Box("B4CPanel", 0.5*BcsPack::getB4CPartThickness(partId), 0.5*BcsPack::getB4CPartHeight(partId), B4CLengthHalf),
          BcsPack::B4CPanelMaterial,
          BcsPack::getB4CPartHorizontalOffset(partId), BcsPack::getB4CPartVerticalOffset(partId), 0,
          lv_pack_box, G4Colour(0, 1, 0), -2, 0, new G4RotationMatrix());
    }
  /// Add Al behing the B4C panel in 2 parts ///
  for (int partId = 0; partId < 2; partId++){
    place(new G4Box("AlPanel", 0.5*BcsPack::getAlPartThickness(partId), 0.5*BcsPack::getAlPartHeight(partId), B4CLengthHalf),
          BcsPack::AlPanelMaterial,
          BcsPack::getAlPartHorizontalOffset(partId),  BcsPack::getAlPartVerticalOffset(partId), 0,
          lv_pack_box, SILVER, -2, 0, new G4RotationMatrix());
    }
  return lv_pack_box;
}

///////////  CREATE CALIBRATION SLIT LOGICAL VOLUME  //////////////////////////
G4LogicalVolume *GeoBCS::createCalibrationMaskLV(CalibMasks::CalibMasksBase calibMask){
  const std::string maskName = "BoronMask-"+calibMask.getName();
  const double maskThicknessHalf = 0.5*calibMask.getThickness();
  const double maskHeightHalf = 0.5*calibMask.getHeight();
  const double maskFullWidthHalf = 0.5*calibMask.getWidth();

  auto lv_calibrationMask = new G4LogicalVolume(new G4Box("EmptyCalibMaskBox", maskThicknessHalf, maskHeightHalf, maskFullWidthHalf), CalibMasks::maskBoxMaterial, "CalibMaskBox");

  double offset = 0.0;
  int i = 0;
  const auto pattern = calibMask.getPattern();
  for(auto part = pattern.begin(); part != pattern.end(); part++,i++ ) {
    const double partWidth = *part;
    if(i % 2 == 0) { //The pattern is: maskPart ,slit, maskPart, slit...
      place(new G4Box(maskName, maskThicknessHalf, maskHeightHalf, 0.5*partWidth),
          CalibMasks::maskMaterial,
          0., 0., -maskFullWidthHalf + offset + 0.5*partWidth,
          lv_calibrationMask, DARKPURPLE, -5, 0, new G4RotationMatrix());
    }
    offset += partWidth;
  }
  return lv_calibrationMask;
}


///////////  CREATE DETECTOR BANK LOGICAL VOLUME  //////////////////////////
G4LogicalVolume *GeoBCS::createBankLV(int bankId){
  const bool larmor2022experiment = getParameterBoolean("larmor_2022_experiment");
  const double strawLength = banks->getStrawLengthByBankId(bankId);

  // const double pack_pack_distance = banks->getPackPackDistance();
  const int numberOfPacks = banks->getNumberOfPacksByBankId(bankId);

  const double packRotation = banks->getPackRotation();

  const double bankSizeXHalf = 0.5* banks->getBankSize(bankId, 0);
  const double bankSizeYHalf = 0.5* banks->getBankSize(bankId, 1);
  const double bankSizeZHalf = 0.5* banks->getBankSize(bankId, 2);

  auto lv_bank = new G4LogicalVolume(new G4Box("EmptyPanelBox", bankSizeZHalf, bankSizeYHalf, bankSizeXHalf),
                                     BcsPack::packBoxFillMaterial, "Bank");

  // override lv_bank to subtract some empty part of front left and right banks where front top and bottom banks would overlap with them
  if (bankId == 8 || bankId == 6) {
    auto fullBankBox = new G4Box("EmptyPanelBox", bankSizeZHalf, bankSizeYHalf, bankSizeXHalf);
    auto bankBoxCut = new G4Box("EmptyPanelBox", 15.0, 70.0, 60.0);
    auto bankBox = new G4SubtractionSolid("EmptyPanelBox", fullBankBox, bankBoxCut, 0, G4ThreeVector(-145.0, -42.0, -425.0));

    lv_bank = new G4LogicalVolume(bankBox, BcsPack::packBoxFillMaterial, "Bank");
  }

  int numberOfPacksForInvertedNumbering = 0;
  if (banks->areTubesInverselyNumbered(bankId)){ //not very nice solution...
    numberOfPacksForInvertedNumbering = numberOfPacks;
  }

  for (int packNumber = 0; packNumber < numberOfPacks; ++packNumber){
    auto lv_pack_box = createPackBoxLV(strawLength, packNumber, numberOfPacksForInvertedNumbering, numberOfPacks);
    place(lv_pack_box,
          banks->getPackPositionInBank(bankId, packNumber, 2), banks->getPackPositionInBank(bankId, packNumber, 1), banks->getPackPositionInBank(bankId, packNumber, 0),
          lv_bank, G4Colour(0, 1, 1), -2, 0, new G4RotationMatrix(0, 0, packRotation));
  }

  const int numberOfBoronMasks = BoronMasks::getNumberOfBoronMasks(bankId);
  for (int maskId = 0; maskId < numberOfBoronMasks; ++maskId){
    const std::string maskName = "BoronMask-"+std::to_string(bankId)+"-"+std::to_string(maskId);
    place(new G4Box(maskName, 0.5*BoronMasks::getSize(bankId, maskId, 2), 0.5*BoronMasks::getSize(bankId, maskId, 1), 0.5*BoronMasks::getSize(bankId, maskId, 0)),
            BoronMasks::maskMaterial,
            banks->getBoronMaskPosition(bankId, maskId, 2), banks->getBoronMaskPosition(bankId, maskId, 1), banks->getBoronMaskPosition(bankId, maskId, 0),
            lv_bank, BLACK, -2, 0, new G4RotationMatrix(0, 0, BoronMasks::getRotation(bankId, maskId)));
  }

  // Add BeamStop to the Rear Bank volume
  const bool withBeamstop = getParameterBoolean("with_beamstop");
  if (bankId == 0 && withBeamstop) {
    const std::string maskName = "BoronMask-Beamstop";
    const double detBankFrontDistance = banks->detectorSystemFrontDistanceFromBankFront(bankId);
    const double verticalPosition = larmor2022experiment ? -banks->getLarmor2022ExperimentBankPositionY() : 0;

    place(new G4Box(maskName, 0.5 * 1*Units::mm, 0.5 * 5*Units::cm, 0.5 * 5*Units::cm),
          BoronMasks::maskMaterial,
          -bankSizeZHalf + detBankFrontDistance - 5*Units::cm, verticalPosition, 0,
          lv_bank, BLACK, -5, 0, new G4RotationMatrix());
  }

  const bool withCalibrationSlits = getParameterBoolean("with_calibration_slits");
  if (withCalibrationSlits && !larmor2022experiment) {
    std::string calibMaskName = "lokiStandard-"+std::to_string(bankId);
    const auto calibMask = CalibMasks::getCalibMask(calibMaskName);
    auto lv_calibrationMask = createCalibrationMaskLV(calibMask);

    place(lv_calibrationMask,
          banks->getCalibMaskPosition(calibMask, bankId, 2), banks->getCalibMaskPosition(calibMask, bankId, 1), banks->getCalibMaskPosition(calibMask, bankId, 0),
          lv_bank, PURPLE, -5, 0, new G4RotationMatrix());
  }

  return lv_bank;
 }

G4LogicalVolume *GeoBCS::createTriangularMaskLV(int maskId){
    //creating special triangular volume by subtraction of 2 boxes
    const double xHalf = BoronMasks::getHalfSizeOfTriangularMask(maskId, 0);
    const double yHalf = BoronMasks::getHalfSizeOfTriangularMask(maskId, 1);
    const double zHalf = BoronMasks::getHalfSizeOfTriangularMask(maskId, 2);
    const double xSideToCut = 2 * xHalf - BoronMasks::getCutPointOfTriangularMask(maskId, 0);
    const double ySideToCut = 2 * yHalf - BoronMasks::getCutPointOfTriangularMask(maskId, 1);
    const double xCutDir = BoronMasks::getCutDirOfTriangularMask(maskId, 0);
    const double yCutDir = BoronMasks::getCutDirOfTriangularMask(maskId, 1);

    auto triangularMaskBox = new G4Box("TriangularMaskBox", xHalf, yHalf, zHalf);
    const double xCutHalf = sqrt(pow(ySideToCut, 2) + pow(xSideToCut, 2)) / 2;
    const double yCutHalf = ySideToCut * xSideToCut / (2 * xCutHalf) / 2;
    const double alpha = atan(ySideToCut / xSideToCut); //rad

    auto triangularMaskCut = new G4Box("TriangularMaskSustractBox", xCutHalf, yCutHalf, zHalf * 1.1);
    G4RotationMatrix *cutRotationMatrix = new G4RotationMatrix(0, 0, -alpha * xCutDir * yCutDir);
    auto cutCentre = G4ThreeVector(xCutDir * (xHalf - xCutHalf * cos(alpha) + yCutHalf * sin(alpha)), yCutDir * (-(ySideToCut - yHalf) + xCutHalf * sin(alpha) + yCutHalf * cos(alpha)), 0);
    auto triangularMask = new G4SubtractionSolid("TriangularMask", triangularMaskBox, triangularMaskCut, cutRotationMatrix, cutCentre);

    const int bankId = BoronMasks::getBankIdOfTriangularMask(maskId);
    const std::string maskName = "BoronMask-triangular-"+std::to_string(bankId)+"-"+std::to_string(maskId);
    auto lv_triangularMask = new G4LogicalVolume(triangularMask, BoronMasks::maskMaterial, maskName);

    return lv_triangularMask;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* GeoBCS::Construct(){
  // this is where we put the entire geometry together, the private functions creating the logical volumes are meant to facilitate the code below
  const double rear_detector_distance = getParameterDouble("rear_detector_distance_m")*Units::m;
  const bool larmor2022experiment = getParameterBoolean("larmor_2022_experiment");
  const int numberOfBanks = larmor2022experiment ? 1 : 9;
  banks = new BcsBanks(rear_detector_distance, numberOfBanks);

  // calculate a value that is big enough to fit your world volume, the "super mother"
  double big_dimension = 1.1*( 1 *Units::m + rear_detector_distance);

  //World volume:
  auto world_material = getParameterMaterial("world_material");
  auto worldvols = place(new G4Box("World", big_dimension, big_dimension, big_dimension), world_material, 0, 0, 0, 0, INVISIBLE);
  auto lvWorld = worldvols.logvol;
  auto pvWorld = worldvols.physvol;

  // Create and place detector banks
  for (int bankId = 0; bankId < banks->getNumberOfBanks(); bankId++){
    auto lv_bank = createBankLV(bankId);

    auto rotation = new G4RotationMatrix();
    rotation->rotateY(banks->getBankRotation(bankId, 1));
    rotation->rotateX(banks->getBankRotation(bankId, 0));
    rotation->rotateZ(banks->getBankRotation(bankId, 2));

    const double verticalBankPosition = !larmor2022experiment ? banks->getBankPosition(bankId, 1) : banks->getLarmor2022ExperimentBankPositionY();

    place(lv_bank, banks->getBankPosition(bankId, 0), verticalBankPosition, banks->getBankPosition(bankId, 2), lvWorld, ORANGE, bankId, 0, rotation);
  }

  // Add 4 triangular boron masks (added to the World instead of the banks)
  if (!larmor2022experiment) {
    for (int maskId = 0; maskId <= 3; maskId++) {
      auto lv_triangularMask = createTriangularMaskLV(maskId);

      const int bankId = BoronMasks::getBankIdOfTriangularMask(maskId);       // 5 or 7
      const double rotateDir = BoronMasks::getCutDirOfTriangularMask(maskId, 1); // +- 1.0

      auto rotation = new G4RotationMatrix();
      rotation->rotateX(banks->getBankRotation(bankId, 2) * rotateDir);

      place(lv_triangularMask,
            banks->getTriangularBoronMaskPosition(maskId, 0),banks->getTriangularBoronMaskPosition(maskId, 1),banks->getTriangularBoronMaskPosition(maskId, 2),
            lvWorld, BLACK, -5, 0, rotation);
    }
  }

  // Add Calibration slit masks for larmor2022experiment, which is outside of the bank
  const bool withCalibrationSlits = getParameterBoolean("with_calibration_slits");
  if (larmor2022experiment && withCalibrationSlits) {
    const auto calibMask = CalibMasks::getCalibMask("larmorCdCalibMask");
    auto lv_calibrationMask = createCalibrationMaskLV(calibMask);
    auto rotation = new G4RotationMatrix();
    rotation->rotateY(banks->getBankRotation(0, 1));
    rotation->rotateX(banks->getBankRotation(0, 0));
    rotation->rotateZ(banks->getBankRotation(0, 2));

    place(lv_calibrationMask,
          banks->getCalibMaskPositionOutsideBank(calibMask, 0, 0), banks->getCalibMaskPositionOutsideBank(calibMask, 0, 1), banks->getCalibMaskPositionOutsideBank(calibMask, 0, 2),
          lvWorld, PURPLE, -5, 0, rotation);
  }

  delete banks;
  return pvWorld;
}


bool GeoBCS::validateParameters() {
  // you can apply conditions to control the sanity of the geometry parameters and warn the user of possible mistakes
  // a nice example: Projects/SingleCell/G4GeoSingleCell/libsrc/GeoB10SingleCell.cc
  double rear_detector_distance = getParameterDouble("rear_detector_distance_m")*Units::m;
  const bool larmor2022experiment = getParameterBoolean("larmor_2022_experiment");
  if(larmor2022experiment) {
    if (rear_detector_distance != 4.099 *Units::m) {
      printf("ERROR: Wrong rear_detector_distance_m value for the larmor_2022_experiment! (It should be 4.099)\n");
      return false;
    }
  }
  else if(rear_detector_distance < 5.0 *Units::m) {
    printf("ERROR: Wrong rear_detector_distance_m value for LOKI! (It should be >=5.0 m)\n");
      return false;
  }
  return true;
}
