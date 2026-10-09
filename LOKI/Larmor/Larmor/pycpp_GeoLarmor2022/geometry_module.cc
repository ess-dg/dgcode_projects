/////////////////////////////////////////
// Declaration of our geometry module: //
/////////////////////////////////////////

#include "G4Interfaces/GeoConstructPyExport.hh"
#include "G4Para.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include <cmath>
#include <string>
#include <cassert>

#include "G4GeoLoki/BcsBanks.hh"
#include "Larmor/Larmor2022Bank.hh"

// The LOKI rear bank test at the Larmor instrument (ISIS) in 2022: the full LOKI rear bank (28 packs, the bank 0 of
// G4GeoLoki/GeoBCSBanks) at 4.099 m, at the beam height of the experiment. Frozen: the bank placement and the tube
// numbering come from Larmor2022Bank (independent of the LOKI geometry code); the volumes are the same as those of
// the rear bank of G4GeoLoki/GeoBCSBanks (with larmor_2022_experiment=True, from which this was extracted).

class GeoLarmor2022 : public G4Interfaces::GeoConstructBase
{
public:
  GeoLarmor2022();
  virtual ~GeoLarmor2022(){}
  virtual G4VPhysicalVolume* Construct();

protected:
  virtual bool validateParameters();
private:
  //Functions
  G4LogicalVolume * createTubeLV(double converter_thickness, double straw_length);
  G4LogicalVolume * createPackBoxLV(double strawLength, int packNumber);
  G4LogicalVolume * createBankLV(const Larmor2022Bank& larmorBank);
  G4LogicalVolume * createCalibrationMaskLV(CalibMasks::CalibMasksBase calibMask);
};

// this line is necessary to be able to declare the geometry in the python simulation script
PYTHON_MODULE( mod ) { GeoConstructPyExport::exportGeo<GeoLarmor2022>(mod, "GeoLarmor2022"); }

namespace {
  // object rotation of a bank (the columns are the bank-local axes in the world)
  G4RotationMatrix toG4Rotation(const BankTransform& transform) {
    const auto& R = transform.rotation;
    return G4RotationMatrix(G4ThreeVector(R[0][0], R[1][0], R[2][0]),
                            G4ThreeVector(R[0][1], R[1][1], R[2][1]),
                            G4ThreeVector(R[0][2], R[1][2], R[2][2]));
  }
  const int bankId = Larmor2022Bank::bankId;
}

////////////////////////////////////////////
// Implementation of our geometry module: //
////////////////////////////////////////////

GeoLarmor2022::GeoLarmor2022()
  : GeoConstructBase("Larmor/GeoLarmor2022"){
  // declare all parameters that can be used from the command line,
  addParameterDouble("rear_detector_distance_m", 4.099, 4.0, 10.0); // default, min, max (must be 4.099)
  addParameterInt("beamstop_id", 0, 0, 5); // id [1-5] from ESS-1178830 Table 4.6 'Selected beamstop sizes' (0 = no beamstop)
  addParameterBoolean("with_calibration_slits", false); // the Cd calibration slit mask in front of the bank
  addParameterBoolean("old_tube_numbering", false); // the tube numbering of ICD v1 (see Larmor2022Bank)

  addParameterString("world_material","G4_Vacuum");
  addParameterString("B4C_panel_material","MAT_B4C:b10_enrichment=0.95");
}

G4LogicalVolume * GeoLarmor2022::createTubeLV(double converterThickness, double strawLength){
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

///////////  CREATE PACK BOX LOGICAL VOLUME  //////////////////////////
G4LogicalVolume *GeoLarmor2022::createPackBoxLV(double strawLength, int packNumber){
  const double packRotation = BcsBanks::getPackRotation();
  // Instead of a rectangular box, a detector pack is encapsulated in parallelepiped, to avoid collision of the corners with the calibraion slits after applying the pack rotation.
  // The PackBoxWidth corresponds to the size of the volume encapsulating the electronics on the sides as well, not just the detectors, but that would cause collision with the calibration slits, so a multiplication factor of 0.799 is applied, to get a volume just large enough to fit in the detectors in the front.
  auto lv_pack_box = new G4LogicalVolume(
    new G4Para("EmptyPackBox", 0.799*0.5*BcsPack::getPackBoxWidth(), 0.5*BcsPack::getPackBoxHeight(), 0.5 * strawLength + BcsPack::getPackBoxIdleLengthOnOneEnd(), packRotation, 0, 0),
    BcsPack::packBoxFillMaterial, "EmptyPackBox");

  /// Add 8 BCS detector tubes ///
  auto lv_front_tube = createTubeLV(BcsTube::getFrontTubeConverterThickness(), strawLength);
  auto lv_back_tube = createTubeLV(BcsTube::getBackTubeConverterThickness(), strawLength);
  G4RotationMatrix* tubeRotationMatrix = new G4RotationMatrix(0, 0, BcsPack::getTubeRotation());

  const bool oldTubeNumbering = getParameterBoolean("old_tube_numbering");
  for (int inPackTubeId = 0; inPackTubeId < 8; inPackTubeId++) {
    place((inPackTubeId % 4 < 2) ? lv_front_tube : lv_back_tube,
          BcsPack::getHorizontalTubeOffset(inPackTubeId), BcsPack::getVerticalTubeOffset(inPackTubeId), 0,
          lv_pack_box, SILVER, Larmor2022Bank::getTubeVolumeNumber(packNumber, inPackTubeId, oldTubeNumbering), 0, tubeRotationMatrix);
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
G4LogicalVolume *GeoLarmor2022::createCalibrationMaskLV(CalibMasks::CalibMasksBase calibMask){
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
G4LogicalVolume *GeoLarmor2022::createBankLV(const Larmor2022Bank& larmorBank){
  const double strawLength = BcsBanks::getStrawLengthByBankId(bankId);
  const int numberOfPacks = Larmor2022Bank::getNumberOfPacks();
  const double packRotation = BcsBanks::getPackRotation();

  // all positions and sizes in the bank volume are in the bank frame (x = depth, y = across the tubes, z = along
  // the tubes), see G4GeoLoki/BcsBanks.hh
  const auto bankHalfSize = BcsBanks::getBankHalfSizeInBankFrame(bankId);

  auto lv_bank = new G4LogicalVolume(new G4Box("EmptyPanelBox", bankHalfSize[0], bankHalfSize[1], bankHalfSize[2]),
                                     BcsPack::packBoxFillMaterial, "Bank");

  for (int packNumber = 0; packNumber < numberOfPacks; ++packNumber){
    auto lv_pack_box = createPackBoxLV(strawLength, packNumber);
    const auto packPosition = BcsBanks::getPackPositionInBankFrame(bankId, packNumber);
    place(lv_pack_box, packPosition[0], packPosition[1], packPosition[2],
          lv_bank, G4Colour(0, 1, 1), -2, 0, new G4RotationMatrix(0, 0, packRotation));
  }

  const int numberOfBoronMasks = BoronMasks::getNumberOfBoronMasks(bankId);
  for (int maskId = 0; maskId < numberOfBoronMasks; ++maskId){
    const std::string maskName = "BoronMask-"+std::to_string(bankId)+"-"+std::to_string(maskId);
    const auto maskSize = BoronMasks::getSizeInBankFrame(bankId, maskId);
    const auto maskPosition = BcsBanks::getBoronMaskPositionInBankFrame(bankId, maskId);
    place(new G4Box(maskName, 0.5*maskSize[0], 0.5*maskSize[1], 0.5*maskSize[2]),
            BoronMasks::maskMaterial,
            maskPosition[0], maskPosition[1], maskPosition[2],
            lv_bank, BLACK, -2, 0, new G4RotationMatrix(0, 0, BoronMasks::getRotation(bankId, maskId)));
  }

  // Add Beamstop to the Rear Bank volume
  const int beamstopId = getParameterInt("beamstop_id");
  if (beamstopId) { // beamstopId==0 means no beamstop
    const std::string maskName = "BoronMask-Beamstop";
    const double detBankFrontDistance = BcsBanks::detectorSystemFrontDistanceFromBankFront(bankId);
    const double distanceFromDetectorFront = 5*Units::cm;

    const double width = BcsBanks::getBeamstopSize(beamstopId, 0);
    const double height = BcsBanks::getBeamstopSize(beamstopId, 1);
    const double thickness = BcsBanks::getBeamstopSize(beamstopId, 2);

    // on the beam axis (world x = y = 0), 5 cm in front of the detector front: the point of the beam axis at that
    // depth, in the bank frame (this compensates the bank elevation)
    const double depth = -bankHalfSize[0] + detBankFrontDistance - distanceFromDetectorFront;
    const BankTransform transform = larmorBank.getBankTransform();
    const auto& R = transform.rotation;
    const auto& t = transform.translation;
    const double z = (depth + R[0][0]*t[0] + R[1][0]*t[1] + R[2][0]*t[2]) / R[2][0]; // n.(P - t) = depth
    const auto local = transform.toLocal({0.0, 0.0, z});
    const G4ThreeVector position(local[0], local[1], local[2]);
    place(new G4Box(maskName, 0.5* thickness, 0.5* height, 0.5* width),
          BoronMasks::maskMaterial,
          position.x(), position.y(), position.z(),
          lv_bank, BLACK, -5, 0, new G4RotationMatrix());
  }

  return lv_bank;
 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* GeoLarmor2022::Construct(){
  // this is where we put the entire geometry together, the private functions creating the logical volumes are meant to facilitate the code below
  const double rear_detector_distance = getParameterDouble("rear_detector_distance_m")*Units::m;
  const Larmor2022Bank larmorBank(rear_detector_distance);

  // calculate a value that is big enough to fit your world volume, the "super mother"
  double big_dimension = 1.1*( 1 *Units::m + rear_detector_distance);

  //World volume:
  auto world_material = getParameterMaterial("world_material");
  auto worldvols = place(new G4Box("World", big_dimension, big_dimension, big_dimension), world_material, 0, 0, 0, 0, INVISIBLE);
  auto lvWorld = worldvols.logvol;
  auto pvWorld = worldvols.physvol;

  // Create and place the detector bank
  auto lv_bank = createBankLV(larmorBank);
  const BankTransform transform = larmorBank.getBankTransform();
  // G4PVPlacement takes the frame rotation, i.e. the inverse of the bank rotation
  auto bankRotation = new G4RotationMatrix(toG4Rotation(transform).inverse());
  place(lv_bank, transform.translation[0], transform.translation[1], transform.translation[2], lvWorld, ORANGE, bankId, 0, bankRotation);

  // Add the calibration slit mask, which is outside of the bank
  if (getParameterBoolean("with_calibration_slits")) {
    const auto calibMask = Larmor2022Bank::getCalibMask();
    auto lv_calibrationMask = createCalibrationMaskLV(calibMask);
    auto rotation = new G4RotationMatrix();
    rotation->rotateY(0.5*M_PI);
    rotation->rotateX(0.0);
    rotation->rotateZ(0.0);

    const auto position = larmorBank.getCalibMaskPosition();
    place(lv_calibrationMask, position[0], position[1], position[2], lvWorld, PURPLE, -5, 0, rotation);
  }

  return pvWorld;
}


bool GeoLarmor2022::validateParameters() {
  // you can apply conditions to control the sanity of the geometry parameters and warn the user of possible mistakes
  // a nice example: Projects/SingleCell/G4GeoSingleCell/libsrc/GeoB10SingleCell.cc
  if (getParameterDouble("rear_detector_distance_m")*Units::m != 4.099 *Units::m) {
    printf("ERROR: Wrong rear_detector_distance_m value for the Larmor 2022 experiment! (It should be 4.099)\n");
    return false;
  }
  return true;
}
