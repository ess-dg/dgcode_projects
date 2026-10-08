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
#include "G4AffineTransform.hh"
#include <algorithm>
#include <stdexcept>
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
  G4LogicalVolume * createPackBoxLV(int bankId, int packNumber);
  G4LogicalVolume * createBankLV(int bankId);
  G4LogicalVolume * createTriangularMaskLV(int maskId);
  G4LogicalVolume * createCalibrationMaskLV(CalibMasks::CalibMasksBase calibMask);
};

// this line is necessary to be able to declare the geometry in the python simulation script
PYTHON_MODULE( mod ) { GeoConstructPyExport::exportGeo<GeoBCS>(mod, "GeoBCSBanks"); }

namespace {
  // object rotation of a bank (the columns are the bank-local axes in the world)
  G4RotationMatrix toG4Rotation(const BankTransform& transform) {
    const auto& R = transform.rotation;
    return G4RotationMatrix(G4ThreeVector(R[0][0], R[1][0], R[2][0]),
                            G4ThreeVector(R[0][1], R[1][1], R[2][1]),
                            G4ThreeVector(R[0][2], R[1][2], R[2][2]));
  }
}

////////////////////////////////////////////
// Implementation of our geometry module: //
////////////////////////////////////////////

GeoBCS::GeoBCS()
  : GeoConstructBase("G4GeoLoki/GeoBCSBanks"){
  // declare all parameters that can be used from the command line,
  addParameterDouble("rear_detector_distance_m", 5.0, 4.0, 10.0); // default, min, max
  addParameterInt("beamstop_id", 0, 0, 5); // id [1-5] from ESS-1178830 Table 4.6 'Selected beamstop sizes' (0 = no beamstop)
  addParameterBoolean("with_calibration_slits", false);

  // bank placements: "nominal-geant4-geometry" or the name of a calibration in G4GeoLoki/data/bank_calibration_<name>.txt
  // (or the path of such a file), see G4GeoLoki/BankCalibration.hh
  addParameterString("bank_calibration", BankCalibration::defaultName);

  addParameterString("world_material","G4_Vacuum");
  addParameterString("B4C_panel_material","MAT_B4C:b10_enrichment=0.95");
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

///////////  CREATE PACK BOX LOGICAL VOLUME  //////////////////////////
G4LogicalVolume *GeoBCS::createPackBoxLV(int bankId, int packNumber){
  const double strawLength = banks->getStrawLengthByBankId(bankId);
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
          lv_pack_box, SILVER, banks->getTubeIdInBank(bankId, packNumber, inPackTubeId), 0, tubeRotationMatrix);
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
  const int numberOfPacks = banks->getNumberOfPacksByBankId(bankId);

  const double packRotation = banks->getPackRotation();

  // all positions and sizes in the bank volume are in the bank frame (x = depth, y = across the tubes, z = along
  // the tubes), see BcsBanks.hh
  const auto bankHalfSize = banks->getBankHalfSizeInBankFrame(bankId);

  auto lv_bank = new G4LogicalVolume(new G4Box("EmptyPanelBox", bankHalfSize[0], bankHalfSize[1], bankHalfSize[2]),
                                     BcsPack::packBoxFillMaterial, "Bank");

  // override lv_bank to subtract some empty part of front left and right banks where front top and bottom banks would overlap with them
  // (the notch is computed from the bank placements in use, see BcsBanks::getBankNotch)
  if (banks->hasBankNotch(bankId)) {
    const auto notch = banks->getBankNotch(bankId);
    if (!notch.isNominal)
      printf("GeoBCSBanks: bank %d: notch enlarged for bank_calibration=%s to the half sizes (%g, %g, %g) mm at (%g, %g, %g) mm (nominal: (15, 70, 60) mm at (-145, -42, -425) mm)\n",
             bankId, banks->getBankCalibration().name().c_str(), notch.halfSize[0], notch.halfSize[1], notch.halfSize[2], notch.centre[0], notch.centre[1], notch.centre[2]);
    auto fullBankBox = new G4Box("EmptyPanelBox", bankHalfSize[0], bankHalfSize[1], bankHalfSize[2]);
    auto bankBoxCut = new G4Box("EmptyPanelBox", notch.halfSize[0], notch.halfSize[1], notch.halfSize[2]);
    auto bankBox = new G4SubtractionSolid("EmptyPanelBox", fullBankBox, bankBoxCut, 0, G4ThreeVector(notch.centre[0], notch.centre[1], notch.centre[2]));

    lv_bank = new G4LogicalVolume(bankBox, BcsPack::packBoxFillMaterial, "Bank");
  }

  for (int packNumber = 0; packNumber < numberOfPacks; ++packNumber){
    auto lv_pack_box = createPackBoxLV(bankId, packNumber);
    const auto packPosition = banks->getPackPositionInBankFrame(bankId, packNumber);
    place(lv_pack_box, packPosition[0], packPosition[1], packPosition[2],
          lv_bank, G4Colour(0, 1, 1), -2, 0, new G4RotationMatrix(0, 0, packRotation));
  }

  const int numberOfBoronMasks = BoronMasks::getNumberOfBoronMasks(bankId);
  for (int maskId = 0; maskId < numberOfBoronMasks; ++maskId){
    const std::string maskName = "BoronMask-"+std::to_string(bankId)+"-"+std::to_string(maskId);
    const auto maskSize = BoronMasks::getSizeInBankFrame(bankId, maskId);
    const auto maskPosition = banks->getBoronMaskPositionInBankFrame(bankId, maskId);
    place(new G4Box(maskName, 0.5*maskSize[0], 0.5*maskSize[1], 0.5*maskSize[2]),
            BoronMasks::maskMaterial,
            maskPosition[0], maskPosition[1], maskPosition[2],
            lv_bank, BLACK, -2, 0, new G4RotationMatrix(0, 0, BoronMasks::getRotation(bankId, maskId)));
  }

  // Add Beamstop to the Rear Bank volume
  const int beamstopId = getParameterInt("beamstop_id");
  if (bankId == 0 && beamstopId) { // beamstopId==0 means no beamstop
    const std::string maskName = "BoronMask-Beamstop";
    const double detBankFrontDistance = banks->detectorSystemFrontDistanceFromBankFront(bankId);
    const double distanceFromDetectorFront = 5*Units::cm;

    const double width = banks->getBeamstopSize(beamstopId, 0);
    const double height = banks->getBeamstopSize(beamstopId, 1);
    const double thickness = banks->getBeamstopSize(beamstopId, 2);

    // on the beam axis (world x = y = 0), 5 cm in front of the detector front: the point of the beam axis at that
    // depth, in the bank frame (this compensates the bank elevation, also of a calibrated bank)
    const double depth = -bankHalfSize[0] + detBankFrontDistance - distanceFromDetectorFront;
    const BankTransform transform = banks->getBankTransform(bankId);
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

  const bool withCalibrationSlits = getParameterBoolean("with_calibration_slits");
  if (withCalibrationSlits) {
    std::string calibMaskName = "lokiStandard-"+std::to_string(bankId);
    const auto calibMask = CalibMasks::getCalibMask(calibMaskName);
    auto lv_calibrationMask = createCalibrationMaskLV(calibMask);

    const auto calibMaskPosition = banks->getCalibMaskPositionInBankFrame(calibMask, bankId);
    place(lv_calibrationMask, calibMaskPosition[0], calibMaskPosition[1], calibMaskPosition[2],
          lv_bank, PURPLE, -5, 0, new G4RotationMatrix());
  }

  // An enlarged notch must not cut into the content of the bank (packs, masks): it is only meant to remove empty
  // space of the bank volume. (Conservative check with the bounding boxes of the daughters in the bank frame.)
  if (banks->hasBankNotch(bankId) && !banks->getBankNotch(bankId).isNominal) {
    const auto notch = banks->getBankNotch(bankId);
    for (size_t i = 0; i < lv_bank->GetNoDaughters(); i++) {
      const auto daughter = lv_bank->GetDaughter(i);
      G4ThreeVector localMin, localMax;
      daughter->GetLogicalVolume()->GetSolid()->BoundingLimits(localMin, localMax);
      const G4RotationMatrix* frameRotation = daughter->GetRotation();
      const G4AffineTransform toBank = frameRotation ? G4AffineTransform(frameRotation, daughter->GetTranslation()) : G4AffineTransform(daughter->GetTranslation());
      G4ThreeVector low(1e99, 1e99, 1e99), high(-1e99, -1e99, -1e99);
      for (int corner = 0; corner < 8; corner++) {
        const G4ThreeVector p = toBank.TransformPoint(G4ThreeVector(corner & 1 ? localMax.x() : localMin.x(),
                                                                    corner & 2 ? localMax.y() : localMin.y(),
                                                                    corner & 4 ? localMax.z() : localMin.z()));
        for (int k = 0; k < 3; k++) {
          low[k] = std::min(low[k], p[k]);
          high[k] = std::max(high[k], p[k]);
        }
      }
      bool apart = false;
      for (int k = 0; k < 3; k++)
        if (high[k] <= notch.centre[k] - notch.halfSize[k] || low[k] >= notch.centre[k] + notch.halfSize[k])
          apart = true;
      if (!apart)
        throw std::runtime_error("GeoBCSBanks: bank_calibration=" + banks->getBankCalibration().name() + ": the notch of bank "
                                 + std::to_string(bankId) + ", enlarged because bank " + std::to_string(banks->getBankNotchNeighbour(bankId))
                                 + " overlaps it, would cut into " + std::string(daughter->GetName())
                                 + ": the two banks collide (check the calibration)");
    }
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
  banks = new BcsBanks(rear_detector_distance, 9, getParameterString("bank_calibration"));

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

    // bank placement (the same transform as used by AimHelper and PixelatedBanks)
    const BankTransform transform = banks->getBankTransform(bankId);
    // G4PVPlacement takes the frame rotation, i.e. the inverse of the bank rotation
    auto rotation = new G4RotationMatrix(toG4Rotation(transform).inverse());

    place(lv_bank, transform.translation[0], transform.translation[1], transform.translation[2], lvWorld, ORANGE, bankId, 0, rotation);
  }

  // Add 4 triangular boron masks (added to the World instead of the banks)
  for (int maskId = 0; maskId <= 3; maskId++) {
    auto lv_triangularMask = createTriangularMaskLV(maskId);

    // mounted on the front of bank 5 or 7: its placement in the bank frame, moved with the bank transform
    const int bankId = BoronMasks::getBankIdOfTriangularMask(maskId);
    const auto inBank = banks->getTriangularBoronMaskPlacementInBankFrame(maskId);
    const BankTransform bank = banks->getBankTransform(bankId);
    const auto position = bank.toGlobal(inBank.position);
    // G4PVPlacement takes the frame rotation: the inverse of the object rotation (bank rotation * rotation in bank)
    const G4RotationMatrix rotationInBank(G4ThreeVector(inBank.rotation[0][0], inBank.rotation[1][0], inBank.rotation[2][0]),
                                          G4ThreeVector(inBank.rotation[0][1], inBank.rotation[1][1], inBank.rotation[2][1]),
                                          G4ThreeVector(inBank.rotation[0][2], inBank.rotation[1][2], inBank.rotation[2][2]));
    auto rotation = new G4RotationMatrix((toG4Rotation(bank) * rotationInBank).inverse());

    place(lv_triangularMask, position[0], position[1], position[2], lvWorld, BLACK, -5, 0, rotation);
  }

  delete banks;
  return pvWorld;
}


bool GeoBCS::validateParameters() {
  // you can apply conditions to control the sanity of the geometry parameters and warn the user of possible mistakes
  // a nice example: Projects/SingleCell/G4GeoSingleCell/libsrc/GeoB10SingleCell.cc
  double rear_detector_distance = getParameterDouble("rear_detector_distance_m")*Units::m;
  const std::string bankCalibration = getParameterString("bank_calibration");
  try {
    BankCalibration::load(bankCalibration);
  }
  catch (const std::exception& error) {
    printf("ERROR: Wrong bank_calibration value: %s\n", error.what());
    return false;
  }
  if(rear_detector_distance < 5.0 *Units::m) {
    printf("ERROR: Wrong rear_detector_distance_m value for LOKI! (It should be >=5.0 m)\n");
      return false;
  }
  return true;
}
