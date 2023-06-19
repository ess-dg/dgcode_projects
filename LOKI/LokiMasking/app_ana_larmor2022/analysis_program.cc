#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Core/FPE.hh"
#include "Utils/ArrayMath.hh"
#include "Utils/NeutronMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <string>
#include <cmath>
#include <iostream>
#include "G4GeoLoki/PixelatedBanks.hh"
#include "LokiMasking/MaskFileCreator.hh"

//Griff analysis. See https://confluence.esss.lu.se/display/DGCODE/Griff for more info.

#ifndef M_PI
#define M_PI 3.14159265358979323846 //  pi
#endif


int main(int argc, char **argv) {
  Core::catch_fpe();
  GriffDataReader dr(argc, argv);

  auto setup = dr.setup();
  auto &geo = setup->geo();
  auto &gen = setup->gen();

  if (geo.getName() != "G4GeoLoki/GeoBCSBanks" && geo.getName() != "G4GeoBCS/GeoLarmorBCSExperiment") {
    printf("Error: Wrong setup for this analysis\n");
    return 1;
  }

  setup->dump();

  GriffAnaUtils::TrackIterator primary_geantinos(&dr);
  primary_geantinos.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  primary_geantinos.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(999));

  GriffAnaUtils::SegmentIterator segments_bank(&dr);
  segments_bank.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Bank"));

  SimpleHists::HistCollection hc;

  const double rearDetectorDistance = setup->geo().getParameterDouble("rear_detector_distance_m") * Units::m;

  auto userData = setup->userData();
  int strawPixelNumber = 512; //hardcoded default
  if(userData.count("analysis_straw_pixel_number")){
    strawPixelNumber = std::stoi(userData["analysis_straw_pixel_number"].c_str());
  }

  const int numberOfPixels = 1568 * strawPixelNumber; //fixed 1568 straws in the geometry; (802816/1605632 for 512/1024)
  
  printf("Rear bank pixel number for analysis: %d in total, %d per straw, \n", numberOfPixels, strawPixelNumber);

  PixelatedBanks banks = PixelatedBanks(rearDetectorDistance, strawPixelNumber);

  auto h_geantino_pixel_enter = hc.book2D("Shows pixels the geantinos entered", strawPixelNumber, 0, strawPixelNumber, numberOfPixels / strawPixelNumber, 0, numberOfPixels / strawPixelNumber, "h_geantino_pixel_enter");
  h_geantino_pixel_enter->setXLabel("Pixel ID along straw");
  h_geantino_pixel_enter->setYLabel("Straw ID");

  auto h_geantino_pixel_enter_masked = hc.book2D("Shows pixels the geantinos entered without transmission through masks", strawPixelNumber, 0, strawPixelNumber, numberOfPixels / strawPixelNumber, 0, numberOfPixels / strawPixelNumber, "h_geantino_pixel_enter_masked");
  h_geantino_pixel_enter_masked->setXLabel("Pixel ID along straw");
  h_geantino_pixel_enter_masked->setYLabel("Straw ID");

  auto h_geantino_xy_bank = hc.book2D("Geantino xy (at bank entering)", 2500, -1250, 1250, 2500, -1250, 1250, "geantino_xy_bank");
       h_geantino_xy_bank->setXLabel("-x [mm]");
       h_geantino_xy_bank->setYLabel("y [mm]");

  auto h_geantino_xy_bank_filtered = hc.book2D("Geantino xy filtered (at bank entering)", 2500, -1250, 1250, 2500, -1250, 1250, "geantino_xy_bank_filtered");
       h_geantino_xy_bank_filtered->setXLabel("-x [mm]");
       h_geantino_xy_bank_filtered->setYLabel("y [mm]");

  auto h_counters = hc.bookCounts("General geantino counters", "geantino_counters");
  auto countTestGeantino = h_counters->addCounter("all_geantino");
  auto countTestGeantinoAbsInMask = h_counters->addCounter("geantino_in_Mask");

  const int indexOffset = 11;
  MaskFileCreator masking("maskFile.xml", numberOfPixels, indexOffset);

  ///const double xWidthVacuumTankEnd = 687.4 *Units::mm;
  ///const double yHeightVacuumTankEnd = 600 *Units::mm;
  ///const double zVacuumTankEnd = 4402 *Units::mm;
  ///const double zBankFront = rearDetectorDistance - banks.detectorSystemFrontDistanceFromBankFront(0) *Units::mm;
  ///const double xBankEnterLimit = (xWidthVacuumTankEnd * 0.5) * (zBankFront / zVacuumTankEnd);
  //const double yBankEnterLimit = (yHeightVacuumTankEnd * 0.5) * (zBankFront / zVacuumTankEnd);
  ////std::cout<<"\n ***** \n xBankEnterLimit: " << xBankEnterLimit << "\n yBankEnterLimit: "<< yBankEnterLimit<< "\n ****** \n";

  const double zVacuumTankEnd = rearDetectorDistance - 90*Units::mm; // "The front of the Loki detector was 90 mm in front of the tank”
  const double zBankFront = rearDetectorDistance - banks.detectorSystemFrontDistanceFromBankFront(0) *Units::mm;

  /////////////////////////////////////////////////////////////////////
  //const double xWidthVacuumTankEnd = 790 *Units::mm;  //overestimated - vacuumtank end
  //const double yHeightVacuumTankEnd = 700 *Units::mm; //overestimated - vacuumtank end
  const double xWidthEmpiricalDetFront = 2*0.27 *Units::m; //illuminated area based on measured data
  const double yHeightEmpiricalDetFront = 2*0.31 *Units::m; //illuminated area based on measured data
  const double xWidthVacuumTankEnd = xWidthEmpiricalDetFront * (zVacuumTankEnd / rearDetectorDistance);
  const double yHeightVacuumTankEnd = yHeightEmpiricalDetFront * (zVacuumTankEnd / rearDetectorDistance);

  const double xBankEnterLimit = (xWidthVacuumTankEnd * 0.5) * (zBankFront / zVacuumTankEnd);
  const double yBankEnterLimit = (yHeightVacuumTankEnd * 0.5) * (zBankFront / zVacuumTankEnd);
  const double yBeamOffset = 0.008 *Units::m; //empirical
  //std::cout<<"\n ***** \n xBankEnterLimit: " << xBankEnterLimit << "\n yBankEnterLimit: "<< yBankEnterLimit<< "\n ****** \n";

  double tmp_x_offset_meters = gen.hasParameterDouble("x_offset_meters") ? gen.getParameterDouble("x_offset_meters") *Units::m : 0.0;
  double tmp_dx_meter = gen.hasParameterDouble("dx_meter") ? gen.getParameterDouble("dx_meter") *Units::m : 0.0;
  if(tmp_x_offset_meters && tmp_dx_meter){
    printf("Error: both 'x_offset_meters' and 'dx_meter' parameters are defined.\n");
    return 1;
  }
  const double xBeamOffset = tmp_x_offset_meters ? tmp_x_offset_meters :tmp_dx_meter;

  while (dr.loopEvents()) {
    while (auto trk_geantino = primary_geantinos.next()) {
      countTestGeantino += 1;

      if(auto bankSegment = segments_bank.next()){
        auto stepFirstInBank = bankSegment->firstStep();
        const double xBankEnter = stepFirstInBank->preGlobalX()/Units::mm;
        const double yBankEnter = stepFirstInBank->preGlobalY()/Units::mm;
        h_geantino_xy_bank->fill(-xBankEnter, yBankEnter, 1);

        //if ( yBankEnterLimit < std::fabs(yBankEnter) ||
        if ( yBankEnter > yBeamOffset + yBankEnterLimit || yBeamOffset - yBankEnterLimit > yBankEnter ||
          xBankEnter > xBeamOffset + xBankEnterLimit || xBeamOffset - xBankEnterLimit > xBankEnter ){
          continue;
        }

        h_geantino_xy_bank_filtered->fill(-xBankEnter, yBankEnter, 1);
      }

      bool geantinoAbsorbed = false;
      for (auto seg = trk_geantino->segmentBegin(); seg != trk_geantino->segmentEnd(); ++seg) {

        if (!geantinoAbsorbed && (seg->volumeName().find("BoronMask-") != std::string::npos || seg->volumeName() == "B4CPanel" || seg->volumeName() == "AlPanel")) {
          countTestGeantinoAbsInMask += 1;
          geantinoAbsorbed = true;
          //break;
        }
        else if (seg->volumeName() == "Converter") {
          const int strawId_conv = seg->volumeCopyNumber(1);
          const int tubeId_conv = seg->volumeCopyNumber(3);
          const int bankId_conv = seg->volumeCopyNumber(5);

          auto step = seg->lastStep();
          const int pixelId = banks.getPixelId(bankId_conv, tubeId_conv, strawId_conv, step->postGlobalX(), step->postGlobalY());

          if (!geantinoAbsorbed && !masking.isPixelEntered(pixelId)) {
            h_geantino_pixel_enter_masked->fill(pixelId % strawPixelNumber, std::floor(pixelId / strawPixelNumber), 1);
            masking.setPixelEntered(pixelId);
          }
          if (!masking.isPixelEnteredAimingCheck(pixelId)) { // check if all pixels are aimed at
            h_geantino_pixel_enter->fill(pixelId % strawPixelNumber, std::floor(pixelId / strawPixelNumber), 1);
            masking.setPixelEnteredAimingCheck(pixelId);
          }
        }
      }
    }
  }   //end of event loop

  masking.createMaskFile();
  hc.saveToFile("bcsloki_masking", true);

  return 0;
}
