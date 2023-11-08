#include "GriffAnaUtils/All.hh"

#include "Core/FPE.hh"
#include "Utils/ArrayMath.hh"
#include "Utils/NeutronMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include "LokiMasking/MaskFileCreator.hh"
#include "G4GeoLoki/PixelatedBanks.hh"
//Griff analysis. See https://confluence.esss.lu.se/display/DGCODE/Griff for more info.

#ifndef M_PI
#define M_PI 3.14159265358979323846 //  pi
#endif


int main(int argc, char **argv) {
  Core::catch_fpe();
  GriffDataReader dr(argc, argv);

  auto setup = dr.setup();
  auto &geo = setup->geo();

  if (geo.getName() != "G4GeoLoki/GeoBCSBanks" && geo.getName() != "G4GeoBCS/GeoLarmorBCSExperiment") {
    printf("Error: Wrong setup for this analysis\n");
    return 1;
  }

  setup->dump();

  GriffAnaUtils::TrackIterator primary_geantinos(&dr);
  primary_geantinos.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  primary_geantinos.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(999));

  SimpleHists::HistCollection hc;

  auto userData = setup->userData();
  PixelatedBanks* banks;
  const double rearDetectorDistance = setup->geo().getParameterDouble("rear_detector_distance_m") *Units::m;
  int strawPixelNumber = 0;
  if(userData.count("analysis_straw_pixel_number")){
    strawPixelNumber = std::stoi(userData["analysis_straw_pixel_number"].c_str());
    banks = new PixelatedBanks(rearDetectorDistance, strawPixelNumber);
  }
  else{ // use default rear bank pixel number
    banks = new PixelatedBanks(rearDetectorDistance);
    strawPixelNumber = banks->getNumberOfPixelsInStraw(0);//NOTE: assuming same number of pixels for each bank
  }

  const int numberOfPixels = banks->getTotalNumberOfPixels();

  auto h_geantino_pixel_enter = hc.book2D("Shows pixels the geantinos entered", strawPixelNumber, 0, strawPixelNumber, numberOfPixels / strawPixelNumber, 0, numberOfPixels / strawPixelNumber, "h_geantino_pixel_enter");
  h_geantino_pixel_enter->setXLabel("Pixel ID along straw");
  h_geantino_pixel_enter->setYLabel("Straw ID");

  auto h_geantino_pixel_enter_masked = hc.book2D("Shows pixels the geantinos entered without transmission through masks", strawPixelNumber, 0, strawPixelNumber, numberOfPixels / strawPixelNumber, 0, numberOfPixels / strawPixelNumber, "h_geantino_pixel_enter_masked");
  h_geantino_pixel_enter_masked->setXLabel("Pixel ID along straw");
  h_geantino_pixel_enter_masked->setYLabel("Straw ID");

  auto h_counters = hc.bookCounts("General geantino counters", "geantino_counters");
  auto countTestGeantino = h_counters->addCounter("all_geantino");
  auto countTestGeantinoAbsInMask = h_counters->addCounter("geantino_in_Mask");

  const int indexOffset = 1; //detector IDs in the IDF (and ICD) file starts from 1 (as opposed to the zero-based numbering in the Geant4 geometry)
  MaskFileCreator masking("maskFile.xml", numberOfPixels, indexOffset);

  while (dr.loopEvents()) {
    while (auto trk_geantino = primary_geantinos.next()) {
      countTestGeantino += 1;

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
          const int pixelId = banks->getPixelId(bankId_conv, tubeId_conv, strawId_conv, step->postGlobalX(), step->postGlobalY());

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
