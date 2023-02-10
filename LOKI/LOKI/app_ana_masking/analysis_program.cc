#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Core/FPE.hh"
#include "Utils/ArrayMath.hh"
#include "Utils/NeutronMath.hh"
#include "SimpleHists/HistCollection.hh"
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include "LOKI/MaskFileCreator.hh"
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
  const double sampleDetectorDistance = setup->geo().getParameterDouble("rear_detector_distance_m") *Units::m;
  int strawPixelNumber = 0;
  if(userData.count("analysis_straw_pixel_number")){
    strawPixelNumber = std::stoi(userData["analysis_straw_pixel_number"].c_str());
    banks = new PixelatedBanks(sampleDetectorDistance, strawPixelNumber);
  }
  else{ // use default rear bank pixel number
    banks = new PixelatedBanks(sampleDetectorDistance);
    strawPixelNumber = banks->getNumberOfPixelsInStraw(0);//NOTE: assuming same number of pixels for each bank
  }

  const int numberOfPixels = banks->getTotalNumberOfPixels();

  auto h_neutron_pixel_geantino = hc.book2D("Show pixels the geantinos entered", strawPixelNumber, 0, strawPixelNumber, numberOfPixels / strawPixelNumber, 0, numberOfPixels / strawPixelNumber, "h_neutron_pixel_geantino");
  h_neutron_pixel_geantino->setXLabel("Pixel ID along straw");
  h_neutron_pixel_geantino->setYLabel("Straw ID");

  auto h_neutron_pixel_geantino_masking = hc.book2D("Show pixels the geantinos entered (no mask transmission)", strawPixelNumber, 0, strawPixelNumber, numberOfPixels / strawPixelNumber, 0, numberOfPixels / strawPixelNumber, "h_neutron_pixel_geantino_masking");
  h_neutron_pixel_geantino_masking->setXLabel("Pixel ID along straw");
  h_neutron_pixel_geantino_masking->setYLabel("Straw ID");

  auto h_neutron_counters = hc.bookCounts("General neutron counters", "neutron_counters");
  auto countTestGeantino = h_neutron_counters->addCounter("all_geantino");
  auto countTestGeantinoAbsInMask = h_neutron_counters->addCounter("geantino_in_Mask");

  const int indexOffset = 11;
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
            h_neutron_pixel_geantino_masking->fill(pixelId % strawPixelNumber, std::floor(pixelId / strawPixelNumber), 1);
            masking.setPixelEntered(pixelId);
          }
          if (!masking.isPixelEnteredAimingCheck(pixelId)) { // check if all pixels are aimed at
            h_neutron_pixel_geantino->fill(pixelId % strawPixelNumber, std::floor(pixelId / strawPixelNumber), 1);
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
