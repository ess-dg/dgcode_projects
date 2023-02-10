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

//Griff analysis. See https://confluence.esss.lu.se/display/DGCODE/Griff for more info.

#ifndef M_PI
#define M_PI 3.14159265358979323846 //  pi
#endif

const int strawPixelNumber = 512;

int getPositionPixelId(const int tubeId, const double positionX){
  const double strawLength = tubeId < 16 ? 1.5*Units::m : 1.2 *Units::m; //0-15 1.5m ; 16-31 1.2 m
  const double pixelLength = strawLength / strawPixelNumber;

  const double strawBegin = - 0.5* strawLength;
  const int invertedPixelId = std::floor((positionX - strawBegin) / pixelLength);

  return (strawPixelNumber - 1) - invertedPixelId; //pixels are numbered in minus x direction
}

int getPixelId(const int tubeId, const int strawId, const double positionX) {
  const int strawPixelOffset = (tubeId * 7 + strawId) * strawPixelNumber;
  const int positionPixelId = getPositionPixelId(tubeId, positionX);
  return strawPixelOffset + positionPixelId;
}

int main(int argc, char **argv) {
  Core::catch_fpe();
  GriffDataReader dr(argc, argv);

  auto setup = dr.setup();
  auto &geo = setup->geo();

  if (geo.getName() != "G4GeoLoki/GeoBCSBanks" && geo.getName() != "G4GeoLoki/GeoLarmorBCSExperiment") {
    printf("Error: Wrong setup for this analysis\n");
    return 1;
  }

  setup->dump();

  GriffAnaUtils::TrackIterator primary_geantinos(&dr);
  primary_geantinos.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  primary_geantinos.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(999));

  SimpleHists::HistCollection hc;

  //const double sampleDetectorDistance = setup->geo().getParameterDouble("rear_detector_distance_m") * Units::m;

  const int numberOfPixels = 4 * 8 * 7 * strawPixelNumber; //hardcoded: 4(pack) * 8(tube) * 7(straw) * 512(pixel)

  auto h_neutron_pixel_geantino = hc.book2D("Show pixels the geantinos entered", strawPixelNumber, 0, strawPixelNumber, numberOfPixels / strawPixelNumber, 0, numberOfPixels / strawPixelNumber, "h_neutron_pixel_geantino");
  h_neutron_pixel_geantino->setXLabel("Pixel ID along straw");
  h_neutron_pixel_geantino->setYLabel("Straw ID");

  auto h_neutron_pixel_geantino_masking = hc.book2D("Show pixels the geantinos entered (no mask transmission)", strawPixelNumber, 0, strawPixelNumber, numberOfPixels / strawPixelNumber, 0, numberOfPixels / strawPixelNumber, "h_neutron_pixel_geantino_masking");
  h_neutron_pixel_geantino_masking->setXLabel("Pixel ID along straw");
  h_neutron_pixel_geantino_masking->setYLabel("Straw ID");

  auto h_counters = hc.bookCounts("General neutron counters", "neutron_counters");
  auto countTestGeantino = h_counters->addCounter("all_geantino");
  auto countTestGeantinoAbsInMask = h_counters->addCounter("geantino_in_Mask");

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

          auto step = seg->lastStep();
          const int pixelId = getPixelId(tubeId_conv, strawId_conv, step->postGlobalX());

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
