#include "GriffAnaUtils/All.hh"

#include "Core/FPE.hh"
#include "Utils/ArrayMath.hh"
#include "Utils/NeutronMath.hh"
#include "SimpleHists/HistCollection.hh"
#include "GriffB10Common/DetHitApproximation.hh"
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include "mcpl.h"

// #include "G4GeoLoki/BcsBanks.hh"
//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griff for more info.

#ifndef M_PI
#define M_PI  3.14159265358979323846  //  pi
#endif

const int IDFdetectorPixelOffset = 11;
const int strawPixelNumber = 512;
//const int strawNumber = 4 * 8 * 7;

int getPositionPixelId(int tubeId, double positionX){
  const double strawLength = tubeId < 16 ? 1.5*Units::m : 1.2 *Units::m; //0-15 1.5m ; 16-31 1.2 m

  const double pixelLength = strawLength / strawPixelNumber;

  const double strawBegin = - 0.5* strawLength;
  const int invertedPixelId = std::floor((positionX - strawBegin) / pixelLength);

  return (strawPixelNumber - 1) - invertedPixelId; //pixels are numbered in minus x direction
}

int getPixelId(int tubeId, int strawId, double positionX) {
  const int strawPixelOffset = (tubeId * 7 + strawId) * strawPixelNumber;
  const int positionPixelId = getPositionPixelId(tubeId, positionX);
  return strawPixelOffset + positionPixelId;
}

int main(int argc, char**argv) {

  Core::catch_fpe();
  GriffDataReader dr(argc,argv);

  auto setup = dr.setup();
  auto &geo = setup->geo();
  //printf("QQQ=============  %s \n", geo.getName().c_str());
  if (geo.getName()!="G4GeoLoki/GeoLarmorBCSExperiment" && geo.getName()!="G4GeoBCS/GeoLarmorBCSExperiment") {
    printf("Error: Wrong setup for this analysis\n");
    return 1;
  }

  auto &gen = setup->gen();

  double sourceSampleDistance = 0*Units::m;
  if (gen.getName()=="G4MCPLPlugins/MCPLGen") {
    sourceSampleDistance = 25.3*Units::m; //From Larmor McStas file
  }
  printf("sourceSampleDistance: %f\n", sourceSampleDistance);

  setup->dump();

  DetHitApproximation hit(&dr,1.2*Units::cm,120*Units::keV,"CountingGas" ); // defined for He3 gas

  GriffAnaUtils::TrackIterator primary_neutrons(&dr);
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));


  mcpl_outfile_t detMcpl;
  mcpl_particle_t *mcplParticle;

  detMcpl = mcpl_create_outfile("detectionEvents_minimal.mcpl");
  mcpl_hdr_add_comment(detMcpl, "Neutrons in this file are actually detection events and userflags are pixel ID's of hits. Created with sb_loki_loki_ana command.");
  mcpl_enable_userflags(detMcpl);
  mcplParticle = mcpl_get_empty_particle(detMcpl);

  SimpleHists::HistCollection hc;

  auto h_neutron_xy_hit = hc.book2D("Neutron xy (hit)", 2500, -1250, 1250, 2500, -1250, 1250, "neutron_xy_hit");
       h_neutron_xy_hit->setXLabel("-x [mm]");
       h_neutron_xy_hit->setYLabel("y [mm]");

  auto h_neutron_Q_hit = hc.book1D("Neutron Q [1/angstrom] (hit)", 250, 0, 2, "neutron_Q_hit");
       h_neutron_Q_hit->setXLabel("Q [1/angstrom]");


  auto h_neutron_counters = hc.bookCounts("General neutron counters","neutron_counters"); /////////////
  auto count_initial_neutrons = h_neutron_counters->addCounter("count_initial_neutrons");
  auto count_neutrons_hit = h_neutron_counters->addCounter("count_neutrons_hit");

/*
  auto count_pixel_gas_vs_det_diff = h_neutron_counters->addCounter("count_pixel_gas_vs_det_diff");
  auto count_pixel_gas_vs_det_diff_more = h_neutron_counters->addCounter("count_pixel_gas_vs_det_diff_more");
  auto count_pixel_conv_vs_det_diff = h_neutron_counters->addCounter("count_pixel_conv_vs_det_diff");
  auto count_pixel_conv_vs_det_diff_more = h_neutron_counters->addCounter("count_pixel_conv_vs_det_diff_more");

  GriffAnaUtils::SegmentIterator segments_gas(&dr); // track or segment filters that will make looping very efficient
  segments_gas.addFilter(new GriffAnaUtils::SegmentFilter_Volume("CountingGas"));
*/

  while (dr.loopEvents())
  {
/*
    double edep = 0;
    while (auto segment = segments_gas.next())
    {
      edep += segment->eDep();
    }
*/
    while (auto neutron = primary_neutrons.next())
    {
      count_initial_neutrons += 1;

      /*
      auto segBegin = neutron->segmentBegin();
      auto stepFirst = segBegin->firstStep();
      const double initialPosition[] = {stepFirst->preGlobalX(), stepFirst->preGlobalY(), stepFirst->preGlobalZ()};
*/

      if (hit.eventHasHit())
      {
        count_neutrons_hit += 1;

        //segments_gas.reset();
        //auto gasSegment = segments_gas.next();
        /// volumeCopyNumber() = CountingGas; volumeCopyNumber(1) = Converter; volumeCopyNumber(2) = straw wall; volumeCopyNumber(3) = EmptyTube;
        /// volumeCopyNumber(4) = TubeWall; volumeCopyNumber(5) = EmptyPackBox; volumeCopyNumber(6) = Bank; volumeCopyNumber(7) = World
        //const int strawId_gas = gasSegment->volumeCopyNumber(2);
        //const int tubeId_gas = gasSegment->volumeCopyNumber(4);


        //const int pixelId_gas = getPixelId(tubeId_gas, strawId_gas, gasSegment->firstStep()->preGlobalX());

        auto segL = neutron->lastSegment();
        /// volumeCopyNumber() = CountingGas; volumeCopyNumber(1) = Converter; volumeCopyNumber(2) = straw wall; volumeCopyNumber(3) = EmptyTube;
        /// volumeCopyNumber(4) = TubeWall; volumeCopyNumber(5) = EmptyPackBox; volumeCopyNumber(6) = Bank; volumeCopyNumber(7) = World
        const int strawId_conv = segL->volumeCopyNumber(1);
        const int tubeId_conv = segL->volumeCopyNumber(3);

        //const int pixelId_conv = getPixelId(tubeId_conv, strawId_conv, segL->firstStep()->preGlobalX());

        const int pixelId_det = getPixelId(tubeId_conv, strawId_conv, hit.eventHitPositionX());

        /*
        //if(pixelId_gas != pixelId_det) {
        if(tubeId_gas != tubeId_conv) {
          printf("%d \t %d\n", tubeId_gas, tubeId_conv);
          count_pixel_gas_vs_det_diff += 1;
          if(abs(pixelId_gas - pixelId_det)>1){
            count_pixel_gas_vs_det_diff_more += 1;
          }
        }

        if (pixelId_conv != pixelId_det) {
          count_pixel_conv_vs_det_diff += 1;
          if (abs(pixelId_conv - pixelId_det) > 1) {
            count_pixel_conv_vs_det_diff_more += 1;
          }
        }
        */
        mcplParticle->time = hit.eventHitTime() / Units::ms;
        mcplParticle->weight = hit.eventHitWeight();
        mcplParticle->userflags = pixelId_det + IDFdetectorPixelOffset;

        mcplParticle->position[0] = 0;  //position_hit[0] / Units::cm; //not used
        mcplParticle->position[1] = 0;  //position_hit[1] / Units::cm; //not used
        mcplParticle->position[2] = 0;  //position_hit[2] / Units::cm; //not used
        mcplParticle->direction[0] = 0; //dummy
        mcplParticle->direction[1] = 0; //dummy
        mcplParticle->direction[2] = 1; //dummy
        mcplParticle->ekin = 0;         //dummy
        mcplParticle->pdgcode = 0;      //dummy

        mcpl_add_particle(detMcpl, mcplParticle);

      } // end hit in event
      //} // end if conversion condition

    } //end of loop over primary neutrons
  }   //end of event loop

  mcpl_close_outfile(detMcpl);
  hc.saveToFile("larmor_det_events", true);

  return 0;
}
