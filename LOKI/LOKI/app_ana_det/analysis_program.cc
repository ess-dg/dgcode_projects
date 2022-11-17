#include "GriffAnaUtils/All.hh"
#include "Core/Units.hh"
#include "Core/FPE.hh"
#include "Utils/ArrayMath.hh"
#include "Utils/NeutronMath.hh"
#include "SimpleHists/HistCollection.hh"
#include "GriffB10Common/DetHitApproximation.hh"
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include "MCPL/mcpl.h"

#include "G4GeoLoki/PixelatedBanks.hh"
//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griff for more info.

#ifndef M_PI
#define M_PI  3.14159265358979323846  //  pi
#endif

int main(int argc, char**argv) {
  
  Core::catch_fpe();
  GriffDataReader dr(argc,argv);

  auto setup = dr.setup();
  auto &geo = setup->geo();
  //printf("QQQ=============  %s \n", geo.getName().c_str());
  if (geo.getName()!="G4GeoLoki/GeoBCSBanks" && geo.getName()!="G4GeoBCS/GeoLarmorBCSExperiment") {
    printf("Error: Wrong setup for this analysis\n");
    return 1;
  }

  /*
  auto &gen = setup->gen();
  double sourceSampleDistance = 0*Units::m;
  if (gen.getName()=="G4MCPLPlugins/MCPLGen" && geo.getName()!="G4GeoBCS/GeoLarmorBCSExperiment") { // . changed to /
    sourceSampleDistance = 23.5966*Units::m; //changed from 22.5, might change //23.5966 from McStas mcdisplay-webgl
  }*/

  setup->dump();

  DetHitApproximation hit(&dr,1.2*Units::cm,120*Units::keV,"CountingGas" ); // defined for He3 gas

  GriffAnaUtils::TrackIterator primary_neutrons(&dr);
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  mcpl_outfile_t detMcpl;
  mcpl_particle_t *mcplParticle;

  detMcpl = mcpl_create_outfile("detectionEvents.mcpl");
  mcpl_hdr_add_comment(detMcpl, "Neutrons in this file are actually detection events and userflags are pixel ID's of hits. Created with ess_loki_loki_ana command.");
  mcpl_enable_userflags(detMcpl);
  mcplParticle = mcpl_get_empty_particle(detMcpl);

  SimpleHists::HistCollection hc;

  auto userData = setup->userData();
  PixelatedBanks* banks;
  const double sampleDetectorDistance = setup->geo().getParameterDouble("rear_detector_distance_m") *Units::m;
  if(userData.count("analysis_straw_pixel_number")){
    const int strawPixelNumber = std::stoi(userData["analysis_straw_pixel_number"].c_str());
    banks = new PixelatedBanks(sampleDetectorDistance, strawPixelNumber);
  }
  else{ // use default rear bank pixel number
    banks = new PixelatedBanks(sampleDetectorDistance);
  }
 
  auto h_neutron_xy_hit = hc.book2D("Neutron xy (hit)", 2500, -1250, 1250, 2500, -1250, 1250, "neutron_xy_hit");
       h_neutron_xy_hit->setXLabel("-x [mm]");
       h_neutron_xy_hit->setYLabel("y [mm]");
  
  auto h_neutron_counters = hc.bookCounts("General neutron counters","neutron_counters"); /////////////
  //auto count_entering_aluminium = h_neutron_counter->addCounter("count_entering_aluminium");
  auto count_initial_neutrons = h_neutron_counters->addCounter("count_initial_neutrons");
  auto count_neutrons_hit = h_neutron_counters->addCounter("count_neutrons_hit");
 

  while (dr.loopEvents()) {  
    
    while (auto neutron = primary_neutrons.next()) {
      count_initial_neutrons += 1;

      if (hit.eventHasHit()) {
        count_neutrons_hit += 1;
        
        auto segL = neutron->lastSegment();
        /// volumeCopyNumber() = CountingGas; volumeCopyNumber(1) = Converter; volumeCopyNumber(2) = straw wall; volumeCopyNumber(3) = EmptyTube;
        /// volumeCopyNumber(4) = TubeWall; volumeCopyNumber(5) = EmptyPackBox; volumeCopyNumber(6) = Bank; volumeCopyNumber(7) = World
        const int strawId_conv = segL->volumeCopyNumber(1);
        const int tubeId_conv = segL->volumeCopyNumber(3);
        const int bankId_conv = segL->volumeCopyNumber(5);
        
        h_neutron_xy_hit->fill(-hit.eventHitPositionX()/Units::mm, hit.eventHitPositionY()/Units::mm, hit.eventHitWeight());

        const int pixelId = banks->getPixelId(bankId_conv, tubeId_conv, strawId_conv, hit.eventHitPositionX(), hit.eventHitPositionY());
        
        mcplParticle->time = hit.eventHitTime()/Units::ms;
        mcplParticle->weight = hit.eventHitWeight();
        mcplParticle->userflags = pixelId; 
        mcplParticle->position[0] = 0; //not used
        mcplParticle->position[1] = 0; //not used
        mcplParticle->position[2] = 0; //not used
        mcplParticle->direction[0] = 0; //dummy
        mcplParticle->direction[1] = 0; //dummy
        mcplParticle->direction[2] = 1; //dummy
        mcplParticle->ekin = 0; //dummy
        mcplParticle->pdgcode = 0; //dummy
        
        mcpl_add_particle(detMcpl, mcplParticle);
      } // end hit in event


    }//end of loop over primary neutrons
  } //end of event loop

  //mcpl_close_outfile(f);
  
  mcpl_close_outfile(detMcpl);

  hc.saveToFile("bcsloki_sans", true);
    
  return 0;
}


