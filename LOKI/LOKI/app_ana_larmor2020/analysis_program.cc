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
#include "LOKI/DetectionFileCreator.hh"

#include "G4GeoLoki/BcsTube.hh"
//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griff for more info.

#ifndef M_PI
#define M_PI  3.14159265358979323846  //  pi
#endif

const int IDFdetectorPixelOffset = 11;
const int strawPixelNumber = 512;
const int strawNumber = 4 * 8 * 7;

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
  bool createDetectionMcplFile = true;
  std::string tmp;
  for (int i=1;i<argc;++i) {
    tmp=argv[i];
    if (tmp=="--nodetfile" ) {
      createDetectionMcplFile = false;
    }
  }
  
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
  else if(gen.getName()=="LOKI.FloodSourceGen/FloodSourceGen"){
    sourceSampleDistance = gen.getParameterDouble("source_sample_distance_meters") *Units::m;
  }
  printf("sourceSampleDistance: %f\n", sourceSampleDistance);

  setup->dump();

  DetHitApproximation hit(&dr,1.2*Units::cm,120*Units::keV,"CountingGas" ); // defined for He3 gas

  GriffAnaUtils::SegmentIterator segments_gas(&dr);
  segments_gas.addFilter(new GriffAnaUtils::SegmentFilter_Volume("CountingGas"));  

  GriffAnaUtils::TrackIterator primary_neutrons(&dr);
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  GriffAnaUtils::SegmentIterator segments_World(&dr);   
  segments_World.addFilter(new GriffAnaUtils::SegmentFilter_Volume("World"));
  segments_World.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  segments_World.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  GriffAnaUtils::SegmentIterator segments_TubeWall(&dr); 
  segments_TubeWall.addFilter(new GriffAnaUtils::SegmentFilter_Volume("TubeWall"));
  segments_TubeWall.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  segments_TubeWall.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  DetectionFileCreator* detectionFile = nullptr;
  if (createDetectionMcplFile == true) {
    detectionFile = new DetectionFileCreator("detectionEvents.mcpl");
  }

  SimpleHists::HistCollection hc;

  const double sampleDetectorDistance = setup->geo().getParameterDouble("rear_detector_distance_m") *Units::m;

  const double tubeRadius = BcsTube::getTubeOuterRadius(); //12.7;

  //float xmin = -53;
  const double ymin = -53; //20+1 tube in negative direction
  //int tubeNr = 40 + 2; // 40 plus 1-1 empty space on both ends to make the plot nice
  //int binsx = 1060;
  const int binsy = 1060/2;
  const int binsz = 160;
  int thetabins = 150/2;
  const double trueThetaMax = 15;
  //float thetamax = 9.0;
  //float dthetamin = -5; //-0.35
  //float dthetamax = 5; //3 
  //int dthetabins = 1000;
  const double zmin = (sampleDetectorDistance - tubeRadius) - 20; //[mm]  //4980 for 5 m sd distance
  const double zmax = (sampleDetectorDistance - tubeRadius) + 140; //[mm //]5140 for 5 m sd distance
  //int zbins = 160;

  auto h_neutron_xy_conv = hc.book2D("Neutron xy (conv)", 2500, -1250, 1250, 2500, -1250, 1250, "neutron_xy_conv");
       h_neutron_xy_conv->setXLabel("-x [mm]");
       h_neutron_xy_conv->setYLabel("y [mm]");
  
  auto h_neutron_xy_hit = hc.book2D("Neutron xy (hit)", 2500, -1250, 1250, 2500, -1250, 1250, "neutron_xy_hit");
       h_neutron_xy_hit->setXLabel("-x [mm]");
       h_neutron_xy_hit->setYLabel("y [mm]");

  auto h_neutron_zy_conv = hc.book2D("Neutron zy (conv)", binsz, zmin, zmax, binsy/2, ymin, -ymin, "neutron_yz_conv");
       h_neutron_zy_conv->setXLabel("z [mm]");
       h_neutron_zy_conv->setYLabel("y [cm]");

  auto h_neutron_zy_hit = hc.book2D("Neutron zy (hit)", binsz, zmin, zmax, binsy/2, ymin, -ymin, "neutron_yz_hit");
       h_neutron_zy_hit->setXLabel("z [mm]");
       h_neutron_zy_hit->setYLabel("y [cm]");

  //
  auto h_neutron_zy_big_conv = hc.book2D("Neutron zy big scale (conv)", 5200, 0, 5200, 2200, -1000, 1200, "neutron_zy_big_conv");
       h_neutron_zy_big_conv->setXLabel("z [mm]");
       h_neutron_zy_big_conv->setYLabel("y [mm]");

  auto h_neutron_zx_big_conv = hc.book2D("Neutron zx big scale (conv)", 5200, 0, 5200, 2500, -1250, 1250, "neutron_zx_big_conv");
       h_neutron_zx_big_conv->setXLabel("z [mm]");
       h_neutron_zx_big_conv->setYLabel("x [mm]");

  auto h_neutron_theta = hc.book1D("Neutron true theta (all events)", thetabins, 0, trueThetaMax, "neutron_true_theta");
       h_neutron_theta->setXLabel("Angle [degree]");

  auto h_neutron_theta_conv = hc.book1D("Neutron theta (conv)", thetabins, 0, trueThetaMax, "neutron_theta_conv");
       h_neutron_theta_conv->setXLabel("Angle [degree]");

  auto h_neutron_theta_hit = hc.book1D("Neutron theta (hit)", thetabins, 0, trueThetaMax, "neutron_theta_hit");
       h_neutron_theta_hit->setXLabel("Angle [degree]");
  
  auto h_neutron_bank_theta_hit = hc.book2D("Neutron Theta for banks [degree](hit)", thetabins, 0, trueThetaMax, 9, 0, 9, "neutron_bank_theta_hit");
       h_neutron_bank_theta_hit->setXLabel("Theta [degree]");
       h_neutron_bank_theta_hit->setYLabel("Bank id");

  //const int numberOfPixels = strawNumber * strawPixelNumber; //banks.getTotalNumberOfPixels();
  //printf("Number of pixels: %d\n", numberOfPixels);
  // auto h_neutron_pixel_hit_count = hc.book1D("Number of hits in pixels (all banks)", numberOfPixels, 0, numberOfPixels, "neutron_pixel_hit_count");
  //      h_neutron_pixel_hit_count->setXLabel("Pixel ID");
  // auto h_neutron_pixel_hit_weight = hc.book1D("Sum weight of hits in pixels (all banks)", numberOfPixels, 0, numberOfPixels, "neutron_pixel_hit_weight");
  //      h_neutron_pixel_hit_weight->setXLabel("Pixel ID");

  

  auto h_neutron_pixel_hit = hc.book2D("Sum weight of hits in pixels (hit)", strawPixelNumber, 0, strawPixelNumber, strawNumber, 0, strawNumber, "h_neutron_pixel_hit");
       h_neutron_pixel_hit->setXLabel("Pixel ID along straw");
       h_neutron_pixel_hit->setYLabel("Straw ID");

/*
  const int numberOfPixels_rear = banks.getNumberOfPixels(0);
  auto h_neutron_pixel_rear_hit = hc.book2D("Sum weight of hits in pixels of rear bank (hit)", 256, 0, 256, numberOfPixels_rear/256, 0, numberOfPixels_rear/256, "h_neutron_pixel_rear_hit");
       h_neutron_pixel_rear_hit->setXLabel("Pixel ID along straw");
       h_neutron_pixel_rear_hit->setYLabel("Straw ID");
*/

  auto h_neutron_segment_number_World = hc.book1D("Number of World segments per neutron", 50, 0, 50, "neutron_segment_number_World");

  auto h_neutron_bankIncidentCounter = hc.book1D("Incident neutron counter for the 9 banks", 9, 0, 9, "neutron_bankIncidentCounter");
       h_neutron_bankIncidentCounter->setXLabel("Bank id");

  auto h_neutron_nonFirstBankIncident = hc.book1D("Incident neutrons for banks (only neutrons which had entered different banks previousl)", 9, 0, 9, "neutron_nonFirstBankIncident");
       h_neutron_nonFirstBankIncident->setXLabel("Bank id");

  auto h_neutron_bankToBankIncident = hc.book2D("Incident neutrons for banks (only neutrons which had entered different banks previousl)", 9, 0, 9, 9, 0, 9,"neutron_bankToBankIncident");
       h_neutron_bankToBankIncident->setXLabel("Bank id (hit)");
       h_neutron_bankToBankIncident->setYLabel("Bank id (previous)");

  auto h_neutron_bankPanelConvCounter = hc.book2D("Neutron conversion counter for 4 panels of 9 banks", 9, 0, 9, 4, 0, 4, "neutron_bankPanelConvCounter");
       h_neutron_bankPanelConvCounter->setXLabel("Bank id");
       h_neutron_bankPanelConvCounter->setYLabel("Panel number");

  auto h_neutron_bankPanelHitCounter = hc.book2D("Neutron hit counter for 4 panels of 9 banks", 9, 0, 9, 4, 0, 4, "neutron_bankPanelHitCounter");
       h_neutron_bankPanelHitCounter->setXLabel("Bank id");
       h_neutron_bankPanelHitCounter->setYLabel("Panel number");

  auto h_neutron_panelConvCounter = hc.book1D("Neutron conversion counter for 4 panels", 4, 0, 4, "neutron_panelConvCounter");
       h_neutron_panelConvCounter->setXLabel("Panel number");

  auto h_neutron_panelHitCounter = hc.book1D("Neutron hit counter for 4 panels", 4, 0, 4, "neutron_panelHitCounter");
       h_neutron_panelHitCounter->setXLabel("Panel number");

  auto h_panel_neutronNr_all = hc.book1D("Number of neutrons entering each panel", 4, 0, 4, "panel_neutronNr_all");

  auto h_panel_lambda = hc.book2D("Incident neutron wavelength for panels", 140, 0, 14, 4, 0, 4, "panel_lambda");
       h_panel_lambda->setXLabel("Wavelength [angstrom]");
       h_panel_lambda->setYLabel("Panel number");

  auto h_bank_lambda = hc.book2D("Incident neutron wavelength for banks", 140, 0, 14, 9, 0, 9, "bank_lambda");
       h_bank_lambda->setXLabel("Wavelength [angstrom]");
       h_bank_lambda->setYLabel("Bank id");
       
  auto h_bank_lambda_hit = hc.book2D("Neutron wavelength for banks (hit)", 140, 0, 14, 9, 0, 9, "bank_lambda_hit");
       h_bank_lambda_hit->setXLabel("Wavelength [angstrom]");
       h_bank_lambda_hit->setYLabel("Bank id");

  auto h_mcpl_lambda_calculated = hc.book1D("Neutron wavelength calculated from TOF in MCPL", 140, 0, 14, "mcpl_lambda_calculated");

  auto h_neutron_Q_hit = hc.book1D("Neutron Q [1/angstrom] (hit)", 250, 0, 2, "neutron_Q_hit");
       h_neutron_Q_hit->setXLabel("Q [1/angstrom]");

  auto h_neutron_bank_Q_hit = hc.book2D("Neutron Q for banks [1/angstrom] (hit)", 250, 0, 2, 9, 0, 9, "neutron_bank_Q_hit");
       h_neutron_bank_Q_hit->setXLabel("Q [1/angstrom]");
       h_neutron_bank_Q_hit->setYLabel("Bank id");
       
       
  
  auto h_neutron_counters = hc.bookCounts("General neutron counters","neutron_counters"); /////////////
  //auto count_entering_aluminium = h_neutron_counter->addCounter("count_entering_aluminium");
  auto count_initial_neutrons = h_neutron_counters->addCounter("count_initial_neutrons");
  auto count_neutrons_entering_aluminium = h_neutron_counters->addCounter("count_neutrons_entering_aluminium");
  auto count_neutrons_converted = h_neutron_counters->addCounter("count_neutrons_converted");
  auto count_neutrons_hit = h_neutron_counters->addCounter("count_neutrons_hit");


  while (dr.loopEvents()) {  
    
    
    while (auto neutron = primary_neutrons.next()) {
      count_initial_neutrons += 1;

      auto segBegin = neutron->segmentBegin();
      auto stepFirst = segBegin->firstStep();
      const double initialPosition[] = {stepFirst->preGlobalX(), stepFirst->preGlobalY(), stepFirst->preGlobalZ()};

      auto dir_true = stepFirst->preMomentumArray();
      const double theta_true = Utils::theta(dir_true)/Units::degree;
      h_neutron_theta->fill(theta_true, neutron->weight());

      if (gen.getName()=="G4MCPLPlugins/MCPLGen"){ //works only for non-zero initial TOF
        double velocity_mcpl_calculated = (sourceSampleDistance / Units::m) / (stepFirst->preTime() / Units::s);
        const double lambda_mcpl_calculated = Utils::neutron_meters_per_second_to_angstrom(velocity_mcpl_calculated);

        h_mcpl_lambda_calculated->fill(lambda_mcpl_calculated, neutron->weight());
      }

      int seg_count_World = 0;
      double seg_length_World = 0.0;
      double neutron_weight = 0.0;
      //int isBackScattered = 0;

      while (auto segment = segments_World.next()) { // loop over all segments of the World
        ++seg_count_World;
        seg_length_World += segment->segmentLength();
        neutron_weight = segment->getTrack()->weight();
        
        /*if (seg_count_World > 1) {
          auto step0 = segment->firstStep();
          double dir_init[3];
          Utils::normalise(step0->preMomentumArray(), dir_init);
          auto pos_init = step0->preGlobalArray();
          
          if (dir_init[2] > 0) {
            //isBackScattered = 1;

            p->position[0] = pos_init[0] / Units::cm; //mcpl file unit
            p->position[1] = pos_init[1] / Units::cm; //mcpl file unit
            p->position[2] = pos_init[2] / Units::cm; //mcpl file unit
            p->direction[0] = dir_init[0];
            p->direction[1] = dir_init[1];
            p->direction[2] = dir_init[2];
            p->ekin = step0->preEKin();
            p->time = step0->preTime();
            p->weight = neutron->weight();
            p->pdgcode = segment->getTrack()->pdgCode();

            mcpl_add_particle(f, p);
          }
        }*/
      }
      if (seg_length_World /*&& isBackScattered == 1*/) {
        h_neutron_segment_number_World->fill(seg_count_World, neutron_weight);
      }

      bool firstBankEnter = true;

      int previousPanelNumber = -42; //anything except [0-6] would do it
      int panelNumber = 0;
      int previousBankNumber = -42; //anything except [0-6] would do it
      int bankNumber = 0;
      while (auto tubeWallSegment = segments_TubeWall.next()) { // loop over all segments of the TubeWall

        bankNumber = (int)tubeWallSegment->volumeCopyNumber(2);
        panelNumber = (int)tubeWallSegment->volumeCopyNumber() % 4;
        neutron_weight = tubeWallSegment->getTrack()->weight();


        if (panelNumber != previousPanelNumber && tubeWallSegment->startEKin()) {
          previousPanelNumber = panelNumber;

          h_panel_neutronNr_all->fill(panelNumber, neutron_weight);

          const double actualEkin = tubeWallSegment->startEKin();
          const double actualLambda = Utils::neutronEKinToWavelength(actualEkin) / Units::angstrom;
          h_panel_lambda->fill(actualLambda, panelNumber, neutron_weight);
        }
        if (bankNumber!= previousBankNumber) { //his includes first enter
          if (!firstBankEnter) {
            h_neutron_nonFirstBankIncident->fill(bankNumber, neutron_weight);
            h_neutron_bankToBankIncident->fill(bankNumber, previousBankNumber, neutron_weight);
          }
          firstBankEnter = false;

          previousBankNumber = bankNumber;
          h_neutron_bankIncidentCounter->fill(bankNumber, neutron_weight);
          
          const double actualEkin = tubeWallSegment->startEKin();
          const double actualLambda = Utils::neutronEKinToWavelength(actualEkin) / Units::angstrom;
          h_bank_lambda->fill(actualLambda, bankNumber, neutron_weight);
        }
      }
      segments_TubeWall.reset(); //it is needed

      if (segments_TubeWall.next()) {
        count_neutrons_entering_aluminium += 1;
      }
      segments_TubeWall.reset(); //it is done for each event automatically


      auto segL = neutron->lastSegment();
      if (segL->volumeName()=="Converter") {
        count_neutrons_converted += 1;
        
        auto stepL = segL->lastStep();
        const double position_conv[3] = {stepL->postGlobalX(), stepL->postGlobalY(), stepL->postGlobalZ()};

        double dir_conv[3];
        Utils::subtract(stepL->postGlobalArray(), stepFirst->preGlobalArray(), dir_conv);
        const double theta_conv = Utils::theta(dir_conv)/Units::degree;
        h_neutron_theta_conv->fill(theta_conv, neutron->weight());
            
        /// volumeCopyNumber() = CountingGas; volumeCopyNumber(1) = Converter; volumeCopyNumber(2) = straw wall; volumeCopyNumber(3) = EmptyTube;
        /// volumeCopyNumber(4) = TubeWall; volumeCopyNumber(5) = EmptyPackBox; volumeCopyNumber(6) = Bank; volumeCopyNumber(7) = World
        const int strawId_conv = segL->volumeCopyNumber(1);
        const int tubeId_conv = segL->volumeCopyNumber(3);
        const int bankId_conv = segL->volumeCopyNumber(5);

        

        h_neutron_zy_conv->fill(position_conv[2]/Units::mm, position_conv[1]/Units::cm, neutron->weight());
        h_neutron_zy_big_conv->fill(position_conv[2]/Units::mm, position_conv[1]/Units::mm, neutron->weight());

        h_neutron_zx_big_conv->fill(position_conv[2]/Units::mm, position_conv[0]/Units::mm, neutron->weight());

        h_neutron_xy_conv->fill(-position_conv[0]/Units::mm, position_conv[1]/Units::mm, neutron->weight());

        const int panelNumber_conv = tubeId_conv % 4;
        h_neutron_bankPanelConvCounter->fill(bankId_conv, panelNumber_conv, neutron->weight());
        h_neutron_panelConvCounter->fill(panelNumber_conv, neutron->weight());

        if (hit.eventHasHit()) {
          count_neutrons_hit += 1;
          const double position_hit[3] = {hit.eventHitPositionX(), hit.eventHitPositionY(), hit.eventHitPositionZ()};
 
          h_neutron_zy_hit->fill(position_hit[2]/Units::mm, position_hit[1]/Units::cm, hit.eventHitWeight());
          h_neutron_xy_hit->fill(-position_hit[0]/Units::mm, position_hit[1]/Units::mm, hit.eventHitWeight());

          const double theta_hit = Utils::theta(hit.eventHitPosition())/Units::degree;
          h_neutron_theta_hit->fill(theta_hit, hit.eventHitWeight());

          const int pixelId = getPixelId(tubeId_conv, strawId_conv, position_hit[0]);
          //h_neutron_pixel_hit_count->fill(pixelId, 1);
          //h_neutron_pixel_hit_weight->fill(pixelId, hit.eventHitWeight());
          h_neutron_pixel_hit->fill(pixelId%strawPixelNumber, std::floor(pixelId/strawPixelNumber), hit.eventHitWeight());

          
          h_neutron_bankPanelHitCounter->fill(bankId_conv, panelNumber_conv, hit.eventHitWeight());
          h_neutron_panelHitCounter->fill(panelNumber_conv, hit.eventHitWeight());

          // if (bankId_conv == 0) {
          //   h_neutron_pixel_rear_hit->fill(pixelId%256, std::floor(pixelId/256), hit.eventHitWeight());
          // }

          //TODO should implement method (in BcsBanks class) to get positionOnWire_hit coordinate. Ask Judit, how it is done in real data reduction.
          const double sampleToExactHitPositionDistance = std::sqrt(std::pow((position_hit[0] - initialPosition[0]), 2) +
                                                                    std::pow((position_hit[1] - initialPosition[1]), 2) +
                                                                    std::pow((position_hit[2] - initialPosition[2]), 2));

          const double tof_hit = hit.eventHitTime()/Units::ms;
          double velocity_calculated = -1;
          if (tof_hit > 0.0) {
            velocity_calculated = ((sampleToExactHitPositionDistance + sourceSampleDistance) / Units::m) / (hit.eventHitTime() / Units::s);
          }
          else {
            printf("Error in hit tof value, tof zero or negative \n");
            return 1;
          }
          const double lambda_hit_calculated = Utils::neutron_meters_per_second_to_angstrom(velocity_calculated);
          const double Q_hit_calculated = 4 * M_PI * sin(0.5 * theta_hit*Units::degree) / lambda_hit_calculated;
          h_neutron_Q_hit->fill(Q_hit_calculated, hit.eventHitWeight());
          h_neutron_bank_Q_hit->fill(Q_hit_calculated, bankId_conv, hit.eventHitWeight());
          h_neutron_bank_theta_hit->fill(theta_hit, bankId_conv, hit.eventHitWeight());

          h_bank_lambda_hit -> fill(lambda_hit_calculated, bankId_conv, hit.eventHitWeight());

          if (createDetectionMcplFile == true) {
            detectionFile->addDetectionEvent(pixelId + IDFdetectorPixelOffset, hit.eventHitTime()/Units::ms);
          }
        } // end hit in event
      } // end if conversion condition
    }//end of loop over primary neutrons
  } //end of event loop

  delete detectionFile;

  hc.saveToFile("larmor_bcs", true);

  return 0;
}
