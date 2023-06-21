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

#include "G4GeoLoki/PixelatedBanks.hh"
//Griff analysis. See https://confluence.esss.lu.se/display/DG/Griff for more info.

#ifndef M_PI
#define M_PI  3.14159265358979323846  //  pi
#endif

double calculateWavelength(const double* initialPosition, const double* position_hit, const double tof, const double preGeant4Distance) {
  const double sampleToPositionDistance = std::sqrt(std::pow((position_hit[0] - initialPosition[0]), 2) +
                                                    std::pow((position_hit[1] - initialPosition[1]), 2) +
                                                    std::pow((position_hit[2] - initialPosition[2]), 2));
  double velocity_calculated = -1;
  if (tof > 0.0) {
    velocity_calculated = ((sampleToPositionDistance + preGeant4Distance)/Units::m) / (tof/Units::s);
  }
  else {
    printf("Error in hit tof value, tof zero or negative \n");
    return 0;
  }
  return Utils::neutron_meters_per_second_to_angstrom(velocity_calculated);
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
  auto userData = setup->userData();

  if (geo.getName()!="G4GeoLoki/GeoBCSBanks" && geo.getName()!="G4GeoBCS/GeoLarmorBCSExperiment") {
    printf("Error: Wrong setup for this analysis\n");
    return 1;
  }
  const double rearDetectorDistance = setup->geo().getParameterDouble("rear_detector_distance_m") *Units::m;
  if (rearDetectorDistance != 4.099*Units::m) {
    printf("Error: rear_detector_distance_m is not correct for the Larmor 2022 experiment!\n");
    return 1;
  }

  auto &gen = setup->gen();

  double preGeant4Distance = 0*Units::m;
  double nominalSamplePosDistance = 25.61*Units::m;
  if(gen.getName()=="LOKI.FloodSourceGen/FloodSourceGen" || gen.getName()=="G4MCPLPlugins/MCPLGen") {
    if(userData.count("nominal_source_sample_distance_meters")) {
      nominalSamplePosDistance = std::stod(userData["nominal_source_sample_distance_meters"].c_str())*Units::m;
    }
    else if (gen.hasParameterDouble("source_sample_distance_meters")) { //old variable name, kept for backward compatibility
      nominalSamplePosDistance = gen.getParameterDouble("source_sample_distance_meters") *Units::m;
    }
    double nominalSamplePosToGeneratorDistance = 0.0; //distance between the nominal sample position and the particle generator
    if(gen.getName()=="G4MCPLPlugins/MCPLGen") {
      nominalSamplePosToGeneratorDistance = 0.2; //default kept for old MCPL files
      if(userData.count("sample_mcpl_distance_m")) {
        nominalSamplePosToGeneratorDistance = std::stod(userData["sample_mcpl_distance_m"].c_str()) *Units::m;
      }
    }
    else if(gen.hasParameterDouble("gen_z_offset_meters")) {
      nominalSamplePosToGeneratorDistance = gen.getParameterDouble("gen_z_offset_meters") *Units::m;
    }
    preGeant4Distance = nominalSamplePosDistance + nominalSamplePosToGeneratorDistance; //approximation, mainly ignoring x and y
  }
  printf("preGeant4Distance: %f\n", preGeant4Distance);

  setup->dump();

  DetHitApproximation hit(&dr,1.2*Units::cm,120*Units::keV,"CountingGas" ); // defined for He3 gas

  GriffAnaUtils::TrackIterator primary_neutrons(&dr);
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  GriffAnaUtils::SegmentIterator segments_bank(&dr);
  segments_bank.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Bank"));
  segments_bank.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  segments_bank.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  GriffAnaUtils::SegmentIterator segments_TubeWall(&dr);
  segments_TubeWall.addFilter(new GriffAnaUtils::SegmentFilter_Volume("TubeWall"));
  segments_TubeWall.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  segments_TubeWall.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  GriffAnaUtils::SegmentIterator segments_StrawWall(&dr);
  segments_StrawWall.addFilter(new GriffAnaUtils::SegmentFilter_Volume("StrawWall"));
  segments_StrawWall.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  segments_StrawWall.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));


  SimpleHists::HistCollection hc;

  PixelatedBanks* banks;
  if(userData.count("analysis_straw_pixel_number")){
    const int strawPixelNumber = std::stoi(userData["analysis_straw_pixel_number"].c_str());
    banks = new PixelatedBanks(rearDetectorDistance, strawPixelNumber);
  }
  else{ // use default rear bank pixel number
    banks = new PixelatedBanks(rearDetectorDistance);
  }

  DetectionFileCreator* detectionFile = nullptr;
  if (createDetectionMcplFile == true) {
    detectionFile = new DetectionFileCreator("detectionEvents.mcpl", userData);
  }

  bool oldTubeNumbering = false;
  if (!geo.hasParameterBoolean("old_tube_numbering") || geo.getParameterBoolean("old_tube_numbering")) {
    oldTubeNumbering = true;
  }

  const double tubeRadius = BcsTube::getTubeOuterRadius(); //12.7;

  const double ymin = -53; //20+1 tube in negative direction
  const int binsy = 1060/2;
  const int binsz = 160;
  int thetabins = 550/2;
  const double trueThetaMax = 55;

  const double zmin = (rearDetectorDistance - tubeRadius) - 20; //[mm]  //4980 for 5 m sd distance
  const double zmax = (rearDetectorDistance - tubeRadius) + 140; //[mm //]5140 for 5 m sd distance

  auto h_neutron_xy_bank = hc.book2D("Neutron xy (at bank entering)", 2500, -1250, 1250, 2500, -1250, 1250, "neutron_xy_bank");
       h_neutron_xy_bank->setXLabel("-x [mm]");
       h_neutron_xy_bank->setYLabel("y [mm]");

  auto h_neutron_xy_bank_filtered = hc.book2D("Neutron xy filtered (at bank entering)", 2500, -1250, 1250, 2500, -1250, 1250, "neutron_xy_bank_filtered");
       h_neutron_xy_bank_filtered->setXLabel("-x [mm]");
       h_neutron_xy_bank_filtered->setYLabel("y [mm]");

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


  auto h_neutron_theta = hc.book1D("Neutron true theta (all events)", thetabins, 0, trueThetaMax, "neutron_true_theta");
       h_neutron_theta->setXLabel("Angle [degree]");

  auto h_neutron_theta_hit = hc.book1D("Neutron theta (hit)", thetabins, 0, trueThetaMax, "neutron_theta_hit");
       h_neutron_theta_hit->setXLabel("Angle [degree]");

  const int numberOfPixels = banks->getNumberOfPixels(0);
  const int numberOfStraws = banks->getNumberOfTubes(0) * 7;
  const int numberOfPixelsPerStraw = numberOfPixels / numberOfStraws;
  printf("Number of pixels per straw: %d \t Total number of pixels: %d \n", numberOfPixelsPerStraw, numberOfPixels);

  auto h_neutron_pixel_hit = hc.book2D("Sum weight of hits in pixels (hit)", numberOfPixelsPerStraw, 0, numberOfPixelsPerStraw, numberOfStraws, 0, numberOfStraws, "h_neutron_pixel_hit");
       h_neutron_pixel_hit->setXLabel("Pixel ID along straw");
       h_neutron_pixel_hit->setYLabel("Straw ID");

  auto h_neutron_layerConvCounter = hc.book1D("Neutron conversion counter for 4 layers", 4, 0, 4, "neutron_layerConvCounter");
       h_neutron_layerConvCounter->setXLabel("Layer number");

  auto h_neutron_layerHitCounter = hc.book1D("Neutron hit counter for 4 layers", 4, 0, 4, "neutron_layerHitCounter");
       h_neutron_layerHitCounter->setXLabel("Layer number");

  auto h_mcpl_lambda_true = hc.book1D("Neutron true wavelength from MCPL", 325, 0, 14, "mcpl_lambda_true");
  auto h_mcpl_lambda = hc.book1D("Neutron wavelength calculated from TOF in MCPL", 325, 0, 14, "mcpl_lambda");

  auto h_incident_lambda_true = hc.book1D("Incident neutron true wavelength", 325, 0, 14, "incident_lambda_true");
  auto h_incident_lambda = hc.book1D("Incident neutron wavelength calculated from TOF", 325, 0, 14, "incident_lambda");


  auto h_layer0_lambda_hit = hc.book1D("Detection neutron wavelength for layer 0", 325, 0, 14, "layer0_lambda_hit");
       h_layer0_lambda_hit->setXLabel("Wavelength [angstrom]");
  auto h_layer1_lambda_hit = hc.book1D("Detection neutron wavelength for layer 1", 325, 0, 14, "layer1_lambda_hit");
       h_layer1_lambda_hit->setXLabel("Wavelength [angstrom]");
  auto h_layer2_lambda_hit = hc.book1D("Detection neutron wavelength for layer 2", 325, 0, 14, "layer2_lambda_hit");
       h_layer2_lambda_hit->setXLabel("Wavelength [angstrom]");
  auto h_layer3_lambda_hit = hc.book1D("Detection neutron wavelength for layer 3", 325, 0, 14, "layer3_lambda_hit");
       h_layer3_lambda_hit->setXLabel("Wavelength [angstrom]");

  auto h_layer0_lambda_true_hit = hc.book1D("Detection neutron true wavelength for layer 0", 325, 0, 14, "layer0_lambda_true_hit");
       h_layer0_lambda_true_hit->setXLabel("Wavelength [angstrom]");
  auto h_layer1_lambda_true_hit = hc.book1D("Detection neutron true wavelength for layer 1", 325, 0, 14, "layer1_lambda_true_hit");
       h_layer1_lambda_true_hit->setXLabel("Wavelength [angstrom]");
  auto h_layer2_lambda_true_hit = hc.book1D("Detection neutron true wavelength for layer 2", 325, 0, 14, "layer2_lambda_true_hit");
       h_layer2_lambda_true_hit->setXLabel("Wavelength [angstrom]");
  auto h_layer3_lambda_true_hit = hc.book1D("Detection neutron true wavelength for layer 3", 325, 0, 14, "layer3_lambda_true_hit");
       h_layer3_lambda_true_hit->setXLabel("Wavelength [angstrom]");

  auto h_layer_incident_lambda = hc.book2D("Incident wavelength for layers", 4, 0, 4, 325, 0, 14, "layer_incident_lambda");
       h_layer_incident_lambda->setXLabel("Layer id");
       h_layer_incident_lambda->setYLabel("Wavelength [angstrom]");

  auto h_layer_incident_lambda_true = hc.book2D("Incident true wavelength for layers", 4, 0, 4, 325, 0, 14, "layer_incident_lambda_true");
       h_layer_incident_lambda_true->setXLabel("Layer id");
       h_layer_incident_lambda_true->setYLabel("Wavelength [angstrom]");

  auto h_layer_first_incident_lambda = hc.book2D("Incident wavelength for layers (first enter only)", 4, 0, 4, 325, 0, 14, "layer_first_incident_lambda");
       h_layer_first_incident_lambda->setXLabel("Layer id");
       h_layer_first_incident_lambda->setYLabel("Wavelength [angstrom]");

  auto h_layer_first_incident_lambda_true = hc.book2D("Incident true wavelength for layers (first enter only)", 4, 0, 4, 325, 0, 14, "layer_first_incident_lambda_true");
       h_layer_first_incident_lambda_true->setXLabel("Layer id");
       h_layer_first_incident_lambda_true->setYLabel("Wavelength [angstrom]");

  auto h_straw_incident_lambda_true = hc.book2D("Incident true wavelength for straws", numberOfStraws, 0, numberOfStraws, 325, 0, 14, "straw_incident_lambda_true");
       h_straw_incident_lambda_true->setXLabel("Straw id");
       h_straw_incident_lambda_true->setYLabel("Wavelength [angstrom]");

  auto h_straw_incident_lambda = hc.book2D("Incident wavelength for straws", numberOfStraws, 0, numberOfStraws, 325, 0, 14, "straw_incident_lambda");
       h_straw_incident_lambda->setXLabel("Straw id");
       h_straw_incident_lambda->setYLabel("Wavelength [angstrom]");

  auto h_straw_lambda_true_hit = hc.book2D("Detection true wavelength for straws ", numberOfStraws, 0, numberOfStraws, 325, 0, 14, "straw_lambda_true_hit");
       h_straw_lambda_true_hit->setXLabel("Straw id");
       h_straw_lambda_true_hit->setYLabel("Wavelength [angstrom]");

  auto h_straw_lambda_hit = hc.book2D("Detection wavelength for straws ", numberOfStraws, 0, numberOfStraws, 325, 0, 14, "straw_lambda_hit");
       h_straw_lambda_hit->setXLabel("Straw id");
       h_straw_lambda_hit->setYLabel("Wavelength [angstrom]");

  auto h_neutron_counters = hc.bookCounts("General neutron counters","neutron_counters"); /////////////
  auto count_initial_neutrons = h_neutron_counters->addCounter("count_initial_neutrons");
  auto count_neutrons_converted = h_neutron_counters->addCounter("count_neutrons_converted");
  auto count_neutrons_hit = h_neutron_counters->addCounter("count_neutrons_hit");


  const double xWidthVacuumTankEnd = 790 *Units::mm;
  const double yHeightVacuumTankEnd = 700 *Units::mm;
  const double zVacuumTankEnd = rearDetectorDistance - 90*Units::mm; // "The front of the Loki detector was 90 mm in front of the tank”
  const double zBankFront = rearDetectorDistance - banks->detectorSystemFrontDistanceFromBankFront(0) *Units::mm;
  const double xBankEnterLimit = (xWidthVacuumTankEnd * 0.5) * (zBankFront / zVacuumTankEnd);
  const double yBankEnterLimit = (yHeightVacuumTankEnd * 0.5) * (zBankFront / zVacuumTankEnd);
  //std::cout<<"\n ***** \n xBankEnterLimit: " << xBankEnterLimit << "\n yBankEnterLimit: "<< yBankEnterLimit<< "\n ****** \n";

  double tmp_x_offset_meters = gen.hasParameterDouble("x_offset_meters") ? gen.getParameterDouble("x_offset_meters") *Units::m : 0.0;
  double tmp_dx_meter = gen.hasParameterDouble("dx_meter") ? gen.getParameterDouble("dx_meter") *Units::m : 0.0;
  if(tmp_x_offset_meters && tmp_dx_meter){
    printf("Error: both 'x_offset_meters' and 'dx_meter' parameters are defined.\n");
    return 1;
  }
  const double xBeamOffset = tmp_x_offset_meters ? tmp_x_offset_meters :tmp_dx_meter;

  while (dr.loopEvents()) {
    while (auto neutron = primary_neutrons.next()) {
      count_initial_neutrons += 1;

      auto segBegin = neutron->segmentBegin();
      auto stepFirst = segBegin->firstStep();
      const double initialPosition[] = {stepFirst->preGlobalX(), stepFirst->preGlobalY(), stepFirst->preGlobalZ()};

      auto dir_true = stepFirst->preMomentumArray();
      const double theta_true = Utils::theta(dir_true)/Units::degree;
      h_neutron_theta->fill(theta_true, neutron->weight());

      const double ekinMcpl = segBegin->startEKin();
      const double lambdaInitial = Utils::neutronEKinToWavelength(ekinMcpl) / Units::angstrom;
      h_mcpl_lambda_true->fill(lambdaInitial, neutron->weight());

      double lambdaMcplCalculated = 0.0;
      if (gen.getName()=="G4MCPLPlugins/MCPLGen"){ //works only for non-zero initial TOF
        const double velocityMcplCalculated = (preGeant4Distance / Units::m) / (stepFirst->preTime() / Units::s);
        lambdaMcplCalculated = Utils::neutron_meters_per_second_to_angstrom(velocityMcplCalculated);
        h_mcpl_lambda->fill(lambdaMcplCalculated, neutron->weight());
      }

      if(auto bankSegment = segments_bank.next()){
        auto stepFirstInBank = bankSegment->firstStep();
        const double xBankEnter = stepFirstInBank->preGlobalX()/Units::mm;
        const double yBankEnter = stepFirstInBank->preGlobalY()/Units::mm;
        h_neutron_xy_bank->fill(-xBankEnter, yBankEnter, neutron->weight());

        //if (xBankEnterLimit < std::fabs(xBankEnter) || yBankEnterLimit < std::fabs(yBankEnter)){
        if ( yBankEnterLimit < std::fabs(yBankEnter) ||
          xBankEnter > xBeamOffset + xBankEnterLimit || xBeamOffset - xBankEnterLimit > xBankEnter ){
          continue;
        }

        h_neutron_xy_bank_filtered->fill(-xBankEnter, yBankEnter, neutron->weight());

        h_incident_lambda_true->fill(lambdaInitial, neutron->weight());
        if (gen.getName()=="G4MCPLPlugins/MCPLGen"){ //works only for non-zero initial TOF
          h_incident_lambda->fill(lambdaMcplCalculated, neutron->weight());
        }
      }
      //else{
      //  std::cout<<"\n No bank segment??";
      //}

      int previousLayerId = -1;
      int firstEnterToLayer[4] = {1, 1, 1, 1};
      while (auto tubeWallSegment = segments_TubeWall.next()) {
        const int tubeId = tubeWallSegment->volumeCopyNumber();
        const int layerId = banks->getTubeLayerId(0, tubeId, oldTubeNumbering);
        const double actualEkin = tubeWallSegment->startEKin();

        auto stepF = tubeWallSegment->firstStep();
        const double position[3] = {stepF->preGlobalX(), stepF->preGlobalY(), stepF->preGlobalZ()};
        const double lambda_calculated = calculateWavelength(initialPosition, position, stepF->preTime(), preGeant4Distance);
        const double actualLambda = Utils::neutronEKinToWavelength(actualEkin)/Units::angstrom;

        if(layerId != previousLayerId && actualEkin) {
          previousLayerId = layerId;
          h_layer_incident_lambda_true->fill(layerId, actualLambda, neutron->weight());
          h_layer_incident_lambda->fill(layerId, lambda_calculated, neutron->weight());
        }
        if(firstEnterToLayer[layerId]) {
          firstEnterToLayer[layerId] = 0;
          h_layer_first_incident_lambda_true->fill(layerId, actualLambda, neutron->weight());
          h_layer_first_incident_lambda->fill(layerId, lambda_calculated, neutron->weight());
        }
      }

      int previousStrawId = -1;
      //int firstEnterToLayer[4] = {1, 1, 1, 1};
      while (auto strawWallSegment = segments_StrawWall.next()) {
        const int locStrawId = strawWallSegment->volumeCopyNumber();
        const int tubeId = strawWallSegment->volumeCopyNumber(2);
        const int globalStrawId = tubeId * 7 + locStrawId;
        const double actualEkin = strawWallSegment->startEKin();

        auto stepF = strawWallSegment->firstStep();
        const double position[3] = {stepF->preGlobalX(), stepF->preGlobalY(), stepF->preGlobalZ()};
        const double lambda_calculated = calculateWavelength(initialPosition, position, stepF->preTime(), preGeant4Distance);
        const double actualLambda = Utils::neutronEKinToWavelength(actualEkin)/Units::angstrom;

        if(globalStrawId != previousStrawId && actualEkin) {
          previousStrawId = globalStrawId;
          h_straw_incident_lambda_true->fill(globalStrawId, actualLambda, neutron->weight());
          h_straw_incident_lambda->fill(globalStrawId, lambda_calculated, neutron->weight());
        }
        // if(firstEnterToLayer[layerId]) {
        //   firstEnterToLayer[layerId] = 0;
        //   h_layer_first_incident_lambda_true->fill(layerId, actualLambda, neutron->weight());
        //   h_layer_first_incident_lambda->fill(layerId, lambda_calculated, neutron->weight());
        // }
      }

      auto segL = neutron->lastSegment();
      if (segL->volumeName()=="Converter") {
        count_neutrons_converted += 1;

        auto stepL = segL->lastStep();
        const double position_conv[3] = {stepL->postGlobalX(), stepL->postGlobalY(), stepL->postGlobalZ()};

        /// volumeCopyNumber() = CountingGas; volumeCopyNumber(1) = Converter; volumeCopyNumber(2) = straw wall; volumeCopyNumber(3) = EmptyTube;
        /// volumeCopyNumber(4) = TubeWall; volumeCopyNumber(5) = EmptyPackBox; volumeCopyNumber(6) = Bank; volumeCopyNumber(7) = World
        const int strawId_conv = segL->volumeCopyNumber(1);
        const int tubeId_conv = segL->volumeCopyNumber(3);
        //const int bankId_conv = segL->volumeCopyNumber(5); only the rear bank is present in the geometry with bankId=0


        h_neutron_zy_conv->fill(position_conv[2]/Units::mm, position_conv[1]/Units::cm, neutron->weight());
        h_neutron_xy_conv->fill(-position_conv[0]/Units::mm, position_conv[1]/Units::mm, neutron->weight());

        const int layerNumber_conv = banks->getTubeLayerId(0, tubeId_conv, oldTubeNumbering);
        h_neutron_layerConvCounter->fill(layerNumber_conv, neutron->weight());

        if (hit.eventHasHit()) {
          count_neutrons_hit += 1;
          const double position_hit[3] = {hit.eventHitPositionX(), hit.eventHitPositionY(), hit.eventHitPositionZ()};

          h_neutron_zy_hit->fill(position_hit[2]/Units::mm, position_hit[1]/Units::cm, hit.eventHitWeight());
          h_neutron_xy_hit->fill(-position_hit[0]/Units::mm, position_hit[1]/Units::mm, hit.eventHitWeight());

          const double theta_hit = Utils::theta(hit.eventHitPosition())/Units::degree;
          h_neutron_theta_hit->fill(theta_hit, hit.eventHitWeight());

          const int pixelId = banks->getPixelId(0, tubeId_conv, strawId_conv, position_hit[0], position_hit[1]);

          h_neutron_pixel_hit->fill(pixelId%numberOfPixelsPerStraw, std::floor(pixelId/numberOfPixelsPerStraw), hit.eventHitWeight());
          h_neutron_layerHitCounter->fill(layerNumber_conv, hit.eventHitWeight());


          const double lambda_hit_calculated = calculateWavelength(initialPosition, position_hit, hit.eventHitTime(), preGeant4Distance);

          const double actualEkin = segL->startEKin();
          const double actualLambda = Utils::neutronEKinToWavelength(actualEkin) / Units::angstrom;

          if (layerNumber_conv == 0){
            h_layer0_lambda_hit->fill(lambda_hit_calculated, hit.eventHitWeight());
            h_layer0_lambda_true_hit->fill(actualLambda, hit.eventHitWeight());
          }
          else if(layerNumber_conv == 1){
            h_layer1_lambda_hit->fill(lambda_hit_calculated, hit.eventHitWeight());
            h_layer1_lambda_true_hit->fill(actualLambda, hit.eventHitWeight());
          }
          else if(layerNumber_conv == 2){
            h_layer2_lambda_hit->fill(lambda_hit_calculated, hit.eventHitWeight());
            h_layer2_lambda_true_hit->fill(actualLambda, hit.eventHitWeight());
          }
          else if(layerNumber_conv == 3){
            h_layer3_lambda_hit->fill(lambda_hit_calculated, hit.eventHitWeight());
            h_layer3_lambda_true_hit->fill(actualLambda, hit.eventHitWeight());
          }

          const int globalStrawId = tubeId_conv * 7 + strawId_conv;
          h_straw_lambda_hit->fill(globalStrawId, lambda_hit_calculated, hit.eventHitWeight());
          h_straw_lambda_true_hit->fill(globalStrawId, actualLambda, hit.eventHitWeight());

          if (createDetectionMcplFile == true) {
            detectionFile->addDetectionEvent(pixelId, hit.eventHitTime()/Units::ms);
          }
        } // end hit in event
      } // end if conversion condition
    }//end of loop over primary neutrons
  } //end of event loop

  delete detectionFile;
  hc.saveToFile("bcsloki_sans", true);
  return 0;
}
