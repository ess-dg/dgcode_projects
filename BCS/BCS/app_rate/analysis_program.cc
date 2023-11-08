#include "GriffDataRead/GriffDataReader.hh"
#include "GriffAnaUtils/TrackIterator.hh"
#include "GriffAnaUtils/TrackFilter_Primary.hh"
#include "GriffAnaUtils/TrackFilter_PDGCode.hh"
#include "GriffAnaUtils/SegmentIterator.hh"
#include "GriffAnaUtils/SegmentFilter_Volume.hh"
#include "SimpleHists/HistCollection.hh"
#include "GriffB10Common/DetHitApproximation.hh"
#include "Utils/ArrayMath.hh"
#include "Utils/NeutronMath.hh"

#include "Core/FPE.hh"
// #include "TMath.h"
// #include <TFile.h>
// #include <TH2Poly.h>
#include <cmath>
#include <iostream>

using namespace std;

// void AddCircleBin(TH2Poly *h2p, double x, double y, double *px, double *py) {

//   double r = 3.5;
//   double a = 0;
//   const int NP = 50; // Number of point to build the current bin
//   double da = 2*TMath::Pi()/NP; // Angle delta
//     for (int i = 0; i<NP; i++) {
//       a     = a+da;
//       px[i] = r*TMath::Cos(a) + x; // px[0];
//       py[i] = r*TMath::Sin(a) + y; //py[0];
//     }
//     h2p->AddBin(NP, px, py);
// }

int main (int argc, char**argv) {

  Core::catch_fpe();

  GriffDataReader dr(argc, argv);
  auto setup = dr.setup();
  auto &geo = setup->geo();
  if (geo.getName() != "G4GeoBCS/GeoBCS") {
    printf("Error: Wrong setup for this analysis\n");
    return 1;
  }
  setup->dump();

  DetHitApproximation hit(&dr, 1.2*Units::cm, 120*Units::keV, "CountingGas");

  //Filters:
  GriffAnaUtils::TrackIterator primary_neutrons(&dr);
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  GriffAnaUtils::SegmentIterator segments_StrawWall(&dr);
  segments_StrawWall.addFilter(new GriffAnaUtils::SegmentFilter_Volume("StrawWall"));
  segments_StrawWall.addFilter(new GriffAnaUtils::TrackFilter_Primary());

  GriffAnaUtils::SegmentIterator segments_TubeWall(&dr);
  segments_TubeWall.addFilter(new GriffAnaUtils::SegmentFilter_Volume("TubeWall"));
  segments_TubeWall.addFilter(new GriffAnaUtils::TrackFilter_Primary());

  GriffAnaUtils::SegmentIterator segments_Converter(&dr);
  segments_Converter.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Converter"));
  segments_Converter.addFilter(new GriffAnaUtils::TrackFilter_Primary());

  SimpleHists::HistCollection hc;

  int tubeCopyNum_min = 0;
  //int tubeCopyNum_max = 47;//Maximum tubecopyNumber is 46 in the simulation but the bin limits go like this:  [0-1[,[1-2[,...,[46,47[ //TODONOW 500
  //int tubeCopyNum_bins = 47; //TODONOW 500
  int tubeCopyNum_max = 500;//Maximum tubecopyNumber is 46 in the simulation but the bin limits go like this:  [0-1[,[1-2[,...,[46,47[ //TODONOW 500
  int tubeCopyNum_bins = 500; //TODONOW 500

  int strawCopyNum_min = 0;
  int strawCopyNum_max = 70;//Maximum strawcopyNumber is 60 in the simulation but the bin limits go like this:  [0-10[,[10-20[,...,[60,70[
  int strawCopyNum_bins = 7;


  auto h_st_peak1 = hc.book2D("h_st_peak1", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_peak1");
  h_st_peak1->setXLabel("straw copy number");
  h_st_peak1->setYLabel("tube copy number");
  auto h_st_peak2 = hc.book2D("h_st_peak2", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_peak2");
  h_st_peak2->setXLabel("straw copy number");
  h_st_peak2->setYLabel("tube copy number");

  auto h_st_peak_conv = hc.book2D("h_st_peak_conv", strawCopyNum_bins,strawCopyNum_min,strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_peak_conv");
  h_st_peak_conv->setXLabel("straw copy number");
  h_st_peak_conv->setYLabel("tube copy number");

  auto h_st_peak_hit = hc.book2D("h_st_peak_hit", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_peak_hit");
  h_st_peak_hit->setXLabel("straw copy number");
  h_st_peak_hit->setYLabel("tube copy number");

  auto h_st_peak_hit_firstPanelDifferentTof = hc.book2D("h_st_peak_hit_firstPanelDifferentTof", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, 3000, 0, 3000, "h_st_peak_hit_firstPanelDifferentTof");
  h_st_peak_hit_firstPanelDifferentTof->setXLabel("straw copy number");
  h_st_peak_hit_firstPanelDifferentTof->setYLabel("ToF & tube number for panel1");

  auto h_st_peak_div1 = hc.book2D("h_st_peak_div1",strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_peak_div1");
  h_st_peak_div1->setXLabel("straw copy number");
  h_st_peak_div1->setYLabel("tube copy number");
  auto h_st_peak_div2 = hc.book2D("h_st_peak_div2",strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_peak_div2");
  h_st_peak_div2->setXLabel("straw copy number");
  h_st_peak_div2->setYLabel("tube copy number");

  auto h_st_peak_div_conv = hc.book2D("h_st_peak_div_conv",strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_peak_div_conv");
  h_st_peak_div_conv->setXLabel("straw copy number");
  h_st_peak_div_conv->setYLabel("tube copy number");

  auto h_st_peak_div_hit = hc.book2D("h_st_peak_div_hit",strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_peak_div_hit");
  h_st_peak_div_hit->setXLabel("straw copy number");
  h_st_peak_div_hit->setYLabel("tube copy number");

  auto h_st_peak_div_hit_firstPanelDifferentTof = hc.book2D("h_st_peak_div_hit_firstPanelDifferentTof", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, 3000, 0, 3000, "h_st_peak_div_hit_firstPanelDifferentTof");
  h_st_peak_div_hit_firstPanelDifferentTof->setXLabel("straw copy number");
  h_st_peak_div_hit_firstPanelDifferentTof->setYLabel("ToF & tube number for panel1");


  auto h_st_global1 = hc.book2D("h_st_global1", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_global1");
  h_st_global1->setXLabel("straw copy number");
  h_st_global1->setYLabel("tube copy number");
  auto h_st_global2 = hc.book2D("h_st_global2", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_global2");
  h_st_global2->setXLabel("straw copy number");
  h_st_global2->setYLabel("tube copy number");

  auto h_st_global_conv = hc.book2D("h_st_global_conv", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_global_conv");
  h_st_global_conv->setXLabel("straw copy number");
  h_st_global_conv->setYLabel("tube copy number");

  auto h_st_global_hit = hc.book2D("h_st_global_hit", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_global_hit");
  h_st_global_hit->setXLabel("straw copy number");
  h_st_global_hit->setYLabel("tube copy number");


  auto h_new_local_average_incident = hc.book2D("h_new_local_average_incident [Hz]", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_new_local_average_incident");
  h_new_local_average_incident->setXLabel("straw copy number");
  h_new_local_average_incident->setYLabel("tube copy number");

  auto h_new_local_average_div_incident = hc.book2D("h_new_local_average_div_incident [Hz]", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_new_local_average_div_incident");
  h_new_local_average_div_incident->setXLabel("straw copy number");
  h_new_local_average_div_incident->setYLabel("tube copy number");

  auto h_new_local_peak_incident = hc.book2D("h_new_local_peak_inciden [kHz]t", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_new_local_peak_incident");
  h_new_local_peak_incident->setXLabel("straw copy number");
  h_new_local_peak_incident->setYLabel("tube copy number");

  auto h_new_local_peak_div_incident = hc.book2D("h_new_local_peak_div_incident [kHz]", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_new_local_peak_div_incident");
  h_new_local_peak_div_incident->setXLabel("straw copy number");
  h_new_local_peak_div_incident->setYLabel("tube copy number");


  auto h_local_peak_incident_firstPanelDifferentTof = hc.book2D("h_local_peak_incident_firstPanelDifferentTof [kHz]", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, 3000, 0, 3000, "h_local_peak_incident_firstPanelDifferentTof");
  h_local_peak_incident_firstPanelDifferentTof->setXLabel("straw copy number");
  h_local_peak_incident_firstPanelDifferentTof->setYLabel("ToF & tube number for panel1");

  auto h_local_peak_div_incident_firstPanelDifferentTof = hc.book2D("h_local_peak_div_incident_firstPanelDifferentTof [kHz]", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, 3000, 0, 3000, "h_local_peak_div_incident_firstPanelDifferentTof");
  h_local_peak_div_incident_firstPanelDifferentTof->setXLabel("straw copy number");
  h_local_peak_div_incident_firstPanelDifferentTof->setYLabel("ToF & tube number for panel1");

  //auto h_neutron_segment_number_StrawWall = hc.book1D("Number of StrawWall segments per neutron",50,0,50,"neutron_segment_number_StrawWall");
  //auto h_neutron_segment_length_StrawWall = hc.book1D("Sum length of StrawWall segments per neutron",300,0,3000,"neutron_segment_length_StrawWall");
  //h_neutron_segment_length_StrawWall->setXLabel("Sum pathlength [um]");

  auto h_panel_lambda = hc.book2D("Neutron wavelength for panels", 5, 0, 5, 100, 0, 15,"h_panel_lambda");
  h_panel_lambda->setXLabel("Panel number");
  h_panel_lambda->setYLabel("Wavelength [angstrom]");

  auto h_panel_lambda_hit = hc.book2D("Hit wavelength for panels", 5, 0, 5, 100, 0, 15, "h_panel_lambda_hit");
  h_panel_lambda_hit->setXLabel("Panel number");
  h_panel_lambda_hit->setYLabel("Wavelength [angstrom]");

  auto h_panel_tof = hc.book2D("Neutron TOF for panels", 5, 0, 5, 110, 0, 110, "h_panel_tof");
  h_panel_tof->setXLabel("Panel number");
  h_panel_tof->setYLabel("TOF [ms]");

  auto h_panel_tof_hit = hc.book2D("Hit TOF for panels", 5, 0, 5, 110, 0, 110, "h_panel_tof_hit");
  h_panel_tof_hit->setXLabel("Panel number");
  h_panel_tof_hit->setYLabel("TOF [ms]");

  auto h_neutron_tof_incident = hc.book1D("Incident neutrons TOF", 110, 0, 110, "neutron_tof");
  h_neutron_tof_incident->setXLabel("TOF [ms]");


  auto h_neutron_lambda_incident = hc.book1D("Incident neutrons lambda", 100, 0, 15, "incident_neutron_lambda");
  h_neutron_lambda_incident->setXLabel("Wavelength [angstrom]");

  auto h_neutron_true_lambda_detected = hc.book1D("Detected neutrons true lambda", 100, 0, 15, "detected_true_lambda");
  h_neutron_true_lambda_detected->setXLabel("Wavelength [angstrom]");

  auto h_neutron_true_lambda_converted = hc.book1D("Converted neutrons true lambda", 100, 0, 15, "converted_true_lambda");
  h_neutron_true_lambda_converted->setXLabel("Wavelength [angstrom]");


  auto h_neutron_tof_straw = hc.book1D("Straw neutrons TOF", 110, 0, 110, "neutron_tof_straw");
  h_neutron_tof_straw->setXLabel("TOF [ms]");

  auto h_neutron_tof_hit = hc.book1D("Hit neutrons TOF", 110, 0, 110, "neutron_tof_hit");
  h_neutron_tof_hit->setXLabel("TOF [ms]");

  float thetamax = 10.0;
  int thetabins = 500;
  auto h_neutron_theta_incident = hc.book1D("Neutron incident theta ", thetabins, 0, thetamax, "neutron_theta_incident");
  h_neutron_theta_incident->setXLabel("Theta [degree]");

  auto xhalf = -60; // 120 cm
  auto binsx = 240 * 5;
  auto yhalf = -60; // 120 cm
  auto binsy = 240 * 5;

  auto h_neutron_xy_incident = hc.book2D("Incident neutron xy [cm]", binsx, xhalf, -xhalf, binsy, yhalf, -yhalf, "neutron_xy_incident");
  h_neutron_xy_incident->setXLabel("x [cm]");
  h_neutron_xy_incident->setYLabel("y [cm]");

  auto h_neutron_xy_hit = hc.book2D("Hit neutron xy [cm]", binsx, xhalf, -xhalf, binsy, yhalf, -yhalf, "neutron_xy_hit");
  h_neutron_xy_hit->setXLabel("x [cm]");
  h_neutron_xy_hit->setYLabel("y [cm]");

  auto h_neutron_xy_div_hit = hc.book2D("Div hit neutron xy [cm]", binsx, xhalf, -xhalf, binsy, yhalf, -yhalf, "neutron_xy_div_hit");
  h_neutron_xy_div_hit->setXLabel("x [cm]");
  h_neutron_xy_div_hit->setYLabel("y [cm]");

  // auto h_zy = hc.book2D("zy",400,498,510,400,-10,10,"zy");
  // h_zy->setXLabel("global z (cm)");
  // h_zy->setYLabel("global y (cm)");

  // auto h_xy_div = hc.book2D("xy_div",300,-50,50,300,-50,50,"xy_div");
  // h_xy_div->setXLabel("global x (cm)");
  // h_xy_div->setYLabel("global y (cm)");

  auto h_neutron_xy_incident_secondPanel = hc.book2D("Panel 2 Incident neutron xy [cm]", binsx, xhalf, -xhalf, binsy, yhalf, -yhalf, "h_neutron_xy_incident_secondPanel");
  h_neutron_xy_incident_secondPanel->setXLabel("x [cm]");
  h_neutron_xy_incident_secondPanel->setYLabel("y [cm]");

  auto h_neutron_xy_incident_thirdPanel = hc.book2D("Panel 3 Incident neutron xy [cm]", binsx, xhalf, -xhalf, binsy, yhalf, -yhalf, "h_neutron_xy_incident_thirdPanel");
  h_neutron_xy_incident_thirdPanel->setXLabel("x [cm]");
  h_neutron_xy_incident_thirdPanel->setYLabel("y [cm]");

  auto h_neutron_xy_incident_fourthPanel = hc.book2D("Panel 4 Incident neutron xy [cm]", binsx, xhalf, -xhalf, binsy, yhalf, -yhalf, "h_neutron_xy_incident_fourthPanel");
  h_neutron_xy_incident_fourthPanel->setXLabel("x [cm]");
  h_neutron_xy_incident_fourthPanel->setYLabel("y [cm]");

  auto h_neutron_xy_incident_fifthPanel = hc.book2D("Panel 5 Incident neutron xy [cm]", binsx, xhalf, -xhalf, binsy, yhalf, -yhalf, "h_neutron_xy_incident_fifthPanel");
  h_neutron_xy_incident_fifthPanel->setXLabel("x [cm]");
  h_neutron_xy_incident_fifthPanel->setYLabel("y [cm]");

  auto h_panel_debugger = hc.book2D("Debugger for panel and tubes", 5, 0, 5, 40, 0, 40,"h_panel_debugger");
  h_panel_debugger->setXLabel("Panel number");
  h_panel_debugger->setYLabel("Tube number");

  auto h_panel_segment_length_Converter = hc.book2D("Sum pathlength of Converter segments per neutron for each panel", 5, 0, 5, 2000, 0, 200, "h_panel_segment_length_Converter");
  h_panel_segment_length_Converter->setXLabel("Panel number");
  h_panel_segment_length_Converter->setYLabel("Sum pathlength [um]");

  /*
    int binsz=520;
    float zmin=9980-5000; //mm
    float zmax=10110-5000; //mm
    int binsy=440;
    float ymin=-11; //cm
    float ymax=11; //cm
    float dz=(zmax-zmin)/binsz ;
    float dy=(ymax-ymin)/binsy ;

    auto h_neutron_zy_hit = hc.book2D("Dummy zy", binsz, zmin, zmax, binsy, ymin, -ymin,"my_Dummy_hist");
    h_neutron_zy_hit->setXLabel("z [mm]");
    h_neutron_zy_hit->setYLabel("y [cm]");
  */

  // TH2Poly *h2p = new TH2Poly();
  // h2p->SetName("test");
  // h2p->SetTitle("");
  // const int NP = 100; // Number of point to build the current bin
  // double px[NP];     // Bin's X positions
  // double py[NP];     // Bin's Y positions
  // AddCircleBin(h2p, 2,3, px,py);

  //int enteredArray [7][47]; //TODONOW 500
  int enteredArray [7][500]; //TODONOW 500

  ///int dummyArray [520][600];
  ///memset(dummyArray, 0, sizeof(dummyArray[0][0]) * binsz * binsy); //set all values zero
  // int myX,myY;
  //  std::fill( &enteredArray[0][0], &enteredArray[0][0] + sizeof(enteredArray), 0 );


  //const int minPeakTof[] = {33, 30, 26, 24, 22}; //config1
  //const int maxPeakTof[] = {36, 34, 28, 27, 26};
  //const int minPeakTof[] = {34, 32, 26, 26, 23}; //config1 try new
  //const int maxPeakTof[] = {35, 33, 27, 27, 24};

  //const int minPeakTof[] = {26, 26, 26, 26, 26}; //just for debugging
  //const int maxPeakTof[] = {27, 27, 27, 27, 27}; //just for debugging

  //config1
  //auto xDivergenceSquare = 4.21666666667 * 4.21666666667;
  //auto yDivergenceSquare = 3.75833333333 * 3.75833333333;
  //const int minPeakTof[] = {26, 26, 26, 26, 26}; //config2
  //const int maxPeakTof[] = {27, 27, 27, 27, 27};

  //config2
  //auto xDivergenceSquare = 2.2 * 2.2; //8.4
  //auto yDivergenceSquare = 2.2 * 2.2; //8.4
  //const int minPeakTof[] = {29, 29, 29, 29, 29}; //config2
  //const int maxPeakTof[] = {30, 30, 30, 30, 30};

  // config3
  auto xDivergenceSquare = 1.58125 * 1.58125; //
  auto yDivergenceSquare = 1.58125 * 1.58125; //
  const int minPeakTof[] = {26, 26, 26, 26, 26}; //
  const int maxPeakTof[] = {27, 27, 27, 27, 27};


  auto radiusDivergenceSquare = xDivergenceSquare + yDivergenceSquare;

  const double freqFactor[] = {(double) 1 / (14 * (maxPeakTof[0] - minPeakTof[0])), // one bin width is 1 ms -> result will be in [kHz] not [Hz]
                               (double) 1 / (14 * (maxPeakTof[1] - minPeakTof[1])),
                               (double) 1 / (14 * (maxPeakTof[2] - minPeakTof[2])),
                               (double) 1 / (14 * (maxPeakTof[3] - minPeakTof[3])),
                               (double) 1 / (14 * (maxPeakTof[4] - minPeakTof[4]))};

  //**************************************************************************************************************************************************************************//
  //*************************************************************************     All neutron      ***************************************************************************//
  //**************************************************************************************************************************************************************************//

  while (dr.loopEvents())
    {
      while (auto neutron = primary_neutrons.next())
        {
          //auto seg0 = neutron->segmentBegin();
          //auto step0 = seg0->firstStep();
          //const auto ekin_true = neutron->segmentBegin()->firstStep()->preEKin();
          //auto trueLambda = Utils::neutronEKinToWavelength(ekin_true) / Units::angstrom;

          const auto xIncident = neutron->segmentBegin()->firstStep()->postGlobalX()/Units::cm;
          const auto yIncident = neutron->segmentBegin()->firstStep()->postGlobalY()/Units::cm;
          const auto weigth_incident = neutron->weight();
          const auto ekinIncident = neutron->segmentBegin()->startEKin();
          const auto lambdaIncident = Utils::neutronEKinToWavelength(ekinIncident)/Units::angstrom;

          //const auto isOutOfDivergence = ((xIncident/Units::cm * xIncident/Units::cm) > xDivergenceSquare || (yIncident/Units::cm * yIncident/Units::cm) > yDivergenceSquare ); //rectangular
          const auto isOutOfDivergence = ( ((xIncident * xIncident) + (yIncident * yIncident))  > radiusDivergenceSquare ); //circle

          //if(0 <= trueLambda && trueLambda < 21.0) //What is that? or Why is it here?
          //  {
          //auto tofIncident = 0.0;

          //**************************************************************************************************************************************************************************//
          //**************************************************************     Tube segments (All neutron)      **********************************************************************//
          //**************************************************************************************************************************************************************************//

          auto previousPanelNr = 100; //anything except [0-6] would do it


          auto firstEnterToDetector = 1;
          auto firstEnterToSecondPanel = 1;
          auto firstEnterToThirdPanel = 1;
          auto firstEnterToFourthPanel = 1;
          auto firstEnterToFifthPanel = 1;


          // auto panel0 = 0;
          // auto panel1 = 0;
          // auto panel2 = 0;
          // auto panel3 = 0;
          // auto panel4 = 0;

          while (auto tubeWallSegment = segments_TubeWall.next()) // loop over all segments of the TubeWall
            {
              // there was a change in storing data to histograms when moving to 40 tubes instead of 7
              // the tube number [0,40[ and straw number goes from  [0,7[ so the unique ID is (tube number)*100 + straw number
              // instead of the old (tube number)*10 + straw number
              auto panelNr = (int)tubeWallSegment->volumeCopyNumber() / 100;

              if(firstEnterToDetector)
                {
                  h_neutron_xy_incident->fill(xIncident, yIncident, weigth_incident);

                  auto thetaIncident = Utils::theta(tubeWallSegment->firstStep()->preMomentumArray());
                  h_neutron_theta_incident->fill(thetaIncident/Units::deg, weigth_incident);

                  auto tofIncident = neutron->segmentBegin()->endTime();
                  h_neutron_tof_incident->fill(tofIncident/Units::ms, weigth_incident);

                  h_neutron_lambda_incident->fill(lambdaIncident, weigth_incident);

                  firstEnterToDetector = 0;
                }


              if( (panelNr == 1) && firstEnterToSecondPanel )
                {
                  h_neutron_xy_incident_secondPanel->fill(tubeWallSegment->firstStep()->preGlobalX()/Units::cm, tubeWallSegment->firstStep()->preGlobalY()/Units::cm, tubeWallSegment->getTrack()->weight());
                  firstEnterToSecondPanel = 0;
                }
              else if( (panelNr == 2) && firstEnterToThirdPanel )
                {
                  h_neutron_xy_incident_thirdPanel->fill(tubeWallSegment->firstStep()->preGlobalX()/Units::cm, tubeWallSegment->firstStep()->preGlobalY()/Units::cm, tubeWallSegment->getTrack()->weight());
                  firstEnterToThirdPanel = 0;
                }
              else if( (panelNr == 3) && firstEnterToFourthPanel )
                {
                  h_neutron_xy_incident_fourthPanel->fill(tubeWallSegment->firstStep()->preGlobalX()/Units::cm, tubeWallSegment->firstStep()->preGlobalY()/Units::cm, tubeWallSegment->getTrack()->weight());
                  firstEnterToFourthPanel = 0;
                }
              else if( (panelNr == 4) && firstEnterToFifthPanel )
                {
                  h_neutron_xy_incident_fifthPanel->fill(tubeWallSegment->firstStep()->preGlobalX()/Units::cm, tubeWallSegment->firstStep()->preGlobalY()/Units::cm, tubeWallSegment->getTrack()->weight());
                  firstEnterToFifthPanel = 0;
                }




              if(panelNr != previousPanelNr )// && tubeWallSegment->startEKin())
                {
                  previousPanelNr = panelNr;
                  auto actualEkin = tubeWallSegment->startEKin();
                  //is there a problem with the Ekin?????
                  auto actualLambda = Utils::neutronEKinToWavelength(actualEkin)/Units::angstrom;
                  h_panel_lambda->fill(panelNr, actualLambda, tubeWallSegment->getTrack()->weight());
                  //h_panel_lambda->fill(panelNr, actualLambda, neutron->weight());
                  h_panel_tof->fill(panelNr, tubeWallSegment->endTime()/Units::ms, tubeWallSegment->getTrack()->weight());
                  //h_panel_tof->fill(panelNr, tubeWallSegment->endTime()/Units::ms, neutron->weight());
                }




              // if( (panelNr == 0 && !panel0) || (panelNr == 1 && !panel1) || (panelNr == 2 && !panel2) || (panelNr == 3 && !panel3) || (panelNr == 4 && !panel4) )// && tubeWallSegment->startEKin())
              //   {
              //     switch (panelNr) {
              //     case 0:
              //       panel0 = 1;
              //       break;
              //     case 1:
              //       panel1 = 1;
              //       break;
              //     case 2:
              //       panel2 = 1;
              //       break;
              //     case 3:
              //       panel3 = 1;
              //       break;
              //     case 4:
              //       panel4 = 1;
              //       break;
              //     }
              //     //auto actualEkin = tubeWallSegment->startEKin();
              //     //is there a problem with the Ekin?????
              //     //auto actualLambda = Utils::neutronEKinToWavelength(actualEkin)/Units::angstrom;
              //     //h_panel_lambda->fill(panelNr, actualLambda, tubeWallSegment->getTrack()->weight());
              //     h_panel_lambda->fill(panelNr, 11, neutron->weight());
              //     //h_panel_tof->fill(panelNr, tubeWallSegment->endTime()/Units::ms, tubeWallSegment->getTrack()->weight());
              //     h_panel_tof->fill(panelNr, tubeWallSegment->endTime()/Units::ms, neutron->weight());

              //     //if(hit.eventHasHit()){
              //     //  h_panel_lambda->fill(panelNr, 8, neutron->weight());
              //     //}
              //   }
            }

          //**************************************************************************************************************************************************************************//
          //**************************************************************     Straw segments (All neutron)      *********************************************************************//
          //**************************************************************************************************************************************************************************//

          memset(enteredArray, 0, sizeof(enteredArray[0][0]) * strawCopyNum_bins * tubeCopyNum_bins); //set all values of enteredArray zero
          auto previousTubeNumberStrawSeg = 666; //anything except [0-39; 100-139; ...; 400-439] would do it
          auto previousStrawNumberStrawSeg = 666; //anything except [0,10,20,30,40,50,60] would do it

          while (auto segment = segments_StrawWall.next()) // loop over all segments of the StrawWall
            {
              auto tubeNumberStrawSeg = segment->volumeCopyNumber(2);
              auto strawNumberStrawSeg = segment->volumeCopyNumber();
              auto panelNumberStrawSeg = (int)tubeNumberStrawSeg / 100;

              auto weightStrawSeg = segment->getTrack()->weight();
              auto tofStrawSeg = segment->endTime();

              h_st_global1->fill(strawNumberStrawSeg, tubeNumberStrawSeg, weightStrawSeg);

              h_neutron_tof_straw->fill(tofStrawSeg/Units::ms, weightStrawSeg);
              if(tofStrawSeg/Units::ms >= minPeakTof[panelNumberStrawSeg] && tofStrawSeg/Units::ms <= maxPeakTof[panelNumberStrawSeg])
                {
                  h_st_peak1->fill(strawNumberStrawSeg, tubeNumberStrawSeg, weightStrawSeg * freqFactor[panelNumberStrawSeg]);
                  if(isOutOfDivergence)
                    {
                      h_st_peak_div1->fill(strawNumberStrawSeg, tubeNumberStrawSeg, weightStrawSeg * freqFactor[panelNumberStrawSeg]);
                    }
                }

              //note: division with 10 is only needed for the straw numbers(because it is 10,20...70, but tube number are OK [1,2,4..39; 100,101...139...;...400,401,...,439;
              if(enteredArray[(int)strawNumberStrawSeg / 10][tubeNumberStrawSeg])
                {
                  h_st_global2->fill(strawNumberStrawSeg, tubeNumberStrawSeg, weightStrawSeg);

                  if(tofStrawSeg/Units::ms >= minPeakTof[panelNumberStrawSeg] && tofStrawSeg/Units::ms <= maxPeakTof[panelNumberStrawSeg])
                    {
                      h_st_peak2->fill(strawNumberStrawSeg, tubeNumberStrawSeg, weightStrawSeg * freqFactor[panelNumberStrawSeg]);
                      if(isOutOfDivergence)
                        {
                          h_st_peak_div2->fill(strawNumberStrawSeg, tubeNumberStrawSeg, weightStrawSeg * freqFactor[panelNumberStrawSeg]);
                        }
                    }
                }
              else{
                enteredArray[strawNumberStrawSeg / 10][tubeNumberStrawSeg] = 1;
              }



              if( (tubeNumberStrawSeg != previousTubeNumberStrawSeg || strawNumberStrawSeg != previousStrawNumberStrawSeg)){ //only when entering a different straw than the previous one

                h_new_local_average_incident->fill(strawNumberStrawSeg, tubeNumberStrawSeg, weightStrawSeg);
                if(isOutOfDivergence){
                  h_new_local_average_div_incident->fill(strawNumberStrawSeg, tubeNumberStrawSeg, weightStrawSeg);
                }

                if(tofStrawSeg/Units::ms >= minPeakTof[panelNumberStrawSeg] && tofStrawSeg/Units::ms <= maxPeakTof[panelNumberStrawSeg]){
                  h_new_local_peak_incident->fill(strawNumberStrawSeg, tubeNumberStrawSeg, weightStrawSeg * freqFactor[panelNumberStrawSeg]);
                  if(isOutOfDivergence){
                    h_new_local_peak_div_incident->fill(strawNumberStrawSeg, tubeNumberStrawSeg, weightStrawSeg * freqFactor[panelNumberStrawSeg]);
                  }
                }

                //Rates for the first panel with 1ms TOF bins (treating 30 TOF bins virtually as 30 panels)
                if(panelNumberStrawSeg == 0 && tofStrawSeg/Units::ms >= 20 && tofStrawSeg/Units::ms < 50)
                  {
                    auto TofPanelNumber =  (int)(tofStrawSeg/Units::ms) - 20; // [0,30[
                    h_local_peak_incident_firstPanelDifferentTof->fill(strawNumberStrawSeg, tubeNumberStrawSeg + 100 * TofPanelNumber, weightStrawSeg * freqFactor[panelNumberStrawSeg]);
                    if(isOutOfDivergence)
                      {
                        h_local_peak_div_incident_firstPanelDifferentTof->fill(strawNumberStrawSeg, tubeNumberStrawSeg + 100 * TofPanelNumber, weightStrawSeg * freqFactor[panelNumberStrawSeg]);
                      }
                  }
                previousTubeNumberStrawSeg = tubeNumberStrawSeg;
                previousStrawNumberStrawSeg = strawNumberStrawSeg;
              }

            }

          //**************************************************************************************************************************************************************************//
          //**************************************************************     Converter segments (All neutron)      *********************************************************************//
          //**************************************************************************************************************************************************************************//

          //Variables to count to sum pathlenght in  of a neutron in
          auto panel0segmentLength = 0.0;
          auto panel1segmentLength = 0.0;
          auto panel2segmentLength = 0.0;
          auto panel3segmentLength = 0.0;
          auto panel4segmentLength = 0.0;

          while (auto segment = segments_Converter.next())// loop over all segments of the Converter
            {
              auto tubeNumberConvSeg = segment->volumeCopyNumber(3);
              auto panelNumberConvSeg = (int)tubeNumberConvSeg / 100;

              switch (panelNumberConvSeg) {
              case 0:
                panel0segmentLength += segment->segmentLength();
                break;
              case 1:
                panel1segmentLength += segment->segmentLength();
                break;
              case 2:
                panel2segmentLength += segment->segmentLength();
                break;
              case 3:
                panel3segmentLength += segment->segmentLength();
                break;
              case 4:
                panel4segmentLength += segment->segmentLength();
                break;
                //   default://DEBUGGING
                //panel4segmentLength += segment->segmentLength();
                //break;
              }
            }

          if(panel0segmentLength) h_panel_segment_length_Converter->fill(0, panel0segmentLength/Units::um, neutron->weight());
          if(panel1segmentLength) h_panel_segment_length_Converter->fill(1, panel1segmentLength/Units::um, neutron->weight());
          if(panel2segmentLength) h_panel_segment_length_Converter->fill(2, panel2segmentLength/Units::um, neutron->weight());
          if(panel3segmentLength) h_panel_segment_length_Converter->fill(3, panel3segmentLength/Units::um, neutron->weight());
          if(panel4segmentLength) h_panel_segment_length_Converter->fill(4, panel4segmentLength/Units::um, neutron->weight());


          //**************************************************************************************************************************************************************************//
          //**************************************************************     Converted neutrons      *******************************************************************************//
          //**************************************************************************************************************************************************************************//

          //int seg_count_StrawWall=0;
          //double seg_length_StrawWall=0;
          // double neutron_weight=0;

          auto segL = neutron->lastSegment();
          if(segL->volumeName()=="Converter") // check that the neutron ends in the converter, this is a condition equivalent to a converted neutron
            {
              // memset(enteredArray, 0, sizeof(enteredArray[0][0]) * strawCopyNum_bins * tubeCopyNum_bins); //set all values of enteredArray zero
              // segments_StrawWall.reset(); //reset SegmentIterator
              // while (auto segment = segments_StrawWall.next()){ // loop over all segments of the StrawWall /////////////////////////////
              //   someHistogram1_conv->fill(segment->volumeCopyNumber(),segment->volumeCopyNumber(2),neutron->weight());// segment->getTrack()->weight()); //should it be 'neutron->weight()'

              //   if(enteredArray[segment->volumeCopyNumber()/10][segment->volumeCopyNumber(2)]){
              //     someHistogram2_conv->fill(segment->volumeCopyNumber(),segment->volumeCopyNumber(2), neutron->weight());//segment->getTrack()->weight());
              //   }
              //   else{
              //     enteredArray[segment->volumeCopyNumber()/10][segment->volumeCopyNumber(2)]=1;;
              //   }
              // }


              // segments_StrawWall.reset(); //reset SegmentIterator
              // while (auto segment = segments_StrawWall.next()){ // loop over all segments of the StrawWall
              //   ++seg_count_StrawWall;
              //   seg_length_StrawWall+=segment->segmentLength();
              //   neutron_weight=segment->getTrack()->weight();
              // }
              // if(seg_length_StrawWall){
              //   h_neutron_segment_length_StrawWall->fill(seg_length_StrawWall/Units::um, neutron_weight);
              //   h_neutron_segment_number_StrawWall->fill(seg_count_StrawWall, neutron_weight);;
              // }

              auto tubeNumberConv = segL->volumeCopyNumber(3);
              auto strawNumberConv = segL->volumeCopyNumber(1);
              auto panelNumberConv = (int)tubeNumberConv / 100;
              auto weightConv = segL->getTrack()->weight();
              auto tofConv = segL->endTime();

              //=/=/=/=/=/=/=//=/=/=/=/=/=/=
              //auto actualLambda = Utils::neutronEKinToWavelength(actualEkin)/Units::angstrom;
              //h_panel_lambda->fill((int)(segL->volumeCopyNumber(3))/100, 9, neutron->weight());
              h_panel_debugger->fill(panelNumberConv, tubeNumberConv%100, 1);
              //=/=/=/=/=/=/==//=/=/=/=/=/=/

              h_neutron_tof_hit->fill(tofConv/Units::ms, weightConv);

              h_neutron_true_lambda_converted->fill(lambdaIncident, weightConv);

              h_st_global_conv->fill(strawNumberConv, tubeNumberConv, weightConv);//TEMPOFF
              if(tofConv/Units::ms >= minPeakTof[panelNumberConv] && tofConv/Units::ms <= maxPeakTof[panelNumberConv])
                {
                  h_st_peak_conv->fill(strawNumberConv, tubeNumberConv, weightConv * freqFactor[panelNumberConv]);
                  if(isOutOfDivergence)
                    {
                      h_st_peak_div_conv->fill(strawNumberConv, tubeNumberConv, weightConv * freqFactor[panelNumberConv]);
                    }
                }

              //**************************************************************************************************************************************************************************//
              //**************************************************************     Neutrons with hit     *********************************************************************************//
              //**************************************************************************************************************************************************************************//

              if(hit.eventHasHit())
                {


                  auto actualEkin = segL->startEKin();
                  auto actualLambda = Utils::neutronEKinToWavelength(actualEkin)/Units::angstrom;
                  h_panel_lambda_hit->fill(panelNumberConv, actualLambda, weightConv);

                  h_panel_tof_hit->fill(panelNumberConv, tofConv/Units::ms, weightConv);

                  h_st_global_hit->fill(strawNumberConv, tubeNumberConv, weightConv);

                  h_neutron_true_lambda_detected->fill(lambdaIncident, weightConv);

                  if(tofConv/Units::ms >= minPeakTof[panelNumberConv] && tofConv/Units::ms <= maxPeakTof[panelNumberConv])
                    {
                      h_st_peak_hit->fill(strawNumberConv, tubeNumberConv, weightConv * freqFactor[panelNumberConv]);
                      if(isOutOfDivergence)
                        {
                          h_st_peak_div_hit->fill(strawNumberConv, tubeNumberConv, weightConv * freqFactor[panelNumberConv]);
                        }
                    }

                  if(panelNumberConv == 0 && tofConv/Units::ms >= 20 && tofConv/Units::ms < 50)
                    {
                      auto TofPanelNumber =  (int)(tofConv/Units::ms) - 20; // [0,30[
                      h_st_peak_hit_firstPanelDifferentTof->fill(strawNumberConv, tubeNumberConv + 100 * TofPanelNumber, weightConv * freqFactor[panelNumberConv]);
                      if(isOutOfDivergence)
                        {
                          h_st_peak_div_hit_firstPanelDifferentTof->fill(strawNumberConv, tubeNumberConv + 100 * TofPanelNumber, weightConv * freqFactor[panelNumberConv]);
                        }
                    }

                  h_neutron_xy_hit->fill(hit.eventHitPositionX()/Units::cm, hit.eventHitPositionY()/Units::cm, weightConv);
                  if(isOutOfDivergence)
                    {
                      h_neutron_xy_div_hit->fill(hit.eventHitPositionX()/Units::cm, hit.eventHitPositionY()/Units::cm, weightConv);
                    }

                  /*
                     if(!dummyArray[(int)((hit.eventHitPositionZ()/Units::mm-zmin)/dz)][(int)((hit.eventHitPositionY()/Units::cm-ymin)/dy)]){//Make a Dummy Histogram for better visualisation
                     if(segL->volumeCopyNumber(1)+100*segL->volumeCopyNumber(3))
                     h_neutron_zy_hit->fill(hit.eventHitPositionZ()/Units::mm,hit.eventHitPositionY()/Units::cm,segL->volumeCopyNumber(1)+100*segL->volumeCopyNumber(3));
                     else
                     h_neutron_zy_hit->fill(hit.eventHitPositionZ()/Units::mm,hit.eventHitPositionY()/Units::cm,1);
                     dummyArray[(int)((hit.eventHitPositionZ()/Units::mm-zmin)/dz)][(int)((hit.eventHitPositionY()/Units::cm-ymin)/dy)]=1;
                     }
                  */

                  // there was a change in storing data to histograms when moving to 40 tubes instead of 7
                  // the tube number [0,40[ and straw number goes from  [0,7[ so the unique ID is (tube number)*100 + straw number
                  // instead of the old (tube number)*10 + straw number
                } //if(hit.eventHasHit()){

            }//if(segL->volumeName()=="Converter")

          //}
          //}
        }//while (auto neutron = primary_neutrons.next())
    }//while (dr.loopEvents())

  hc.saveToFile("rates_bcs.shist", true);
  // TFile f("bcs.root","RECREATE");
  // h2p->Write();
  // f.Close();

  return 0;
}
