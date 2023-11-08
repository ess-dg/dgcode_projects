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
#include <random>
#include "RandUtils/RandHelper.hh"

#ifndef M_PI
#define M_PI  3.14159265358979323846  //  pi
#endif

class signal
{
private:
  //static double dth_limit_min;
  static double dth_limit_max;
  //static double dphi_limit_min;
  static double dphi_limit_max;
  //static double dx_limit_min;
  static double dx_limit_max;
  //static double dy_limit_min;
  static double dy_limit_max;

public:
  //static double dtof_limit_min;
  static double dtof_limit_max; //setter function would be cleaner...
  //static double dlambda_limit_min;
  static double dlambda_limit_max;  //setter function would be cleaner...
  //static double dQ_limit_min;
  static double dQ_limit_max;

  static int isSignal_dth     (double dth)     {return (int) (-dth_limit_max     < dth     && dth     < dth_limit_max);}
  static int isSignal_dphi    (double dphi)    {return (int) (-dphi_limit_max    < dphi    && dphi    < dphi_limit_max);}
  static int isSignal_dx      (double dx)      {return (int) (-dx_limit_max      < dx      && dx      < dx_limit_max);}
  static int isSignal_dy      (double dy)      {return (int) (-dy_limit_max      < dy      && dy      < dy_limit_max);}
  static int isSignal_dtof    (double dtof)    {return (int) (-dtof_limit_max    < dtof    && dtof    < dtof_limit_max);}
  static int isSignal_dlambda (double dlambda) {return (int) (-dlambda_limit_max < dlambda && dlambda < dlambda_limit_max);}
  static int isSignal_dQ      (double dQ)      {return (int) (-dQ_limit_max      < dQ      && dQ      < dQ_limit_max);}
};

//double signal::dth_limit_min = 3 * -0.0345141928622;
double signal::dth_limit_max = 3 * 0.0345141928622;

//double signal::dphi_limit_min = 3 * -0.542488871037;
double signal::dphi_limit_max = 3 * 0.542488871037;

//double signal::dx_limit_min = 3 * -0.291508708298;
double signal::dx_limit_max = 3 * 0.291508708298;

//double signal::dy_limit_min = -0.3725; //outer radius of the converter layer (=inner radius of the staw) [cm]
double signal::dy_limit_max = 0.3725; // 3 * 0.347029558939;

//These 3 are wavelength dependent so set later
//double signal::dtof_limit_min = 3 * -0.000685914823486;
double signal::dtof_limit_max = 0.0;

//double signal::dlambda_limit_min = 3 * -0.001340988271;
double signal::dlambda_limit_max = 0.0;

//double signal::dQ_limit_min = 3 * -0.00210002220634;
double signal::dQ_limit_max = 0.0;



void setCoordinatesOfStrawCenterByTubeAndStrawCopyNumber(double* coordinateArrayToFill ,int tubeNumber, int strawNumber, GriffDataReader& dr){
  auto panelNumber = (int)tubeNumber / 100;
  tubeNumber = tubeNumber%100;
  strawNumber = (int)strawNumber/10;

  auto setup = dr.setup();

  const int numberOfTubes = setup->geo().getParameterInt("number_of_tubes");
  const double straw_radius = setup->geo().getParameterDouble("straw_radius_mm") *Units::mm;

  //position of straw centers in a tube before making a phisical volume van a ratation applied. [mm]
  double strawCenterBeforeRotation[7][2] = {{-2 * straw_radius, 0},
                                          {-straw_radius, -straw_radius * tan(M_PI/3)},
                                          {-straw_radius, straw_radius * tan(M_PI/3)},
                                          {0, 0},
                                          {straw_radius, -straw_radius * tan(M_PI/3)},
                                          {straw_radius, straw_radius * tan(M_PI/3)},
                                          {2 * straw_radius, 0}};

  //apply a rotation of alpha degrees on the center coordinates
  double strawCenter[7][2] = {{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}};
  auto alpha = 50 * (M_PI/180); //# 50 degree in radian
  for(auto i=0; i<=6; i++){
    strawCenter[i][0] = strawCenterBeforeRotation[i][0] * 1.0 * cos(alpha) - strawCenterBeforeRotation[i][1] * sin(alpha);
    strawCenter[i][1] = strawCenterBeforeRotation[i][0] * sin(alpha) + strawCenterBeforeRotation[i][1] * cos(alpha);
  }

    const double sample_detector_distance = setup->geo().getParameterDouble("sample_detector_distance_m") *Units::m;
    const double tube_box = setup->geo().getParameterDouble("tube_box_mm") *Units::mm;
    const double panel_offset = setup->geo().getParameterDouble("panel_offset_mm") *Units::mm;

    double sum_offset = 0.0;// #shift of the adjacent panels

    for(auto i=0; i<panelNumber; i++) {
      sum_offset += panel_offset;
      if(sum_offset >= tube_box){// # the offset cannot be larger than a tube_box
        sum_offset -= tube_box;// # we want the panels to cover each other
      }
    }


    auto yTubeCoord = (tubeNumber * tube_box - (numberOfTubes - 1) * 0.5 * tube_box + sum_offset);// #coordinates of the tube center
    auto zTubeCoord = sample_detector_distance + panelNumber * tube_box;


    //X coordinate is not the responsibility of this function because of the smearing along the wire
    //static double center[] = {0.0, yTubeCoord + strawCenter[strawNumber][0], zTubeCoord + strawCenter[strawNumber][1]};

    coordinateArrayToFill[1] = yTubeCoord + strawCenter[strawNumber][0];
    coordinateArrayToFill[2] = zTubeCoord + strawCenter[strawNumber][1];

//return center;
}

double smearCoordinateX(double xCoord, GriffDataReader &dr ){
  double fwhm = 6.0;//[mm]
  double standardDeviation = fwhm / 2.355;
  //double variance = standardDeviation * standardDeviation;

  ////static std::default_random_engine generator(time(0));
  // static std::default_random_engine generator(dr.seed());
  //std::normal_distribution<double> distribution(xCoord, fwhm);
  //return distribution(generator);

  //TK: Updated code below to use RandUtils::RandHelper instead of internal
  //NCrystal RandUtils (keeping dubious strategy on seeding via seed of *first*
  //GriffDataReader).
  static RandUtils::RandHelper rng{dr.seed()};

  return xCoord + rng.genNorm() * standardDeviation;
}


/// /// TODOGEANTINO
bool getPrimaryNeutronAndGeantino(GriffDataReader& dr,
                                                  const GriffDataRead::Track*& trk_neutron,
                                                  const GriffDataRead::Track*& trk_geantino)
{
  if (dr.nPrimaryTracks()!=2)
    return false;
  trk_neutron = dr.primaryTrackBegin();
  trk_geantino = trk_neutron +1;
  if (trk_neutron->pdgCode()!=2112)
    std::swap(trk_neutron,trk_geantino);
  if (trk_neutron->pdgCode()!=2112)
    return false;
  if (trk_geantino->pdgCode()!=999)
    return false;
  return true;
}
/// TRY IT


int main(int argc, char**argv) {

  Core::catch_fpe(); // framework tool fo catching floating point exceptions

  GriffDataReader dr(argc,argv); // the object that will read the griff file (output of simulation)
  auto setup = dr.setup();
  auto &geo = setup->geo();
  if (geo.getName() != "G4GeoBCS/GeoBCS") {
    printf("Error: Wrong setup for this analysis\n"); // check that you are running the correct analysis for the correct geometry
    return 1;
  }
  auto &gen = setup->gen();
  double source_to_sample_distance = 0*Units::m;
  // might be useful to check what generator you used, can access its parameters
  // or helps you define new ones, relevant to the analysis
  //if(gen.getName() == "G4MCPLPlugins.MCPLGen")
  //  source_to_sample_distance = 22.5*Units::m; //TODO heh?

  double genLambda = gen.getParameterDouble("neutron_wavelength_aangstrom");
  printf ("%f \n", genLambda);
  if(genLambda <= 1.0){//0.6
    signal::dtof_limit_max = 0.000565;
    signal::dlambda_limit_max = 0.00045;
    signal::dQ_limit_max = 3 * 0.0064933;
  }
  else if(1.0 < genLambda && genLambda <= 2.0){//1.8
    signal::dtof_limit_max = 0.001695;
    signal::dlambda_limit_max = 0.00134;
    signal::dQ_limit_max = 3 * 0.0021000;
  }
  else if(2.0 < genLambda && genLambda <= 4.0){//3
    signal::dtof_limit_max = 0.002825;
    signal::dlambda_limit_max = 0.00224;
    signal::dQ_limit_max = 3 * 0.0012500;
  }
  else if(4.0 < genLambda && genLambda <= 6.0){//5
    signal::dtof_limit_max = 0.004708;
    signal::dlambda_limit_max = 0.00373;
    signal::dQ_limit_max = 3 * 0.0007389;
  }
  else if(10.0 <= genLambda){//11
    signal::dtof_limit_max = 0.01063;
    signal::dlambda_limit_max = 0.00820;
    signal::dQ_limit_max = 3 * 0.0003338;
  }
  else{
    printf("Error: Problem with the genarator wavelength \n");
    return 1;
  }

  const double sample_detector_distance = setup->geo().getParameterDouble("sample_detector_distance_m") *Units::m;
  const double tube_radius = setup->geo().getParameterDouble("tube_radius_mm") *Units::mm;
  //auto &gen = setup->gen();
  setup->dump();

  DetHitApproximation hit(&dr, 1.2*Units::cm, 120*Units::keV, "CountingGas"); // define your hit approximation (detected neutron)

  GriffAnaUtils::SegmentIterator segments_gas(&dr); // track or segment filters that will make looping very efficient
  segments_gas.addFilter(new GriffAnaUtils::SegmentFilter_Volume("CountingGas"));

  GriffAnaUtils::TrackIterator primary_neutrons(&dr);
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  primary_neutrons.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));


  // GriffAnaUtils::SegmentIterator segments_World(&dr);
  // segments_World.addFilter(new GriffAnaUtils::SegmentFilter_Volume("World"));
  // segments_World.addFilter(new GriffAnaUtils::TrackFilter_Primary());

  GriffAnaUtils::SegmentIterator segments_TubeWall(&dr);
  segments_TubeWall.addFilter(new GriffAnaUtils::SegmentFilter_Volume("TubeWall"));
  segments_TubeWall.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  segments_TubeWall.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  // GriffAnaUtils::SegmentIterator segments_EmptyTube(&dr);
  // segments_EmptyTube.addFilter(new GriffAnaUtils::SegmentFilter_Volume("EmptyTube"));
  // segments_EmptyTube.addFilter(new GriffAnaUtils::TrackFilter_Primary());

  GriffAnaUtils::SegmentIterator segments_StrawWall(&dr);
  segments_StrawWall.addFilter(new GriffAnaUtils::SegmentFilter_Volume("StrawWall"));
  segments_StrawWall.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  segments_StrawWall.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  GriffAnaUtils::SegmentIterator segments_Converter(&dr);
  segments_Converter.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Converter"));
  segments_Converter.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  segments_Converter.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  GriffAnaUtils::SegmentIterator segments_CountingGas(&dr);
  segments_CountingGas.addFilter(new GriffAnaUtils::SegmentFilter_Volume("CountingGas"));
  segments_CountingGas.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  segments_CountingGas.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));

  SimpleHists::HistCollection hc; // define the object that will hold/collect all your histograms

  float xmin = -53;
  float ymin = -53; //20+1 tube in negative direction
  int tubeNr = 40 + 2; // 40 plus 1-1 empty space on both ends to make the plot nice
  int binsx = 1060;
  int binsy = 1060/2;
  int binsz = 160;
  int thetabins = 450;
  float trueThetaMax = 6; //8.5 if covering the ful panel
  float thetamax = 9.0;
  float dthetamin = -5; //-0.35
  float dthetamax = 5; //3
  int dthetabins = 1000;
  float zmin = (sample_detector_distance - tube_radius) - 20; //[mm]  //4980 for 5 m sd distance
  float zmax = (sample_detector_distance - tube_radius) + 140; //[mm //]5140 for 5 m sd distance
  int zbins = 160;


  // first string is the title of the histogram, last one how it will appear in the collection browser

  ////////incident neutron beam///////

  auto h_neutron_x = hc.book1D("Neutron true x (all events)", binsx/2, xmin, -xmin, "neutron_true_x");
       h_neutron_x->setXLabel("x [cm]");
  auto h_neutron_y = hc.book1D("Neutron true y (all events)", binsy/2, ymin, -ymin, "neutron_true_y");
       h_neutron_y->setXLabel("y [cm]");
  auto h_neutron_xy = hc.book2D("Neutron true xy (all events)", binsx/2, xmin, -xmin, binsy/2, ymin, -ymin, "neutron_true_xy");
       h_neutron_xy->setXLabel("x [cm]");
       h_neutron_xy->setYLabel("y [cm]");

  auto h_neutron_y_firstEnter =  hc.book1D("Neutron y first enter (all events)", tubeNr, ymin, -ymin, "neutron_y_firstEnter");
       h_neutron_xy->setXLabel("x [cm]");

  auto h_panel_neutronNr_all = hc.book1D("Number of neutrons entering each panel", 5, 0, 4, "panel_neutronNr_all");

  auto h_neutron_theta = hc.book1D("Neutron true theta (all events)", thetabins/2, 0, trueThetaMax, "neutron_true_theta");
       h_neutron_theta->setXLabel("Angle [degree]");
  auto h_neutron_true_phi = hc.book1D("Neutron true phi (all events)", 360,-180,180,"neutron_true_phi");
       h_neutron_true_phi->setXLabel("Angle [degree]");

  auto h_neutron_lambda = hc.book1D("Neutron true wavelength (all events)", 1300, 0, 13, "neutron_lambda");
       h_neutron_lambda->setXLabel("Wavelength [angstrom]");
  auto h_neutron_ekin = hc.book1D("Neutron true ekin (all events)", 500, 0, 100,"neutron_ekin");
  auto h_neutron_tof = hc.book1D("Neutron true tof at sample position(all events)", 100, 0, 1, "neutron_tof");
       h_neutron_tof->setXLabel("TOF [ms]");
  auto h_neutron_Q = hc.book1D("Neutron true Q (all events)", 350, 0, 0.35, "neutron_Q");


  //auto h_neutron_true_theta_enter_straw = hc.book1D("Neutron true theta (enter Straw)", thetabins/2, 0, thetamax, "neutron_true_theta_enter_straw");
  //     h_neutron_true_theta_enter_straw->setXLabel("Angle [degree]");

  auto h_panel_lambda = hc.book2D("Wavelength for panels", 5, 0, 5, 325, 0, 13, "panel_lambda");
       h_panel_lambda->setXLabel("Panel number");
       h_panel_lambda->setYLabel("Wavelength [angstrom]");

  ////////conversion///////////////
  auto h_neutron_x_conv = hc.book1D("Neutron x (converted)", binsx/2, xmin, -xmin, "neutron_x_conv");
       h_neutron_x_conv->setXLabel("x [cm]");
  auto h_neutron_dx_conv = hc.book1D("Neutron dx (conv) for detected neutrons", 800, xmin-10, -xmin+10, "neutron_dx_conv");
       h_neutron_dx_conv->setXLabel("dx [cm]");
  auto h_neutron_dx_conv_center = hc.book1D("Neutron dx (conv) Center for detected neutrons", 100, -0.01 , 0.01, "neutron_dx_conv_center");
       h_neutron_dx_conv_center->setXLabel("dx [cm]");
  auto h_neutron_y_conv = hc.book1D("Neutron y (converted)", binsy/2, ymin, -ymin, "neutron_y_conv");
       h_neutron_y_conv->setXLabel("y [mm]");
  auto h_neutron_dy_conv = hc.book1D("Neutron dy (conv) for detected neutrons", binsy/2, ymin, -ymin, "neutron_dy_conv");
       h_neutron_dy_conv->setXLabel("dy [cm]");
  auto h_neutron_z_conv = hc.book1D("Neutron z (converted)", binsz, zmin, zmax, "neutron_z_conv");
       h_neutron_z_conv->setXLabel("z [mm]");

  auto h_neutron_xy_conv = hc.book2D("Neutron xy (converted)", binsx/2, xmin, -xmin, binsy/2, ymin, -ymin, "neutron_xy_conv");
       h_neutron_xy_conv->setXLabel("x [cm]");
       h_neutron_xy_conv->setYLabel("y [cm]");
  auto h_neutron_zy_conv = hc.book2D("Neutron zy (converted)", binsz/2, zmin, zmax, binsy/2, ymin, -ymin, "neutron_yz_conv");
       h_neutron_zy_conv->setXLabel("z [mm]");
       h_neutron_zy_conv->setYLabel("y [cm]");
  auto h_neutron_zx_conv = hc.book2D("Neutron zx (converted)", binsz/2, zmin, zmax, binsx/2, xmin, -xmin, "neutron_zx_conv");
       h_neutron_zx_conv->setXLabel("z [mm]");
       h_neutron_zx_conv->setYLabel("x [cm]");

  auto h_neutron_phi_conv = hc.book1D("Neutron phi [degree] (conv)", 360, -180, 180, "neutron_phi_conv");
       h_neutron_phi_conv->setXLabel("Angle [degree]");
  auto h_neutron_dphi_conv = hc.book1D("Neutron dphi (conv)", 4*360, -180, 180, "neutron_dphi_conv");
       h_neutron_dphi_conv->setXLabel("dphi [degree]");
  auto h_neutron_theta_conv = hc.book1D("Neutron theta (converted)", thetabins/2, 0, thetamax, "neutron_theta_conv");
       h_neutron_theta_conv->setXLabel("Angle [degree]");
  //auto h_neutron_true_theta_conv = hc.book1D("Neutron true theta (converted)", thetabins/2, 0, trueThetaMax, "neutron_true_theta_conv");
  //     h_neutron_true_theta_conv->setXLabel("Angle [degree]");
  //auto h_neutron_true_phi_conv = hc.book1D("Neutron true phi (converted)", 360, -180, 180, "neutron_true_phi_conv");
  //     h_neutron_true_phi_conv->setXLabel("Angle [degree]");

  //auto h_neutron_phivstheta_conv = hc.book2D("Neutron phi vs theta [degree] (conv)", thetabins/2, 0, thetamax, 360, -180, 180, "neutron_phivstheta_conv");
  //     h_neutron_phivstheta_conv->setXLabel("theta [degree]");
  //     h_neutron_phivstheta_conv->setYLabel("phi [degree]");

  // auto h_neutron_dthetavstheta_conv = hc.book2D("Neutron dtheta vs theta [degree] (conv)", thetabins, 0, thetamax, dthetabins, dthetamin, dthetamax, "neutron_dthetavstheta_conv");
  //      h_neutron_dthetavstheta_conv->setXLabel("theta [degree]");
  //      h_neutron_dthetavstheta_conv->setYLabel("dtheta [degree]");

  // auto h_neutron_dthetavsthetatrue_conv = hc.book2D("Neutron dtheta vs theta true [degree] (conv)", thetabins, 0, trueThetaMax, dthetabins, dthetamin, dthetamax, "neutron_dthetavsthetatrue_conv");
  //      h_neutron_dthetavsthetatrue_conv->setXLabel("theta true [degree]");
  //      h_neutron_dthetavsthetatrue_conv->setYLabel("dtheta [degree]");

  //auto h_neutron_dthetavsz_conv = hc.book2D("Neutron dtheta vs z (converted)", zbins/2, zmin, zmax, dthetabins/10, dthetamin, dthetamax, "neutron_dthetavsz_conv");
  //     h_neutron_dthetavsz_conv->setXLabel("z [mm]");
  //     h_neutron_dthetavsz_conv->setYLabel("dtheta [degree]");

  auto h_neutron_dth_conv = hc.book1D("Neutron thetaConv - true theta (converted)", dthetabins/2, dthetamin, dthetamax, "neutron_dth_conv");
       h_neutron_dth_conv->setXLabel("dTheta [degree]");
  //auto h_neutron_dy_conv = hc.book1D("Neutron yConv - yTrue (converted)", 300, -10, 10, "neutron_dy_conv");

  auto h_neutron_ekin_conv = hc.book1D("Neutron true ekin (converted)", 500, 0, 100, "neutron_ekin_conv");
       h_neutron_ekin_conv->setXLabel("Energy [meV]");
  auto h_neutron_Q_conv = hc.book1D("Neutron Q (converted)", 500, 0, 0.20, "neutron_Q_conv");
  auto h_neutron_dQ_conv = hc.book1D("Neutron dQ (converted)", 1400, -0.9, 0.9, "neutron_dQ_conv");//kmilanTODO
  //auto h_neutron_dQoverQ_conv = hc.book1D("Neutron dQ/Q (converted)", 3000, -0.0000001, 0.0000001, "neutron_dQoverQ_conv");//kmilanTODO

  //auto h_neutron_dthVelocity_conv = hc.book1D("Neutron thetaVelocityConv - thetaVelicityInit (converted)", 200, -10, 190, "neutron_dthVelocity_conv");

  auto h_neutron_lambdaGeant_conv = hc.book1D("Neutron wavelength just before conversion (from Geant)", 130,0,13, "neutron_lambdaGeant_conv");//kmilanTODO
       h_neutron_lambdaGeant_conv->setXLabel("Wavelength [A]");

  auto h_neutron_tof_conv = hc.book1D("Neutron tof (conv)", 500, 2.0, 3.0, "neutron_tof_conv");
       h_neutron_tof_conv->setXLabel("TOF [ms]");

  auto h_neutron_tofProjected_conv = hc.book1D("Neutron tof calculated from Projected conv coord", 500, 2.0, 3, "neutron_tofProjected_conv");
       h_neutron_tofProjected_conv->setXLabel("TOF [ms]");

  auto h_neutron_dtof_conv = hc.book1D("Neutron dtof (tof_conv - (calculated from Projected conv coord))", 1000, -0.4, 0.5, "neutron_dtof_conv"); // -0.03, 0.17
       h_neutron_dtof_conv->setXLabel("dTOF [ms]");
  auto h_neutron_lambda_conv = hc.book1D("Neutron wavelength (calculated from tofConv)",130,0,13,"neutron_lambda_conv");//kmilanTODO
       h_neutron_lambda_conv->setXLabel("Wavelength [A]");
  auto h_neutron_dlambda_conv = hc.book1D("Neutron wavelength-true wavelength (converted)",1000, -0.3, 0.4,"neutron_dlambda_conv");//kmilanTODO
       h_neutron_dlambda_conv->setXLabel("Wavelength [A]");

  // auto h_neutron_dlambdaoverlambda_conv = hc.book1D("Neutron dlambda/lambda (converted)",3000,-0.01,0.05,"neutron_dlambdaoverlambda_conv"); //kmilanTODO


  /////////hit///////////////
  auto h_neutron_x_hit = hc.book1D("Neutron x (hit)", binsx/2, xmin, -xmin, "neutron_x_hit");
       h_neutron_x_hit->setXLabel("x [cm]");
  auto h_neutron_dx_hit = hc.book1D("Neutron dx (hit)", 800, xmin-10, -xmin+10, "neutron_dx_hit");
       h_neutron_dx_hit->setXLabel("dx [cm]");
  auto h_neutron_dx_noSmearing_hit = hc.book1D("Neutron dx without smearing (hit)", 800, xmin-10, -xmin+10, "neutron_dx_noSmearing_hit");
       h_neutron_dx_noSmearing_hit->setXLabel("dx [cm]");
  auto h_neutron_dx_mask1cm_hit = hc.book1D("Neutron dx with 1cm bins (hit)", 125, xmin-9.5, -xmin+9.5 , "neutron_dx_mask1cm_hit");
       h_neutron_dx_mask1cm_hit->setXLabel("dx [cm]");
  auto h_neutron_dx_mask1p75cm_hit = hc.book1D("Neutron dx with 1.75cm bins (hit)", 73, xmin-(10+1.75/2), -xmin+(10+1.75/2), "neutron_dx_mask1p75cm_hit");
       h_neutron_dx_mask1p75cm_hit->setXLabel("dx [cm]");


  //Auto h_neutron_true_x_hit = hc.book1D("Neutron true x [cm] (hit)", binsx, xmin, -xmin, "neutron_true_x_hit");
  auto h_neutron_y_hit = hc.book1D("Neutron y (hit)", binsy/2, ymin, -ymin, "neutron_y_hit");
       h_neutron_y_hit->setXLabel("y [cm]");
  auto h_neutron_dy_hit = hc.book1D("Neutron dy (hit)", 2*binsy, ymin, -ymin, "neutron_dy_hit");
       h_neutron_dy_hit->setXLabel("dy [cm]");
  auto h_neutron_dy_mask1cm_hit = hc.book1D("Neutron dy with 1cm bins (hit)", 125, ymin-9.5, -ymin+9.5, "neutron_dy_mask1cm_hit");
       h_neutron_dy_mask1cm_hit->setXLabel("dy [cm]");
  auto h_neutron_dy_mask0p75cm_hit = hc.book1D("Neutron dy with 0.75 cm bins (hit)", 169, -63-(0.75/2), 63+(0.75/2), "neutron_dy_mask0p75cm_hit");
       h_neutron_dy_mask0p75cm_hit->setXLabel("dy [cm]");
  auto h_neutron_dy_mask0p8cm_hit = hc.book1D("Neutron dy with 0.8 cm bins (hit)", 157, -62.4-(0.8/2), 62.4+(0.8/2), "neutron_dy_mask0p8cm_hit");
       h_neutron_dy_mask0p8cm_hit->setXLabel("dy [cm]");
  //auto h_neutron_true_y_hit = hc.book1D("Neutron true y [cm] (hit)",binsy, ymin, -ymin, "neutron_true_y_hit");
  auto h_neutron_z_hit = hc.book1D("Neutron z (hit)", 10*binsz, zmin,zmax, "neutron_z_hit");
       h_neutron_z_hit->setXLabel("z [mm]");

  auto h_neutron_z_hitMinusConv = hc.book1D("Neutron z (hit-conv)", 800, -4, 4, "neutron_z_hitMinusConv");
       h_neutron_z_hitMinusConv->setXLabel("z [mm]");
  //auto h_neutron_true_z_hit = hc.book1D("Neutron true z [mm] (hit)", binsz, zmin, zmax, "neutron_true_z_hit"); //If this is the position of the source than it should be around 0

 //auto h_neutron_y_hit_4efficiency =  hc.book1D("Neutron y for efficiency (hit)", tubeNr, ymin, -ymin, "neutron_y_hit_4efficiency");
  //     h_neutron_xy->setXLabel("x [cm]");

  auto h_panel_neutronNr_hit = hc.book1D("Number of hits in each panel", 5, 0, 4, "panel_neutronNr_hit");

  auto h_neutron_xy_hit = hc.book2D("Neutron xy (hit)", binsx/2, xmin, -xmin, binsy/2, ymin, -ymin,"neutron_xy_hit");
       h_neutron_xy_hit->setXLabel("x [cm]");
       h_neutron_xy_hit->setYLabel("y [cm]");

  auto h_neutron_xy_signalOnly_hit = hc.book2D("Neutron xy only if within signal limits (hit)", binsx/2, xmin, -xmin, binsy/2, ymin, -ymin,"neutron_xy_signalOnly_hit");
       h_neutron_xy_signalOnly_hit->setXLabel("x [cm]");
       h_neutron_xy_signalOnly_hit->setYLabel("y [cm]");

  // auto h_neutron_xy_trueDir_hit = hc.book2D("Neutron xy fromOriginalDirProjectedToHitPositionZ (hit)", binsx, xmin, -xmin, binsy, ymin, -ymin,"h_neutron_xy_trueDir_hit");
  //      h_neutron_xy_trueDir_hit->setXLabel("x [cm]");
  //      h_neutron_xy_trueDir_hit->setYLabel("y [cm]");

  auto h_neutron_zy_hit = hc.book2D("Neutron zy (hit)", binsz, zmin, zmax, binsy/2, ymin, -ymin, "neutron_yz_hit");
       h_neutron_zy_hit->setXLabel("z [mm]");
       h_neutron_zy_hit->setYLabel("y [cm]");
  auto h_neutron_zx_hit = hc.book2D("Neutron zx (hit)", binsz, zmin, zmax, binsx/2, xmin, -xmin, "neutron_zx_hit");
       h_neutron_zx_hit->setXLabel("z [mm]");
       h_neutron_zx_hit->setYLabel("x [cm]");


  auto h_neutron_theta_hit = hc.book1D("Neutron theta [degree] (hit)", thetabins/2, 0, thetamax, "neutron_theta_hit");
  //auto h_neutron_true_theta_hit = hc.book1D("Neutron true theta [degree] (hit)", thetabins/2, 0, trueThetaMax, "neutron_true_theta_hit");
  auto h_neutron_phi_hit = hc.book1D("Neutron phi  (hit)", 360, -180, 180, "neutron_phi_hit");
       h_neutron_phi_hit->setXLabel("phi [degree]");
  auto h_neutron_dphi_hit = hc.book1D("Neutron dphi (hit)", 4*360, -180, 180, "neutron_dphi_hit");
       h_neutron_dphi_hit->setXLabel("dphi [degree]");
  //auto h_neutron_dphi_hit_noSmearing = hc.book1D("Neutron dphi (hit_noSmearing)", 2*360, -180, 180, "neutron_dphi_hit_noSmearing");
  //     h_neutron_dphi_hit_noSmearing->setXLabel("dphi [degree]");
  auto h_neutron_true_phi_hit = hc.book1D("Neutron true phi (hit)", 360, -180, 180, "neutron_true_phi_hit");

  auto h_neutron_phivstheta_hit = hc.book2D("Neutron phi vs theta [degree] (hit)", thetabins/2, 0, thetamax, 360, -180, 180, "neutron_phivstheta_hit");
       h_neutron_phivstheta_hit->setXLabel("theta [degree]");
       h_neutron_phivstheta_hit->setYLabel("phi [degree]");

  //auto h_neutron_dthetavstheta_hit = hc.book2D("Neutron dtheta vs theta [degree] (hit)", thetabins/2, 0, thetamax, dthetabins/10, dthetamin, dthetamax, "neutron_dthetavstheta_hit");
  //     h_neutron_dthetavstheta_hit->setXLabel("theta [degree]");
  //     h_neutron_dthetavstheta_hit->setYLabel("dtheta [degree]");

  //auto h_neutron_dthetavsdLambda_hit = hc.book2D("Neutron dtheta vs dLambda (converted)", 800, -0.1, 0.3, dthetabins/10, dthetamin, dthetamax, "neutron_dthetavsdLambda_hit");
  //     h_neutron_dthetavsdLambda_hit->setXLabel("dlambda [A]");
  //     h_neutron_dthetavsdLambda_hit->setYLabel("dtheta [degree]");


  // auto h_neutron_dthetavsthetatrue_hit = hc.book2D("Neutron dtheta vs theta true [degree] (hit)", thetabins, 0, trueThetaMax, dthetabins, dthetamin, dthetamax, "neutron_dthetavsthetatrue_hit");
  //      h_neutron_dthetavsthetatrue_hit->setXLabel("theta true [degree]");
  //      h_neutron_dthetavsthetatrue_hit->setYLabel("dtheta [degree]");

  auto h_neutron_dthetavsz_hit = hc.book2D("Neutron dtheta vs z (hit)", zbins, zmin, zmax, dthetabins/10, dthetamin, dthetamax, "neutron_dthetavsz_hit");
       h_neutron_dthetavsz_hit->setXLabel("z [mm]");
       h_neutron_dthetavsz_hit->setYLabel("dtheta [degree]");

  auto h_neutron_dthetavsPanel_hit = hc.book2D("Neutron dtheta vs Panel (hit)",  5, 0, 4, dthetabins/10, dthetamin, dthetamax, "neutron_dthetavsPanel_hit");
       h_neutron_dthetavsPanel_hit->setXLabel("panel ");
       h_neutron_dthetavsPanel_hit->setYLabel("dtheta [degree]");

       //  auto h_neutron_dthetavsConvLambda_hit = hc.book2D("Neutron dtheta vs Wavelength at conversion point (hit)",  60, 0.5, 12.5, dthetabins, dthetamin, dthetamax, "h_neutron_dthetavsConvLambda_hit");
       //h_neutron_dthetavsConvLambda_hit->setXLabel("Wavelength [A]");
       //h_neutron_dthetavsConvLambda_hit->setYLabel("dtheta [degree]");

  auto h_neutron_dth_hit = hc.book1D("Neutron theta - true theta (hit)", dthetabins/2, dthetamin, dthetamax, "neutron_dth_hit");
       h_neutron_dth_hit->setXLabel("dTheta [degree]");
  auto h_neutron_dth_noSmearing_hit = hc.book1D("Neutron theta - true theta noSmearing (hit)", dthetabins/2, dthetamin, dthetamax, "neutron_dth_noSmearing_hit");
       h_neutron_dth_noSmearing_hit->setXLabel("dTheta [degree]");
  //auto h_neutron_dy_hit = hc.book1D("Neutron y - yTrue (hit)", 300, -10, 10, "neutron_dy_hit");

  //auto h_neutron_dthetaovertheta_hit = hc.book1D("Neutron dtheta/theta (hit)", 300, -100, 1.1, "neutron_dthetaovertheta_hit");//kmilanTODO
  auto h_neutron_tof_hit = hc.book1D("Neutron tof [ms] (hit)", 500, 2.0, 3, "neutron_tof_hit");
       h_neutron_tof_hit->setXLabel("TOF [ms]");

  auto h_neutron_tofProjected_hit = hc.book1D("Neutron tof calculated from Projected hit coord (hit)", 500, 2.0, 3, "neutron_tofProjected_hit");
       h_neutron_tofProjected_hit->setXLabel("TOF [ms]");

  auto h_neutron_dtof_hit = hc.book1D("Neutron dtof (tof_hit - (calculated from Projected hit coord)) (hit)", 1000, -0.4, 0.5, "neutron_dtof_hit"); //-0.03, 0.17
       h_neutron_dtof_hit->setXLabel("dTOF [ms]");

  //auto h_neutron_tofHit_minus_tofConv = hc.book1D("Neutron tof_hit - tof_conv)", 1000, -0.01, 0.01, "h_neutron_tofHit_minus_tofCnv");
  //     h_neutron_tofHit_minus_tofConv->setXLabel("dTOF [ms]");

  auto h_neutron_Q_hit = hc.book1D("Neutron Q [1/angstrom] (hit)", 200, 0, 0.4, "neutron_Q_hit");
       h_neutron_Q_hit->setXLabel("Q [1/angstrom]");
  auto h_neutron_Q_signalOnly_hit = hc.book1D("Neutron Q signalOnly [1/angstrom] (hit)", 200, 0, 0.4, "neutron_Q_signalOnly_hit");
       h_neutron_Q_signalOnly_hit->setXLabel("Q [1/angstrom]");
  auto h_neutron_dQ_hit = hc.book1D("Neutron dQ (hit)", 1400, -0.9, 0.9, "neutron_dQ_hit");
       h_neutron_dQ_hit->setXLabel("Q [1/angstrom]");


 //auto h_neutron_IQ_hit = hc.book1D("Neutron IQ [1/angstrom] (hit)", 300, 0, 0.3, "neutron_IQ_hit");
       //h_neutron_IQ_hit->setXLabel("Q [1/angstrom]");

  //auto h_neutron_dQ_noSmearing_hit = hc.book1D("Neutron dQ noSmearing (hit)", 800, -0.3, 0.3, "neutron_dQ_noSmearing_hit");
  //     h_neutron_dQ_noSmearing_hit->setXLabel("Q [1/angstrom]");
  //auto h_neutron_dQoverQ_hit = hc.book1D("Neutron dQ/Q (hit)", 300, -40, 2, "neutron_dQoverQ_hit");//kmilanTODO
  auto h_neutron_lambdaFromTof_hit = hc.book1D("Neutron wavelength calc. from TOF(hit)", 130, 0, 13, "neutron_lambdaFromTof_hit");//kmilanTODO
  auto h_neutron_dlambda_hit = hc.book1D("Neutron wavelength-true wavelength (hit)", 1000, -0.3, 0.4, "neutron_dlambda_hit");//kmilanTODO
       h_neutron_dlambda_hit->setXLabel("Wavelength [A]");
  //auto h_neutron_dlambda_noSmearing_hit = hc.book1D("Neutron wavelength-true wavelength noSmearing (hit)", 100, -0.1, 0.3, "neutron_dlambda_noSmearing_hit");//kmilanTODO
  //     h_neutron_dlambda_noSmearing_hit->setXLabel("Wavelength [A]");
  //auto h_neutron_dlambdaoverlambda_hit = hc.book1D("Neutron dlambda/lambda (hit)", 300, -0.01, 0.05, "neutron_dlambdaoverlambda_hit"); //kmilanTODO

  auto h_neutron_dekin_hit = hc.book1D("Neutron ekin true - detected (hit)", 300, -1, 200, "neutron_dekin_hit");
       h_neutron_dekin_hit->setXLabel("Energy [meV]");
  auto h_neutron_dekin_hit_center = hc.book1D("Neutron ekin true - detected (hit) Center renge", 150, -0.5, 1, "neutron_dekin_hit_center");
       h_neutron_dekin_hit_center->setXLabel("Energy [meV]");
  auto h_neutron_dekin_hit_center2 = hc.book1D("Neutron ekin true - detected (hit) Center renge", 100, -0.05, 0.05, "neutron_dekin_hit_center2");
       h_neutron_dekin_hit_center2->setXLabel("Energy [meV]");

  //auto h_neutron_x_FP = hc.book1D("Neutron x (first pass)",binsx,xmin,-xmin,"neutron_x_FP");
  //auto h_neutron_y_FP = hc.book1D("Neutron y (first pass)",binsy,ymin,-ymin,"neutron_y_FP");
  //auto h_neutron_theta_FP = hc.book1D("Neutron theta (first pass)",thetabins,0,thetamax,"neutron_theta_FP");
  //auto h_neutron_xy_FP = hc.book2D("Neutron xy (first pass)",binsx,xmin,-xmin,binsy,ymin,-ymin,"neutron_xy_FP");
  //h_neutron_xy_FP->setXLabel("x [cm]");
  //h_neutron_xy_FP->setYLabel("y [cm]");
  //auto h_neutron_Q_FP = hc.book1D("Neutron Q from first pass in gas",350,0,0.35,"neutron_Q_FP");
  //  auto h_dQ = hc.book1D("(Qi-Qf)/Qi",200,-2,2,"dQ");
  //auto h_neutron_ekin_FP = hc.book1D("Neutron ekin (first pass)",500,0,100,"neutron_ekin_FP");
  //auto h_neutron_theta_migration = hc.book2D("Neutron theta migration",200,0,20,200,0,20,"neutron_theta_migration");
  //auto h_neutron_deltaTheta = hc.book1D("Neutron dTheta",400,-5,5,"neutron_rel_dTheta");
  //auto h_neutron_tof_FP = hc.book1D("Neutron theta (first pass)",300,0,6,"neutron_theta_FP");
  //auto h_neutron_tof_conv = hc.book1D("Neutron theta (converted)",300,0,6,"neutron_theta_conv");

  auto h_edep = hc.book1D("Energy deposited in counting gas [keV]", 230, 0, 2300, "energyDeposition");
  auto h_neutron_nsegments = hc.book1D("Number of segments per neutron", 130, 0, 130, "neutron_nsegments");



  auto h_edep_alpha = hc.book1D("Energy deposited in counting gas by alpha [keV]", 230, 0, 2300, "energyDeposition_alpha");
  auto h_edep_lithium = hc.book1D("Energy deposited in counting gas by lithium [keV]", 230, 0, 2300, "energyDeposition_lithium");
  auto h_edep_alphalithium = hc.book1D("Energy deposited in counting gas by alpha or lithium [keV]", 230, 0, 2300, "energyDeposition_alphalithium");
  auto h_edep_electrons = hc.book1D("Energy deposited in counting gas by electrons [keV]", 230, 0, 2300, "energyDeposition_electrons");
  auto h_edep_gammas = hc.book1D("Energy deposited in counting gas by gammas [keV]", 230, 0, 2300, "energyDeposition_gammas");
  auto h_edep_else = hc.book1D("Energy deposited in counting gas by anything else [keV]", 230, 0, 2300, "energyDeposition_else");



  //auto h_neutron_tofExactHitPosition = hc.book1D("h_neutron_tofExactHitPosition", 500, 2.0, 3, "neutron_tofExactHitPosition");
  //auto h_neutron_dtofExactHitPosition = hc.book1D("h_neutron_dtofExactHitPosition", 100, -0.03, 0.17, "neutron_dtofExactHitPosition");

  //auto h_neutron_tofWireNoSmearPosition = hc.book1D("h_neutron_tofWireNoSmearPosition", 500, 2.0, 3, "neutron_tofWireNoSmearPosition");
  //auto h_neutron_dtofWireNoSmearPosition = hc.book1D("h_neutron_dtofWireNoSmearPosition", 100, -0.03, 0.17, "neutron_dtofWireNoSmearPosition");

  //auto h_neutron_tofWireWithSmearPosition = hc.book1D("h_neutron_tofWireWithSmearPosition", 500, 2.0, 3, "neutron_tofWireWithSmearPosition");
  auto h_neutron_dtofWireWithSmearPosition = hc.book1D("h_neutron_dtofWireWithSmearPosition", 1000, -0.4, 0.5, "neutron_dtofWireWithSmearPosition");//--------------------------------
  ////////////////segment counters///////////////

  //auto h_neutron_segment_number_World = hc.book1D("Number of World segments per neutron", 50, 0, 50, "neutron_segment_number_World");
  auto h_neutron_segment_number_TubeWall = hc.book1D("Number of TubeWall segments per neutron", 50, 0, 50, "neutron_segment_number_TubeWall");
  //auto h_neutron_segment_number_EmptyTube = hc.book1D("Number of EmptyTube segments per neutron", 50, 0, 50, "neutron_segment_number_EmptyTube");
  auto h_neutron_segment_number_StrawWall = hc.book1D("Number of StrawWall segments per neutron", 50, 0, 50, "neutron_segment_number_StrawWall");
  auto h_neutron_segment_number_Converter = hc.book1D("Number of Converter segments per neutron", 50, 0, 50, "neutron_segment_number_Converter");
  auto h_neutron_segment_number_CountingGas = hc.book1D("Number of CountingGas segments per neutron", 50, 0, 50,"neutron_segment_number_CountingGas");
  //TODO change them to COUNTER histogram

  //auto h_neutron_segment_length_World = hc.book1D("Sum length of World segments per neutron", 300, 30, 50, "neutron_segment_length_World");
  //     h_neutron_segment_length_World->setXLabel("Sum pathlength [cm]");
  auto h_neutron_segment_length_TubeWall = hc.book1D("Sum length of TubeWall segments per neutron", 600, 0, 100, "neutron_segment_length_TubeWall");
       h_neutron_segment_length_TubeWall->setXLabel("Sum pathlength [mm]");
  // auto h_neutron_segment_length_EmptyTube = hc.book1D("Sum length of EmptyTube segments per neutron", 120, 0, 60, "neutron_segment_length_EmptyTube");
  //     h_neutron_segment_length_EmptyTube->setXLabel("Sum pathlength [mm]");
  auto h_neutron_segment_length_StrawWall = hc.book1D("Sum length of StrawWall segments per neutron", 500, 0, 5, "neutron_segment_length_StrawWall");
       h_neutron_segment_length_StrawWall->setXLabel("Sum pathlength [mm]");
  auto h_neutron_segment_length_Converter = hc.book1D("Sum length of Converter segments per neutron", 200, 0, 0.1, "neutron_segment_length_Converter");
       h_neutron_segment_length_Converter->setXLabel("Sum pathlength [mm]");
  auto h_neutron_segment_length_CountingGas = hc.book1D("Sum length of CountingGas segments per neutron", 200, 0, 150, "neutron_segment_length_CountingGas");
       h_neutron_segment_length_CountingGas->setXLabel("Sum pathlength [mm]");


  auto h_counter = hc.bookCounts("Count detection events considered SIGNAL (within the 3 sigma) for diferent histograms","signal_counters"); //qualitative description
  auto countSignal_dth = h_counter->addCounter("dth");
  auto countSignal_dth_noSmear = h_counter->addCounter("dth_noSmear");
  auto countSignal_dlambda = h_counter->addCounter("dlambda");
  auto countSignal_dlambda_noSmear = h_counter->addCounter("dlambda_noSmear");
  auto countSignal_dlambdaANDdth = h_counter->addCounter("dlambdaANDdth");
  auto countSignal_dQ= h_counter->addCounter("dQ");
  auto countSignal_dQ_noSmear = h_counter->addCounter("dQ_noSmear");
  auto countSignal_dx = h_counter->addCounter("dx");
  auto countSignal_dx_noSmear = h_counter->addCounter("dx_noSmear");
  auto countSignal_dy = h_counter->addCounter("dy");
  auto countSignal_dtof = h_counter->addCounter("dtof");
  auto countSignal_dphi = h_counter->addCounter("dphi");


  auto h_counter_rest = hc.bookCounts("Count detection events considered SIGNAL (within the 3 sigma) for diferent histograms","signal_counters_rest"); //qualitative description
  auto countTestGeantino = h_counter_rest->addCounter("geantinoTest");
  auto countSignal_dlambdaConv = h_counter_rest->addCounter("dlambdaConv");
  auto countSignal_dtofConv = h_counter_rest->addCounter("dtofConv");
  auto countSignal_dQConv = h_counter_rest->addCounter("dQConv");

  // auto h_neutron_segment_number_loop = hc.book1D("Number of Converter segments per neutron with loop",50,0,50,"neutron_segment_number_loop");
  //auto h_neutron_segment_length_loop = hc.book1D("Sum length of Converter segments per neutron with loop",1000,0,100,"neutron_segment_length_loop");
  //h_neutron_segment_length_loop->setXLabel("Sum pathlength [um]");
  //////////===========================================================================================////

  int tubeCopyNum_min = 0;
  int tubeCopyNum_max = 500;//Maximum tubecopyNumber is 46 in the simulation but the bin limits go like this:  [0-1[,[1-2[,...,[46,47[ //TODONOW 500
  int tubeCopyNum_bins = 500; //TODONOW 500
  int strawCopyNum_min = 0;
  int strawCopyNum_max = 70;//Maximum strawcopyNumber is 60 in the simulation but the bin limits go like this:  [0-10[,[10-20[,...,[60,70[
  int strawCopyNum_bins = 7;

  auto h_st_countingGas_hit = hc.book2D("h_st_countingGas_hit", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_countingGas_hit");
       h_st_countingGas_hit->setXLabel("straw copy number");
       h_st_countingGas_hit->setYLabel("tube copy number");

  auto h_st_tubeWall_hit = hc.book1D("h_st_tubeWall_hit", tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_tubeWall_hit");
       h_st_tubeWall_hit->setXLabel("tube copy number");

  auto h_st_strawWall_hit = hc.book2D("h_st_strawWall_hit", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_st_strawWall_hit");
       h_st_strawWall_hit->setXLabel("straw copy number");
       h_st_strawWall_hit->setYLabel("tube copy number");


  const double x_panel_min = -52;
  const double x_panel_max = 52;
  const double y_panel_min = -52;
  const double y_panel_max = 52;


  /////////////////////////////////////////// Efficiency ////////////////////////////////////////////////////

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

  /// TODOGEANTINO
  auto h_geantino_global1 = hc.book2D("h_geantino_global1", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_geantino_global1");
       h_geantino_global1->setXLabel("straw copy number");
       h_geantino_global1->setYLabel("tube copy number");

  auto h_geantino_global2 = hc.book2D("h_geantino_global2", strawCopyNum_bins, strawCopyNum_min, strawCopyNum_max, tubeCopyNum_bins, tubeCopyNum_min, tubeCopyNum_max, "h_geantino_global2");
       h_geantino_global2->setXLabel("straw copy number");
       h_geantino_global2->setYLabel("tube copy number");
  // /// Todogeantino End

  int enteredArray [7][500];
  /// TODOGEANTINO
  int enteredArrayGeantino [7][500];

  //////////===========================================================================================////


  while (dr.loopEvents()) {



    ////////////////////////////////////////////////////// Segment Iterators //////////////////////////////////////////////////////////

    double edep_alpha = 0.0;
    double edep_lithium = 0.0;
    double edep_alphalithium = 0.0;
    double edep_electrons = 0.0;
    double edep_gammas = 0.0;
    double edep_else = 0.0;

    double edep = 0;
    double w = 1. ;
    while (auto segment = segments_gas.next()){ // loop over all segments of the counting gas  and collect the energy deposition of anything, that's how you get a distribution that should be close enough to a pulse height spectrum
      // under certain experimental conditions
      edep += segment->eDep();
      //assert((segment->getTrack()->weight()==neutron->weight())&&"different weights within segments of the same track");
      w = segment->getTrack()->weight();


      auto pdg = segment->getTrack()->pdgCode();
      if (pdg==1000020040/*alpha*/||pdg==1000030070/*Li7[0.0]*/) {
        edep_alphalithium += segment->eDep();

        if (pdg==1000020040){
          edep_alpha += segment->eDep();
        }
        else{//pdg==1000030070/*Li7[0.0]*/
          edep_lithium += segment->eDep();
        }
      }
      else if (pdg==11||pdg==-11) {
        edep_electrons += segment->eDep();
      }
      else if (pdg==22||pdg==-22) {
        edep_gammas += segment->eDep();
      }
      else{
        edep_else += segment->eDep();
      }
    }

    if(edep) h_edep->fill(edep/Units::keV,w);

    if(edep_alpha) h_edep_alpha->fill(edep_alpha/Units::keV,w);
    if(edep_lithium) h_edep_lithium->fill(edep_lithium/Units::keV,w);
    if(edep_alphalithium) h_edep_alphalithium->fill(edep_alphalithium/Units::keV,w);
    if(edep_electrons) h_edep_electrons->fill(edep_electrons/Units::keV,w);
    if(edep_gammas) h_edep_gammas->fill(edep_gammas/Units::keV,w);
    if(edep_else) h_edep_else->fill(edep_else/Units::keV,w);

    /////////////////////////////////////////////////////////// Loop over primary neutrons //////////////////////////////////////////////////////////

    // w and neutron->weight() are the same

    /// TODO comment while cycle (begin and end also)
    while (auto neutron = primary_neutrons.next()){ // loop over primary neutrons generated by the gun
      h_neutron_nsegments->fill(neutron->nSegments(), neutron->weight());

      // true means what the Monte Carlo generates, the ideal neutron
      auto seg0 = neutron->segmentBegin();
      auto step0 = seg0->firstStep();

      auto dir_true = step0->preMomentumArray();
      const double phi_true = Utils::phi(dir_true)/Units::deg;
      const double theta_true = Utils::theta(dir_true)/Units::deg;
      const double ekin_true = step0->preEKin();
      const double lambda_true = Utils::neutronEKinToWavelength(ekin_true)/Units::angstrom;

      // Q doesn't make sense for the FlexGen
      const double Q_true = 4 * M_PI * sin(0.5 * theta_true *Units::deg) / Utils::neutronEKinToWavelength(ekin_true);

      const double positionExact_source[] = {step0->preGlobalX(), step0->preGlobalY(), step0->preGlobalZ()};

      h_neutron_x->fill(positionExact_source[0]/Units::cm, neutron->weight());
      h_neutron_y->fill(positionExact_source[1]/Units::cm, neutron->weight());
      //h_neutron_z->fill(step0->preGlobalZ()/Units::cm, neutron->weight());

      h_neutron_xy->fill(positionExact_source[0]/Units::cm, positionExact_source[1]/Units::cm, neutron->weight());

      const double xIncident = step0->postGlobalX()/Units::cm;
      const double yIncident = step0->postGlobalY()/Units::cm;
      if(xIncident >= x_panel_min && xIncident <= x_panel_max &&  yIncident >= y_panel_min && yIncident <= y_panel_max){
        h_neutron_y_firstEnter->fill(yIncident, neutron->weight());
      }

      h_neutron_theta->fill(theta_true, neutron->weight());//TODO true?
      h_neutron_true_phi->fill(phi_true, neutron->weight());

      h_neutron_lambda->fill(lambda_true, neutron->weight());
      h_neutron_ekin->fill(ekin_true/Units::meV, neutron->weight());
      h_neutron_tof->fill(seg0->startTime()/Units::ms, neutron->weight()); //step0->preTime()
      h_neutron_Q->fill(Q_true*Units::angstrom, neutron->weight());


      /// TODOGEANTINO
      const GriffDataRead::Track* trk_neutron;
      const GriffDataRead::Track* trk_geantino;
      //if (!getPrimaryNeutronAndGeantino(dr,trk_neutron,trk_geantino)) {
      if (getPrimaryNeutronAndGeantino(dr,trk_neutron,trk_geantino)) {

        countTestGeantino += 1;

        memset(enteredArrayGeantino, 0, sizeof(enteredArray[0][0]) * strawCopyNum_bins * tubeCopyNum_bins); //set all values of enteredArray zero

        for (auto seg = trk_geantino->segmentBegin(); seg!=trk_geantino->segmentEnd(); ++seg){

          auto geantino_weight=seg->getTrack()->weight();
          if (seg->volumeName()=="StrawWall"){

            auto tubeNumberStrawSegGeantino = seg->volumeCopyNumber(2);
            auto strawNumberStrawSegGeantino = seg->volumeCopyNumber();

            h_geantino_global1->fill(strawNumberStrawSegGeantino, tubeNumberStrawSegGeantino, geantino_weight);

            if(enteredArrayGeantino[(int)strawNumberStrawSegGeantino / 10][tubeNumberStrawSegGeantino])
              {
                h_geantino_global2->fill(strawNumberStrawSegGeantino, tubeNumberStrawSegGeantino, geantino_weight);
              }
            else{
              enteredArrayGeantino[strawNumberStrawSegGeantino / 10][tubeNumberStrawSegGeantino] = 1;
            }

          }
        }
        //}
      }
      /// TODOGEANTINO END



      // loop over neutron segments, until you hit the first
      // segment in the StrawWall (first pass=FP)
      //auto segCuEnd = neutron->segmentEnd();
      //auto segCu = seg0;
      //for (++segCu; segCu < segCuEnd; ++segCu) {
      //  if (segCu->volumeName() == "StrawWall"){
      //    h_neutron_true_theta_enter_straw->fill(theta_true, neutron->weight());
      //    break;
      //  }
      //}
      //if (segCu==segCuEnd) {
      //  //Neutron never reaches gas
      //   continue;
      // }

      auto previousPanelNr = 10; //anything except [0-6] would do it
      auto panelNr = 0;
      auto actualEkin = 0.0;
      auto actualLambda = 0.0;
      //auto panelStart = 5000; //Start of the first panal in Z direction [mm] //TODO RU SURE?
      //auto panelWidth = 25.4; // length of a panel in Z direction [mm]

      while (auto tubeWallSegment = segments_TubeWall.next()){ // loop over all segments of the TubeWall

        ///////
        //panelNr = (int)((tubeWallSegment->firstStep()->preGlobalZ()/Units::mm - panelStart) / panelWidth) ; // AND do the "same" for HIT. create 2 1DHists.
        panelNr = (int)tubeWallSegment->volumeCopyNumber() / 100;

        if(panelNr != previousPanelNr && tubeWallSegment->startEKin()){
          previousPanelNr = panelNr;

          h_panel_neutronNr_all->fill(panelNr, tubeWallSegment->getTrack()->weight());

          actualEkin=tubeWallSegment->startEKin();
          //is there a problem with the Ekin?????
          actualLambda = Utils::neutronEKinToWavelength(actualEkin)/Units::angstrom;
          h_panel_lambda->fill(panelNr, actualLambda, tubeWallSegment->getTrack()->weight());
        }
        /////////
      }
      segments_TubeWall.reset();//it is needed


      memset(enteredArray, 0, sizeof(enteredArray[0][0]) * strawCopyNum_bins * tubeCopyNum_bins); //set all values of enteredArray zero
      int previousTubeNumberStrawSeg = 666; //anything except [0-39; 100-139; ...; 400-439] would do it
      int previousStrawNumberStrawSeg = 666; //anything except [0,10,20,30,40,50,60] would do it

      while (auto segment = segments_StrawWall.next()){ // loop over all segments of the StrawWall

        auto neutron_weight=segment->getTrack()->weight();

        auto tubeNumberStrawSeg = segment->volumeCopyNumber(2);
        auto strawNumberStrawSeg = segment->volumeCopyNumber();
        h_st_global1->fill(strawNumberStrawSeg, tubeNumberStrawSeg, neutron_weight);

        //note: division with 10 is only needed for the straw numbers(because it is 10,20...70, but tube number are OK [1,2,4..39; 100,101...139...;...400,401,...,439;
        if(enteredArray[(int)strawNumberStrawSeg / 10][tubeNumberStrawSeg])
          {
            h_st_global2->fill(strawNumberStrawSeg, tubeNumberStrawSeg, neutron_weight);
          }
        else{
          enteredArray[strawNumberStrawSeg / 10][tubeNumberStrawSeg] = 1;
        }

        if( (tubeNumberStrawSeg != previousTubeNumberStrawSeg || strawNumberStrawSeg != previousStrawNumberStrawSeg)){ //only when entering a different straw than the previous one
          h_new_local_average_incident->fill(strawNumberStrawSeg, tubeNumberStrawSeg, neutron_weight);
        }
        previousTubeNumberStrawSeg = tubeNumberStrawSeg;
        previousStrawNumberStrawSeg = strawNumberStrawSeg;

      }
      segments_StrawWall.reset();//it is needed


      //auto step1 = seg0->firstStep();
      //////////////////////////////////////////////////////////////////////TEMP TRY FOR ABSORPTION/////////////////////////////////////////////////////////////////////////////////////////////////
      double neutron_weight = 0;
      int seg_count_TubeWall = 0;
      int seg_count_StrawWall = 0;
      double seg_length_TubeWall = 0;
      double seg_length_StrawWall = 0;


      while (auto segment = segments_TubeWall.next()){ // loop over all segments of the TubeWall
          ++seg_count_TubeWall;
          seg_length_TubeWall+=segment->segmentLength();
          neutron_weight=segment->getTrack()->weight();
        }
        if(seg_length_TubeWall){
          h_neutron_segment_length_TubeWall->fill(seg_length_TubeWall/Units::mm, neutron_weight);
          h_neutron_segment_number_TubeWall->fill(seg_count_TubeWall, neutron_weight);;
        }

        while (auto segment = segments_StrawWall.next()){ // loop over all segments of the StrawWall
          ++seg_count_StrawWall;
          seg_length_StrawWall+=segment->segmentLength();
          neutron_weight=segment->getTrack()->weight();
        }
        if(seg_length_StrawWall){
          h_neutron_segment_length_StrawWall->fill(seg_length_StrawWall/Units::mm, neutron_weight);
          h_neutron_segment_number_StrawWall->fill(seg_count_StrawWall, neutron_weight);;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      auto segL = neutron->lastSegment();
      if(segL->volumeName() == "Converter") { // check that the neutron ends in the converter, this is a condition equivalent to a converted neutron

        //segment counters
        //int seg_count_World = 0;
//        int seg_count_TubeWall = 0;
        //int seg_count_EmptyTube = 0;
//       int seg_count_StrawWall = 0;
        int seg_count_Converter = 0;
        int seg_count_CountingGas = 0;

        //double seg_length_World = 0; //Units::cm;
//        double seg_length_TubeWall = 0;
        //double seg_length_EmptyTube = 0;
//        double seg_length_StrawWall = 0;
        double seg_length_Converter = 0;
        double seg_length_CountingGas = 0;



//        double neutron_weight = 0;
        // auto segL = neutron->lastSegment();
        //if(segL->volumeName()=="Converter") { // check that the neutron ends in the converter, this is a condition equivalent to a converted neutron


        // while (auto segment = segments_World.next()){ // loop over all segments of the World
        //   ++seg_count_World;
        //   seg_length_World+=segment->segmentLength();
        //   neutron_weight=segment->getTrack()->weight();
        // }
        // if(seg_length_World){
        //   h_neutron_segment_length_World->fill(seg_length_World/Units::cm, neutron_weight);
        //   h_neutron_segment_number_World->fill(seg_count_World, neutron_weight);
        // }
 /*
        while (auto segment = segments_TubeWall.next()){ // loop over all segments of the TubeWall
          ++seg_count_TubeWall;
          seg_length_TubeWall+=segment->segmentLength();
          neutron_weight=segment->getTrack()->weight();
        }
        if(seg_length_TubeWall){
          h_neutron_segment_length_TubeWall->fill(seg_length_TubeWall/Units::mm, neutron_weight);
          h_neutron_segment_number_TubeWall->fill(seg_count_TubeWall, neutron_weight);;
        }

        // while (auto segment = segments_EmptyTube.next()){ // loop over all segments of the EmptyTube
        //   ++seg_count_EmptyTube;
        //   seg_length_EmptyTube+=segment->segmentLength();
        //   neutron_weight=segment->getTrack()->weight();
        // }
        // if(seg_length_EmptyTube){
        //   h_neutron_segment_length_EmptyTube->fill(seg_length_EmptyTube/Units::mm, neutron_weight);
        //   h_neutron_segment_number_EmptyTube->fill(seg_count_EmptyTube, neutron_weight);;
        // }


        while (auto segment = segments_StrawWall.next()){ // loop over all segments of the StrawWall
          ++seg_count_StrawWall;
          seg_length_StrawWall+=segment->segmentLength();
          neutron_weight=segment->getTrack()->weight();
        }
        if(seg_length_StrawWall){
          h_neutron_segment_length_StrawWall->fill(seg_length_StrawWall/Units::mm, neutron_weight);
          h_neutron_segment_number_StrawWall->fill(seg_count_StrawWall, neutron_weight);;
        }
 */
        while (auto segment = segments_Converter.next()){ // loop over all segments of the Converter
          ++seg_count_Converter;
          seg_length_Converter+=segment->segmentLength();
          neutron_weight=segment->getTrack()->weight();
        }
        if(seg_length_Converter){
          h_neutron_segment_length_Converter->fill(seg_length_Converter/Units::mm, neutron_weight);
          h_neutron_segment_number_Converter->fill(seg_count_Converter, neutron_weight);;
        }

        while (auto segment = segments_CountingGas.next()){ // loop over all segments of the CountingGas
          ++seg_count_CountingGas;
          seg_length_CountingGas+=segment->segmentLength();
          neutron_weight=segment->getTrack()->weight();
        }
        if(seg_length_CountingGas){
          h_neutron_segment_length_CountingGas->fill(seg_length_CountingGas/Units::mm, neutron_weight);
          h_neutron_segment_number_CountingGas->fill(seg_count_CountingGas, neutron_weight);;
        }
        //}

        /////////////////////////////////////////////////////////// End of Segment Iterators //////////////////////////////////////////////////////////




        auto stepL = segL->lastStep();
        assert(step0 && stepL && "step info required");

        double dir_conv[3];
        Utils::subtract(stepL->postGlobalArray(), step0->preGlobalArray(), dir_conv);

        //preMomentumArrray
        //double dirVelocity_conv[3];
        //Utils::subtract(stepL->postGlobalArray(), stepL->preGlobalArray(), dirVelocity_conv);

        //auto lastStepDir_cov = stepL->preMomentumArray();
        //const double thetaVelocityDir_conv = Utils::theta(lastStepDir_cov);
        //h_neutron_dthVelocity_conv->fill(thetaVelocityDir_conv/Units::deg - theta_true , neutron->weight());


        const double theta_conv = Utils::theta(dir_conv);
        const double phi_conv = Utils::phi(dir_conv)/Units::deg;
        //const double theta_conv = Utils::theta(stepL->postGlobalArray());
        const double ekin_conv = stepL->preEKin(); // postEKin() returns 0 since the neutron stops
        const double lambdaGeant_conv =  Utils::neutronEKinToWavelength(ekin_conv)/Units::angstrom;
        const double Q_conv = 4 * M_PI * sin(0.5 * theta_conv) / Utils::neutronEKinToWavelength(ekin_conv);
        const double tof_conv = stepL->postTime()/Units::ms;

        //const double positionExact_conv[] = {stepL->postGlobalX(), stepL->postGlobalY(), stepL->postGlobalZ()};
        const int tubeNumberConv = segL->volumeCopyNumber(3);
        const int strawNumberConv = segL->volumeCopyNumber(1);
        double position_conv[3] = {stepL->postGlobalX(), stepL->postGlobalY(), stepL->postGlobalZ()};

        const double rho_conv = position_conv[2]/Units::cm * tan(theta_true *Units::deg);
        //const double rho = positionOnWire_hit[2]/Units::cm * tan(theta_true *Units::deg);
        const double xFromOriginalDirProjectedToConvPositionZ = rho_conv * cos(phi_true *Units::deg) + positionExact_source[0]/Units::cm; //[cm]
        const double yFromOriginalDirProjectedToConvPositionZ = rho_conv * sin(phi_true *Units::deg) + positionExact_source[1]/Units::cm; //[cm]


        h_neutron_x_conv->fill(position_conv[0]/Units::cm, neutron->weight());
        h_neutron_y_conv->fill(position_conv[1]/Units::cm, neutron->weight());
        h_neutron_z_conv->fill(position_conv[2]/Units::mm, neutron->weight());

        h_neutron_xy_conv->fill(position_conv[0]/Units::cm, position_conv[1]/Units::cm, neutron->weight());
        h_neutron_zy_conv->fill(position_conv[2]/Units::mm, position_conv[1]/Units::cm, neutron->weight());
        h_neutron_zx_conv->fill(position_conv[2]/Units::mm, position_conv[0]/Units::cm, neutron->weight());
        /*
        h_neutron_xy_conv->fill(stepL->postGlobalX()/Units::cm, stepL->postGlobalY()/Units::cm, neutron->weight());
        h_neutron_zy_conv->fill(stepL->postGlobalZ()/Units::mm, stepL->postGlobalY()/Units::cm, neutron->weight());
        h_neutron_zx_conv->fill(stepL->postGlobalZ()/Units::mm, stepL->postGlobalX()/Units::cm, neutron->weight());
        */
        h_neutron_theta_conv->fill(theta_conv/Units::deg, neutron->weight());
        //h_neutron_true_theta_conv->fill(theta_true, neutron->weight());
        h_neutron_phi_conv->fill(phi_conv, neutron->weight());
        //h_neutron_true_phi_conv->fill(phi_true, neutron->weight());
        double dphi_conv = phi_conv - phi_true; //NOT const because it needs correction
        dphi_conv += ((dphi_conv<-180)-(dphi_conv>180)) * 360; //correct dphi_conv because of the 'periodic boundary condition' (0==360)
        h_neutron_dphi_conv->fill(dphi_conv, neutron->weight());

        //h_neutron_phivstheta_conv->fill(theta_conv/Units::deg, phi_conv, neutron->weight());
        //h_neutron_dthetavsthetatrue_conv->fill(theta_true, theta_conv/Units::deg - theta_true, neutron->weight());
        //h_neutron_dthetavstheta_conv->fill(theta_conv/Units::deg, theta_conv/Units::deg - theta_true, neutron->weight());
        //h_neutron_dthetavsz_conv->fill(position_conv[2]/Units::mm, theta_conv/Units::deg - theta_true, neutron->weight());

        h_neutron_dth_conv->fill(theta_conv/Units::deg - theta_true, neutron->weight());
        //h_neutron_dy_conv->fill(position_conv[1]/Units::cm - positionExact_source[1]/Units::cm, neutron->weight()); //TODO This is WRONG

        h_neutron_ekin_conv->fill(ekin_conv/Units::meV, neutron->weight());
        h_neutron_lambdaGeant_conv->fill(lambdaGeant_conv, neutron->weight());

        h_neutron_Q_conv->fill(Q_conv*Units::angstrom, neutron->weight());
        const double dQ_conv = Q_conv*Units::angstrom - Q_true*Units::angstrom;
        h_neutron_dQ_conv->fill(dQ_conv, neutron->weight());
        countSignal_dQConv += signal::isSignal_dQ(dQ_conv);
        //h_neutron_dQoverQ_conv->fill((Q_conv*Units::angstrom - Q_true*Units::angstrom) / Q_conv*Units::angstrom, neutron->weight());





        //const double rho_conv = position_conv[2] *Units::mm * tan(theta_true *Units::deg);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //const double xFromOriginalDirProjectedToConvPositionZ = rho_conv * cos(phi_true *Units::deg);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //const double yFromOriginalDirProjectedToConvPositionZ = rho_conv * sin(phi_true *Units::deg);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        //const double sourceToProjectedConvPositionDistance = sqrt( xFromOriginalDirProjectedToConvPositionZ *Units::mm * xFromOriginalDirProjectedToConvPositionZ *Units::mm +
        //                                                           yFromOriginalDirProjectedToConvPositionZ *Units::mm * yFromOriginalDirProjectedToConvPositionZ *Units::mm +
        //                                                           position_conv[2] *Units::mm * position_conv[2] *Units::mm );//[mm]
        const double sourceToProjectedConvPositionDistance = position_conv[2] *Units::mm / cos(theta_true *Units::deg);
        const double TofProjetcedConv = sourceToProjectedConvPositionDistance/*[mm]*/ / Utils::neutron_angstrom_to_meters_per_second(lambda_true) /*[m/s]*/; //[ms]


        h_neutron_tof_conv->fill(tof_conv, hit.eventHitWeight());
        h_neutron_tofProjected_conv->fill(TofProjetcedConv, neutron->weight());
        const double dtof_conv = tof_conv - TofProjetcedConv;
        h_neutron_dtof_conv->fill(dtof_conv, neutron->weight());
        countSignal_dtofConv += signal::isSignal_dtof(dtof_conv);




        const double sample_to_convPoint_distance = sqrt( position_conv[0] *Units::mm * position_conv[0] *Units::mm +
                                                          position_conv[1] *Units::mm * position_conv[1] *Units::mm +
                                                          position_conv[2] *Units::mm * position_conv[2] *Units::mm );//[mm]
        const double velocity_conv = ((source_to_sample_distance + sample_to_convPoint_distance)/Units::m) / (tof_conv *Units::ms /Units::second);
        const double lambda_conv = Utils::neutron_meters_per_second_to_angstrom(velocity_conv);

        h_neutron_lambda_conv->fill(lambda_conv, neutron->weight());
        const double dlambda_conv = lambda_conv - lambda_true;
        h_neutron_dlambda_conv->fill(dlambda_conv, neutron->weight());
        countSignal_dlambdaConv += signal::isSignal_dlambda(dlambda_conv);

        h_st_global_conv->fill(strawNumberConv, tubeNumberConv, neutron->weight());//TEMPOFF

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////// HIT ////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(hit.eventHasHit()) { // check that there is a hit in the event, i.e. a neutron has been detected

          const double dx_conv = position_conv[0]/Units::cm - xFromOriginalDirProjectedToConvPositionZ;
          h_neutron_dx_conv->fill(dx_conv, neutron->weight());
          h_neutron_dx_conv_center->fill(dx_conv, neutron->weight());
          const double dy_conv = position_conv[1]/Units::cm - yFromOriginalDirProjectedToConvPositionZ;
          h_neutron_dy_conv->fill(dy_conv, neutron->weight());

          //const double positionExact_hit[] = {stepL->postGlobalX(), stepL->postGlobalY(), stepL->postGlobalZ()};

          double positionOnWire_hit[3];
          setCoordinatesOfStrawCenterByTubeAndStrawCopyNumber(positionOnWire_hit, tubeNumberConv, strawNumberConv, dr); //tube and straw number is the same for conversion and hit (copy constructor could be used also..)
          //double positionOnWire_hit[3] = {stepL->postGlobalX(), stepL->postGlobalY(), stepL->postGlobalZ()}; //TEMP for debugging

          positionOnWire_hit[0] = smearCoordinateX(hit.eventHitPositionX(), dr);

          //const double generatorToExactHitPositionDistance = sqrt( hit.eventHitPositionX() * hit.eventHitPositionX() + //only works for pointSource at origo USE step0->preGlobalArray()
          //                                                    hit.eventHitPositionY() * hit.eventHitPositionY() +
          //                                                    hit.eventHitPositionZ() * hit.eventHitPositionZ() );

          const double sampleToWireNoSmearPositionDistance = sqrt( hit.eventHitPositionX() * hit.eventHitPositionX() + //only works for pointSource at origo
                                                                   positionOnWire_hit[1] * positionOnWire_hit[1] +
                                                                   positionOnWire_hit[2] * positionOnWire_hit[2] );

          const double sampleToWireWithSmearPositionDistance = sqrt( positionOnWire_hit[0] * positionOnWire_hit[0] + //only works for pointSource at origo
                                                                     positionOnWire_hit[1] * positionOnWire_hit[1] +
                                                                     positionOnWire_hit[2] * positionOnWire_hit[2] );

          //const double theta_hit = Utils::theta(hit.eventHitPosition());
          const double theta_hit = Utils::theta(positionOnWire_hit)/Units::deg;
          const double theta_hit_noSmearing = Utils::theta(hit.eventHitPosition())/Units::deg;

          const double phi_hit = Utils::phi(positionOnWire_hit)/Units::deg;
          const double phi_hit_noSmearing = Utils::phi(hit.eventHitPosition())/Units::deg;
          const double tof_hit = hit.eventHitTime()/Units::ms;

          auto panelNumber_conv = (int)tubeNumberConv / 100; //panel number is the same for conversion and hit

          double velocity = -1;
          double velocityNoSmearing = -1;
          if(tof_hit){
            velocity = ((sampleToWireWithSmearPositionDistance + source_to_sample_distance)/Units::m) / (hit.eventHitTime()/Units::second); //source_to_sample_distance=0
            velocityNoSmearing = ((sampleToWireNoSmearPositionDistance + source_to_sample_distance)/Units::m) / (hit.eventHitTime()/Units::second);
          }
          else{
              printf("Error in hit tof value, tof zero or negative \n");
              return 1;
            }

          const double lambdaHit = Utils::neutron_meters_per_second_to_angstrom(velocity);
          const double lambdaHitNoSmearing = Utils::neutron_meters_per_second_to_angstrom(velocityNoSmearing);
          const double Q_hit = 4 * M_PI * sin(0.5 * theta_hit*Units::deg) / lambdaHit;
          const double Q_hit_noSmearing = 4 * M_PI * sin(0.5 * theta_hit_noSmearing*Units::deg) / lambdaHitNoSmearing;

          const double rho = hit.eventHitPositionZ()/Units::cm * tan(theta_true *Units::deg);
          //const double rho = positionOnWire_hit[2]/Units::cm * tan(theta_true *Units::deg);
          const double xFromOriginalDirProjectedToHitPositionZ = rho * cos(phi_true *Units::deg) + positionExact_source[0]/Units::cm; //[cm]
          const double yFromOriginalDirProjectedToHitPositionZ = rho * sin(phi_true *Units::deg) + positionExact_source[1]/Units::cm; //[cm]

          h_neutron_x_hit->fill(positionOnWire_hit[0]/Units::cm, hit.eventHitWeight());
          const double dx_hit = positionOnWire_hit[0]/Units::cm - xFromOriginalDirProjectedToHitPositionZ;
          h_neutron_dx_hit->fill(dx_hit, hit.eventHitWeight()); //debug reproduce h_neutron_x_hit
          h_neutron_dx_mask1cm_hit->fill(dx_hit, hit.eventHitWeight());
          h_neutron_dx_mask1p75cm_hit->fill(dx_hit, hit.eventHitWeight());

          countSignal_dx += signal::isSignal_dx(dx_hit);
          //h_neutron_true_x_hit->fill(positionExact_source[0]/Units::cm, hit.eventHitWeight());
          h_neutron_dx_noSmearing_hit->fill(hit.eventHitPositionX()/Units::cm - xFromOriginalDirProjectedToHitPositionZ, hit.eventHitWeight()); //debug reproduce h_neutron_x_hit
          countSignal_dx_noSmear += signal::isSignal_dx(hit.eventHitPositionX()/Units::cm - xFromOriginalDirProjectedToHitPositionZ);

          h_neutron_y_hit->fill(positionOnWire_hit[1]/Units::cm, hit.eventHitWeight());
          const double dy_hit = positionOnWire_hit[1]/Units::cm - yFromOriginalDirProjectedToHitPositionZ;

          h_neutron_dy_hit->fill(dy_hit, hit.eventHitWeight()); //debug reproduce h_neutron_y_hit
          h_neutron_dy_mask1cm_hit->fill(dy_hit, hit.eventHitWeight());
          h_neutron_dy_mask0p75cm_hit->fill(dy_hit, hit.eventHitWeight());
          h_neutron_dy_mask0p8cm_hit->fill(dy_hit, hit.eventHitWeight());
          countSignal_dy += signal::isSignal_dy(dy_hit);
          //h_neutron_true_y_hit->fill(positionExact_source[1]/Units::cm, hit.eventHitWeight());

          h_neutron_z_hit->fill(positionOnWire_hit[2]/Units::mm, hit.eventHitWeight());
          //h_neutron_true_z_hit->fill(positionExact_source[2]/Units::mm, hit.eventHitWeight());

          h_neutron_z_hitMinusConv->fill(hit.eventHitPositionZ()/Units::mm - position_conv[2]/Units::mm, hit.eventHitWeight());//---------------------------------------------------------------------------------

          //h_neutron_y_hit_4efficiency->fill(positionOnWire_hit[1]/Units::cm, hit.eventHitWeight());

          h_neutron_xy_hit->fill(positionOnWire_hit[0]/Units::cm, positionOnWire_hit[1]/Units::cm, hit.eventHitWeight());

          if(signal::isSignal_dx(dx_hit) && signal::isSignal_dy(dy_hit)) {
            h_neutron_xy_signalOnly_hit->fill(positionOnWire_hit[0]/Units::cm, positionOnWire_hit[1]/Units::cm, hit.eventHitWeight());
            h_neutron_Q_signalOnly_hit->fill(Q_hit, hit.eventHitWeight());
          }



          //h_neutron_xy_trueDir_hit->fill(xFromOriginalDirProjectedToHitPositionZ, yFromOriginalDirProjectedToHitPositionZ, hit.eventHitWeight());
          h_neutron_zy_hit->fill(positionOnWire_hit[2]/Units::mm, positionOnWire_hit[1]/Units::cm, hit.eventHitWeight());
          h_neutron_zx_hit->fill(positionOnWire_hit[2]/Units::mm, positionOnWire_hit[0]/Units::cm, hit.eventHitWeight());

          h_panel_neutronNr_hit->fill(panelNumber_conv, hit.eventHitWeight());

          h_neutron_theta_hit->fill(theta_hit, hit.eventHitWeight());
          //h_neutron_true_theta_hit->fill(theta_true, hit.eventHitWeight());
          const double dth_hit = theta_hit - theta_true;
          h_neutron_dth_hit->fill(dth_hit, hit.eventHitWeight());
          countSignal_dth += signal::isSignal_dth(dth_hit);
          h_neutron_dth_noSmearing_hit->fill(theta_hit_noSmearing  - theta_true, hit.eventHitWeight());
          countSignal_dth_noSmear += signal::isSignal_dth(theta_hit_noSmearing  - theta_true);

          h_neutron_phi_hit->fill(phi_hit, hit.eventHitWeight());
          h_neutron_true_phi_hit->fill(phi_true, hit.eventHitWeight());
          double dphi_hit = phi_hit - phi_true; //NOT const because it needs correction
          dphi_hit += ((dphi_hit<-180)-(dphi_hit>180)) * 360; //correct dphi_hit because of the 'periodic boundary condition' (0==360)
          h_neutron_dphi_hit->fill(dphi_hit, hit.eventHitWeight());
          countSignal_dphi += signal::isSignal_dphi(dphi_hit);

          double dphi_hit_noSmearing = phi_hit_noSmearing - phi_true; //NOT const because it needs correction
          dphi_hit_noSmearing += ((dphi_hit_noSmearing<-180)-(dphi_hit_noSmearing>180)) * 360; //correct dphi_hit because of the 'periodic boundary condition' (0==360)
          //h_neutron_dphi_hit_noSmearing->fill(dphi_hit_noSmearing, hit.eventHitWeight());


          h_neutron_phivstheta_hit->fill(theta_hit, phi_hit, hit.eventHitWeight());
          //h_neutron_dthetavstheta_hit->fill(theta_hit, theta_hit-theta_true, hit.eventHitWeight());
          //h_neutron_dthetavsthetatrue_hit->fill(theta_true, theta_hit-theta_true, hit.eventHitWeight());
          h_neutron_dthetavsz_hit->fill(positionOnWire_hit[2]/Units::mm, theta_hit-theta_true, hit.eventHitWeight());

          h_neutron_dthetavsPanel_hit->fill(panelNumber_conv, theta_hit-theta_true, hit.eventHitWeight());
          //h_neutron_dthetavsConvLambda_hit->fill(lambda_conv, theta_hit-theta_true, hit.eventHitWeight());




          h_neutron_tof_hit->fill(tof_hit, hit.eventHitWeight());
          //     const double sourceToProjectedHitPositionDistance = sqrt( xFromOriginalDirProjectedToHitPositionZ *Units::cm * xFromOriginalDirProjectedToHitPositionZ *Units::cm +
          //                                                          yFromOriginalDirProjectedToHitPositionZ *Units::cm * yFromOriginalDirProjectedToHitPositionZ *Units::cm +
          //                                                          positionOnWire_hit[2] *Units::mm * positionOnWire_hit[2] *Units::mm );//[mm]
          const double sourceToProjectedHitPositionDistance = hit.eventHitPositionZ() / cos(theta_true *Units::deg); // positionOnWire_hit[2] / cos(theta_true *Units::deg);
          const double tofProjetcedHit = sourceToProjectedHitPositionDistance/*[mm]*/ / Utils::neutron_angstrom_to_meters_per_second(lambda_true) /*[m/s]*/; //[ms]
          h_neutron_tofProjected_hit->fill(tofProjetcedHit, hit.eventHitWeight());
          const double dtof_hit = tof_hit - tofProjetcedHit;
          h_neutron_dtof_hit->fill(dtof_hit, hit.eventHitWeight());
          countSignal_dtof += signal::isSignal_dtof(dtof_hit);

          ///

          //const double tofExactHitPosition = generatorToExactHitPositionDistance/*[mm]*/ / Utils::neutron_angstrom_to_meters_per_second(lambda_true) /*[m/s]*/; //[ms]
          //h_neutron_tofExactHitPosition->fill(tofExactHitPosition, hit.eventHitWeight());
          //h_neutron_dtofExactHitPosition->fill(tof_conv - tofExactHitPosition, hit.eventHitWeight());




          //const double tofNoSmearPosition = sampleToWireNoSmearPositionDistance/*[mm]*/ / Utils::neutron_angstrom_to_meters_per_second(lambda_true) /*[m/s]*/; //[ms]
          //h_neutron_tofWireNoSmearPosition->fill(tofNoSmearPosition, hit.eventHitWeight());
          //h_neutron_dtofWireNoSmearPosition->fill(tof_conv - tofNoSmearPosition, hit.eventHitWeight());


          const double tofWithSmearPosition = sampleToWireWithSmearPositionDistance/*[mm]*/ / Utils::neutron_angstrom_to_meters_per_second(lambda_true) /*[m/s]*/; //[ms]
          //h_neutron_tofWireWithSmearPosition->fill(tofWithSmearPosition, hit.eventHitWeight());
          h_neutron_dtofWireWithSmearPosition->fill(tof_hit - tofWithSmearPosition, hit.eventHitWeight());



          h_neutron_Q_hit->fill(Q_hit, hit.eventHitWeight());
          const double dQ_hit = Q_hit - Q_true * Units::angstrom;
          h_neutron_dQ_hit->fill(dQ_hit, hit.eventHitWeight());
          countSignal_dQ += signal::isSignal_dQ(dQ_hit);
          //h_neutron_dQ_noSmearing_hit->fill(Q_hit_noSmearing - Q_true * Units::angstrom, hit.eventHitWeight());
          countSignal_dQ_noSmear += signal::isSignal_dQ(Q_hit_noSmearing - Q_true * Units::angstrom);

          //h_neutron_dthetaovertheta_hit->fill((theta_hit - theta_true) / theta_hit, hit.eventHitWeight());
          //h_neutron_dQoverQ_hit->fill((Q_hit-Q_true * Units::angstrom) / Q_hit, hit.eventHitWeight());
          //h_neutron_dlambdaoverlambda_hit->fill((lambda - Utils::neutronEKinToWavelength(ekin_true)/Units::angstrom) / lambdaHit, hit.eventHitWeight());

          h_neutron_lambdaFromTof_hit->fill(lambdaHit, hit.eventHitWeight());
          const double dlambda_hit = lambdaHit - lambda_true;
          h_neutron_dlambda_hit->fill(dlambda_hit, hit.eventHitWeight());
          countSignal_dlambda += signal::isSignal_dlambda(dlambda_hit);
          //h_neutron_dlambda_noSmearing_hit->fill(lambdaHitNoSmearing - lambda_true, hit.eventHitWeight());
          countSignal_dlambda_noSmear += signal::isSignal_dlambda(lambdaHitNoSmearing - lambda_true);

          //h_neutron_dthetavsdLambda_hit->fill(dlambda_hit, theta_hit-theta_true, hit.eventHitWeight());
          countSignal_dlambdaANDdth += signal::isSignal_dlambda(dlambda_hit) * signal::isSignal_dth(dth_hit); //if fulfills both criteria

          const double ekin_hit = Utils::neutron_angstrom_to_meV(lambdaHit); //Utils::neutron_angstrom_to_meV
          h_neutron_dekin_hit->fill(ekin_true/Units::meV - ekin_hit, hit.eventHitWeight());
          h_neutron_dekin_hit_center->fill(ekin_true/Units::meV - ekin_hit, hit.eventHitWeight());
          h_neutron_dekin_hit_center2->fill(ekin_true/Units::meV - ekin_hit, hit.eventHitWeight());

          const int tubeNumber = segL->volumeCopyNumber(3);
          const int strawNumber = segL->volumeCopyNumber(1);
          h_st_countingGas_hit->fill(strawNumber, tubeNumber, neutron->weight());


          //h_neutron_tofHit_minus_tofConv->fill(tof_hit - tof_conv, hit.eventHitWeight());

          h_st_global_hit->fill(strawNumberConv, tubeNumberConv, hit.eventHitWeight());
        }
      } // end of loop over converted neutrons/last segment in converter
      else if(segL->volumeName() == "TubeWall") { // check that the neutron ends in the Aluminium, this is a condition equivalent to a converted neutron
        const int tubeNumber = segL->volumeCopyNumber();

        h_st_tubeWall_hit->fill(tubeNumber, neutron->weight());
      }
      else if(segL->volumeName() == "StrawWall") { // check that the neutron ends in the Aluminium, this is a condition equivalent to a converted neutron
        const int tubeNumber = segL->volumeCopyNumber(2);
        const int strawNumber = segL->volumeCopyNumber();

        h_st_strawWall_hit->fill(strawNumber, tubeNumber, neutron->weight());
      }
    }//end of loop over primary neutrons

  } //end of event loop

  hc.saveToFile("bcs", true);

  return 0;
}
