#include "GriffDataRead/GriffDataReader.hh"
#include "GriffDataRead/DumpObj.hh"
#include "Utils/Format.hh"
#include "Utils/ArrayMath.hh"
#include "Units/Units.hh"
#include <cmath>
//Update March 2014: ROOT hists were not actually used. If needed, migrate to SimpleHists instead.
// #include "TH1D.h"
//#include "TFile.h"
// #include "TF1.h"

#include "GriffAnaUtils/SegmentIterator.hh"
#include "GriffAnaUtils/TrackFilter_PDGCode.hh"
#include "GriffAnaUtils/TrackFilter_Primary.hh"
#include "GriffAnaUtils/SegmentFilter_Volume.hh"
// #include "GriffAnaUtils/StepFilter_EKin.hh"
#include <cassert>

int main (int argc,char**argv) {

  ////////////////////////////////
  //Load data and extract setup:

  GriffDataReader dr(argc,argv);
  auto setup = dr.setup();
  auto & geo = setup->geo();
  auto & gen = setup->gen();
  if (gen.getName()!="G4StdGenerators/FlexGen"||geo.getName()!="G4SimPlaneScatter/GeoPlane") {
    printf("Error: Wrong setup for this analysis\n");
    return 1;
  }
  //todo: check that nothing is randomized (utility function in GriffAnaUtils
  //
  // if (gen.getParameterBoolean("randomize_energy"))
  //   return 1;

  const double eini = gen.getParameterDouble("fixed_energy_eV")*Units::eV;
  const double theta_ini = gen.getParameterDouble("fixed_polarangle_deg")*Units::deg;
  const double thickness = geo.getParameterDouble("thickness_mm")*Units::mm;
  //const double transverse_extent = geo.getParameterDouble("transverse_extent_m")*Units::m;
  const std::string material_name = geo.getParameterString("material");

  std::string title;
  Utils::string_format(title,"%g meV %s with theta=%g degree passing through %g cm of %s [%s]",
                       eini/(0.001*Units::eV),
                       gen.getParameterString("particleName").c_str(),//fixme: make sure it is not pdgCode based
                       theta_ini/Units::deg,
                       thickness/Units::cm,
                       geo.getParameterString("material").c_str(),
                       setup->metaData("G4PhysicsList").c_str());
  setup->dump();
  printf("TITLE: %s\n",title.c_str());

  //Filters:
  GriffAnaUtils::SegmentIterator valid_segments(&dr);
  valid_segments.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));
  valid_segments.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  valid_segments.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Plane"));

  //hists:
  // TH1D h_forward_deltatheta("forward_deltatheta",("forward_deltatheta [degrees] "+title).c_str(),1000,-190.0,190.0);

  //file
  std::string filename = "neutronsim.txt";
  Utils::string_format(filename,"neutronsim_slab_%gcm_material_%s.txt",
                       thickness/Units::cm,
                       material_name.c_str());


  FILE * file = fopen(filename.c_str(),"w");


  //analysis loop:
  unsigned neutrondiesinside(0);
  unsigned nbadexitvol(0);
  unsigned nneutrons(0);
  unsigned nforward(0);
  unsigned nUnscattered(0);
  //unsigned nforward_unscattered(0);
  unsigned nback(0);
  bool first=true;

  double common_inidir[3] = {0,0,0};
  double common_inipos[3] = {0,0,0};

  double unscat_exitpos[3] = {0,0,0};
  double unscat_exitdir[3] = {0,0,0};

  uint64_t nScattersTot(0);
  uint64_t nScattersEvts(0);

  while (dr.loopEvents()) {
    // if (nneutrons>100)
    //   break;
    ++nneutrons;
    unsigned nvalidsegments(0);
    while (auto segment = valid_segments.next()){
      ++nvalidsegments;
      if (nvalidsegments!=1) {
        printf("Error: More than one valid segments seen in event!\n");
        return 1;
      }
      if (segment->iSegment()!=0) {
        printf("Error: Segment is not first segment!\n");
        return 1;
      }
      auto nextSeg = segment->getNextSegment();
      if (!nextSeg) {
        //printf("neutron never makes it out of the plane!\n");
        ++neutrondiesinside;
        continue;
      }
      bool forward(nextSeg->volumeName()=="ForwardScatterPlane");
      if (!forward&&nextSeg->volumeName()!="BackScatterPlane") {
        if (nextSeg->volumeName()!="World")
          printf("WARNING: Unexpected exit volume: %s\n",nextSeg->volumeName().c_str());
        ++nbadexitvol;
        continue;
      }
      if (forward)
        ++nforward;
      else
        ++nback;
      //GriffDataRead::dump(segment);printf("\n");

      //const double unitz[3] = { 0, 0, 1 };
      //const double minusunitz[3] = { 0, 0, -1 };

      //if only Transportation and hadElastic present, do:
      //NSCATTERS INITIALX Y Z DIRX Y Z EXIT X Y Z DIR X Y Z
      assert(segment->lastStep()->postProcessDefinedStep()=="Transportation");

      //count interaction types:
      int nScatters(0);
      int nTransports(0);
      int nOthers(0);
      for (auto step = segment->stepBegin();step!=segment->stepEnd();++step) {
        if (step->stepLength()==0.0)
          continue;
        //        printf("steplength %g\n",step->stepLength());
        //printf("%s %s\n",step->preProcessDefinedStep().c_str(),step->postProcessDefinedStep().c_str());
        if (step->postProcessDefinedStep()=="hadElastic") ++nScatters;
        else if (step->postProcessDefinedStep()=="Transportation") ++nTransports;
        else ++nOthers;
      }
      if (nOthers>0) {
        printf("WARNING: Ignoring event with exiting neutron but non-hadElastic interaction\n");
        --nneutrons;
        continue;
      }
      assert(nTransports==1);

      auto inipos = segment->firstStep()->preGlobalArray();
      auto finpos = nextSeg->firstStep()->preGlobalArray();
      double inidir[3];
      double findir[3];
      Utils::normalise(segment->firstStep()->preMomentumArray(),inidir);
      Utils::normalise(nextSeg->firstStep()->preMomentumArray(),findir);
      assert(inipos[0]==0&&inipos[1]==0&&inipos[2]==0);

      if (first) {
        common_inipos[0] = inipos[0];
        common_inipos[1] = inipos[1];
        common_inipos[2] = inipos[2];
        common_inidir[0] = inidir[0];
        common_inidir[1] = inidir[1];
        common_inidir[2] = inidir[2];
        first=false;
        printf("\n");
        //FIXME: Put some meta-data here!!! (like inipos, inidir, plane thickness, ...
        fprintf(file,"#\n");
        fprintf(file,"# Simulation of neutrons scattering through a single slab\n");
        fprintf(file,"# Created by Thomas Kittelmann (ESS) with Geant4 in the dgcode framework\n");
        fprintf(file,"#\n");
        fprintf(file,"# All positions are given in units of centimetres.\n");
        //neutron energy, scatter type (material), ...
        fprintf(file,"#\n");
        fprintf(file,"# The slab extends from z=0 to z=%f and has very large transverse dimensions\n",thickness/Units::cm);
        fprintf(file,"# Slab material: %s\n",material_name.c_str());
        fprintf(file,"# Neutron energy: %g meV\n",eini/(0.001*Units::eV));
        fprintf(file,"#\n");
        fprintf(file,"# All neutrons in this file starts at:\n");
        fprintf(file,"#                (%14.8f,%14.8f,%14.8f)\n",inipos[0]/Units::cm,inipos[1]/Units::cm,inipos[2]/Units::cm);
        fprintf(file,"# All neutrons in this file starts with direction:\n");
        fprintf(file,"#                (%14.8f,%14.8f,%14.8f)\n",inidir[0],inidir[1],inidir[2]);
        fprintf(file,"#\n");
        fprintf(file,"# The columns mean:\n");
        fprintf(file,"#    BACKS: 0 if the neutron exited through the back face, 1 if through the forward face\n");
        fprintf(file,"#    NSCAT: Number of scatterings taking place before the neutrons exit the slab\n");
        fprintf(file,"#    EXIT_POS: Position where the neutron exits the slab\n");
        fprintf(file,"#    EXIT_DIR: Direction of the neutron when it exits the slab\n");
        fprintf(file,"#    WEIGHT: How many simulated neutrons undertook this history\n");
        fprintf(file,"#\n");

        //Units are cm.
        fprintf(file,"#  WEIGHT BACK NSCAT     EXIT_POS_X     EXIT_POS_Y     EXIT_POS_Z     EXIT_DIR_X     EXIT_DIR_Y     EXIT_DIR_Z\n");
      };//endif first
      if (!(common_inipos[0]==inipos[0]&&
            common_inipos[1]==inipos[1]&&
            common_inipos[2]==inipos[2]&&
            common_inidir[0]==inidir[0]&&
            common_inidir[1]==inidir[1]&&
            common_inidir[2]==inidir[2] )) {
        printf("ERROR: Initial position or direction not fixed!\n");
        return 1;
      }
      if (nScatters==0&&forward&&nTransports==1) {
        ++nUnscattered;
        if (nUnscattered==1) {
          unscat_exitpos[0]=finpos[0];
          unscat_exitpos[1]=finpos[1];
          unscat_exitpos[2]=finpos[2];
          unscat_exitdir[0]=findir[0];
          unscat_exitdir[1]=findir[1];
          unscat_exitdir[2]=findir[2];
        }
        continue;//write all at the end
      }

      nScattersTot += nScatters;
      ++nScattersEvts;
      fprintf(file," %8i %4i %5i %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n",
              1,
              (forward?0:1),
              nScatters,
              finpos[0]/Units::cm,finpos[1]/Units::cm,finpos[2]/Units::cm,
              findir[0],findir[1],findir[2]
              );
      // if (forward) {
      //   //const double exit_costheta = Utils::costheta(nextSeg->firstStep()->preMomentumArray(),unitz);
      //   const double delta_costheta = Utils::costheta(nextSeg->firstStep()->preMomentumArray(),
      //                                                 segment->firstStep()->preMomentumArray());

      //   if (delta_costheta>0.999999)
      //     ++nforward_unscattered;
      //   //printf("%g\n",delta_costheta);
      //   h_forward_deltatheta.Fill(acos(delta_costheta));
      // }

    }
  }

  fprintf(file," %8i %4i %5i %14.8f %14.8f %14.8f %14.8f %14.8f %14.8f\n",
          nUnscattered,
          (true?0:1),
          0,
          unscat_exitpos[0]/Units::cm,unscat_exitpos[1]/Units::cm,unscat_exitpos[2]/Units::cm,
          unscat_exitdir[0],unscat_exitdir[1],unscat_exitdir[2]
          );

  printf("Wrote output to file %s\n",filename.c_str());
  fclose(file);

  printf("%g %% of neutrons ended up inside the plane\n",neutrondiesinside*100.0/nneutrons);
  printf("%g %% of neutrons exited through forward direction\n",nforward*100.0/nneutrons);
  printf("%g %% of neutrons exited through backward direction\n",nback*100.0/nneutrons);
  printf("%g %% of neutrons exited through side or other bad direction\n",nbadexitvol*100.0/nneutrons);

  printf("%g number of scatterings in average exiting event\n",double(nScattersTot)/(nUnscattered+nScattersEvts));

  // TFile f("hist_plane.root","RECREATE");
  // h_forward_deltatheta.Write();
  // f.Close();

  return 0;

}
