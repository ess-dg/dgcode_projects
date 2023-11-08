#include "GriffDataRead/GriffDataReader.hh"
#include "GriffAnaUtils/SegmentIterator.hh"
#include "GriffAnaUtils/TrackFilter_PDGCode.hh"
#include "GriffAnaUtils/TrackFilter_Primary.hh"
#include "GriffAnaUtils/SegmentFilter_Volume.hh"

#include "SimpleHists/HistCollection.hh"
#include "Utils/Format.hh"
#include "Utils/ArrayMath.hh"
#include "Units/Units.hh"

int main (int argc,char**argv) {

  ////////////////////////////////
  //Open data file(s):

  GriffDataReader dr(argc,argv);

  ////////////////////////////////
  //Verify that we run on the right kind of files:

  auto setup = dr.setup();
  auto & geo = setup->geo();
  auto & gen = setup->gen();
  if (gen.getName()!="G4StdGenerators/FlexGen"||geo.getName()!="G4SimPlaneScatter/GeoPlane") {
    printf("Error: Wrong setup for this analysis\n");
    return 1;
  }

  /////////////////////////////////
  //Book a "collection" of histograms (actually just 1) that we will fill while
  //looping over the file:

  SimpleHists::HistCollection hc;

  std::string title;
  Utils::string_format(title,"Scatterings in %g cm %s plates",geo.getParameterString("material").c_str());
  auto h_scat = hc.book1D(title,100,0.0,180,"scat");
  h_scat->setXLabel("degree");

  /////////////////////////////////
  //Data filters (we want primary neutrons in the volume named "Plane"):
  GriffAnaUtils::SegmentIterator valid_segments(&dr);
  valid_segments.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(2112));
  valid_segments.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  valid_segments.addFilter(new GriffAnaUtils::SegmentFilter_Volume("Plane"));

  while (dr.loopEvents()) {
    while (auto segment = valid_segments.next()){
      double scat_angle = Utils::theta(segment->firstStep()->preMomentumArray(),segment->lastStep()->postMomentumArray());
      if (scat_angle!=scat_angle)
        continue;//protect against NaN (in case of null vectors)

      //Only fill when actual scattering occured, otherwise the first bin will dominate wildly:
      if (fabs(scat_angle)>0.00001*Units::deg)
        h_scat->fill(scat_angle/Units::deg);
    }
  }

  printf("Saving histograms\n");
  hc.saveToFile("planescatter.shist",true);
  return 0;

}
