#ifndef GriffB10Common_DetectionFileCreator_hh
#define GriffB10Common_DetectionFileCreator_hh

#include "MCPL/mcpl.h"
#include <iostream>

class DetectionFileCreator {
public:
  DetectionFileCreator(const char* fileName);
  ~DetectionFileCreator();

  void addDetectionEvent(const int pixelId, const double tof);

private:
  std::string m_fileName;
  mcpl_outfile_t detMcpl;
  mcpl_particle_t *mcplParticle;

  void initiateMcplParticle();
};

#endif
