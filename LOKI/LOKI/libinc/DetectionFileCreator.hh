#ifndef LOKI_DetectionFileCreator_hh
#define LOKI_DetectionFileCreator_hh

#include "mcpl.h"
#include "GriffDataRead/Setup.hh"
#include <iostream>

class DetectionFileCreator {
public:
  DetectionFileCreator(const char* fileName, GriffDataRead::StrMap& userData);
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
