#include "LOKI/DetectionFileCreator.hh"


DetectionFileCreator::DetectionFileCreator(const char* fileName, GriffDataRead::StrMap& userData): m_fileName(fileName) {
  this->detMcpl = mcpl_create_outfile(fileName);
  mcpl_hdr_set_srcname(this->detMcpl, "LOKI/DetectionFileCreator Class");
  mcpl_hdr_add_comment(this->detMcpl, "Particles in this file are actually detection events, and ekins are pixel ID's. Created with the DetectionFileCreator class.");
  mcpl_enable_universal_pdgcode(this->detMcpl, 99); //for file size reduction
  mcpl_enable_universal_weight(this->detMcpl, 1); //for file size reduction
  this->mcplParticle = mcpl_get_empty_particle(this->detMcpl);
  initiateMcplParticle();

  for(GriffDataRead::StrMap::iterator iter = userData.begin(); iter != userData.end(); ++iter) {
    std::string key = iter->first;
    std::string value = iter->second;
    mcpl_hdr_add_data(this->detMcpl, key.c_str(), value.size(), value.c_str());
  }
}

DetectionFileCreator::DetectionFileCreator(const char* fileName): m_fileName(fileName) {
  this->detMcpl = mcpl_create_outfile(fileName);
  mcpl_hdr_set_srcname(this->detMcpl, "LOKI DetectionFileCreator Class");
  mcpl_hdr_add_comment(this->detMcpl, "Particles in this file are actually detection events, and ekins are pixel ID's. Created with the DetectionFileCreator class.");
  mcpl_enable_universal_pdgcode(this->detMcpl, 99); //for file size reduction
  mcpl_enable_universal_weight(this->detMcpl, 1); //for file size reduction
  this->mcplParticle = mcpl_get_empty_particle(this->detMcpl);
  initiateMcplParticle();
}

DetectionFileCreator::~DetectionFileCreator() {
  mcpl_close_outfile(this->detMcpl);
  std::cout << "Created detection event file: " << this->m_fileName << std::endl;
}

void DetectionFileCreator::initiateMcplParticle() {
  this->mcplParticle->position[0] = 0;
  this->mcplParticle->position[1] = 0;
  this->mcplParticle->position[2] = 0;
  this->mcplParticle->direction[0] = 0;
  this->mcplParticle->direction[1] = 0;
  this->mcplParticle->direction[2] = 1;
}

void DetectionFileCreator::addDetectionEvent(const int pixelId, const double tof) {
  this->mcplParticle->ekin = pixelId;
  this->mcplParticle->time = tof;
  mcpl_add_particle(this->detMcpl, this->mcplParticle);
}
