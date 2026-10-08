// Test-only Griff analysis for test_geantino_aiming: for every event, prints the
// straws the primary geantino crossed (their CountingGas volume, which contains the
// pixel centres; the Converter is only the thin shell around it), with the range of
// pixel ids between its entry and exit points, and whether it had crossed an absorber
// (boron masks, B4C/Al panels, beamstop) before. One line per event:
//   EVENT <index> <bank>/<tube>/<straw>:<first pixel>-<last pixel>[@<absorber>] ...
// ('@<absorber>' = the name of the first absorber volume crossed before this straw). The pixel ids come from the
// copy numbers and the 3D PixelatedBanks::getPixelId, as in the analysis programs
// (local index clamped to the straw for the entry/exit points on the straw ends).

#include "GriffAnaUtils/All.hh"
#include "G4GeoLoki/PixelatedBanks.hh"
#include "Units/Units.hh"
#include <iostream>
#include <string>
#include <algorithm>

int main(int argc, char **argv) {
  GriffDataReader dr(argc, argv);
  auto setup = dr.setup();
  auto &geo = setup->geo();
  auto userData = setup->userData();
  const double rear = geo.getParameterDouble("rear_detector_distance_m") * Units::m;
  const int n = std::stoi(userData["analysis_straw_pixel_number"]);
  const std::string bankCalibration = geo.hasParameterString("bank_calibration") ? geo.getParameterString("bank_calibration") : BankCalibration::nominalName;
  PixelatedBanks banks(rear, n, 9, bankCalibration);

  GriffAnaUtils::TrackIterator geantinos(&dr);
  geantinos.addFilter(new GriffAnaUtils::TrackFilter_Primary());
  geantinos.addFilter(new GriffAnaUtils::TrackFilter_PDGCode(999));

  // diagnostic: Converter segments whose exit point (the point the masking analysis uses)
  // gives a pixel index outside the straw (-1 or N: an id of the neighbouring straw)
  long nConverterSegments = 0, nExitOutsideStraw = 0;

  while (dr.loopEvents()) {
    std::cout << "EVENT " << dr.eventIndexInCurrentFile();
    while (auto trk = geantinos.next()) {
      std::string absorber;
      for (auto seg = trk->segmentBegin(); seg != trk->segmentEnd(); ++seg) {
        const std::string vol = seg->volumeName();
        if (vol == "Converter") {
          const int straw = seg->volumeCopyNumber(1), tube = seg->volumeCopyNumber(3), bank = seg->volumeCopyNumber(5);
          const int base = banks.getBankPixelOffset(bank) + (tube * 7 + straw) * n;
          auto last = seg->lastStep();
          const int id = banks.getPixelId(bank, tube, straw, last->postGlobalX(), last->postGlobalY(), last->postGlobalZ());
          nConverterSegments++;
          if (id < base || id >= base + n) nExitOutsideStraw++;
        }
        if (vol.find("BoronMask-") != std::string::npos || vol == "B4CPanel" || vol == "AlPanel" || vol.find("Beamstop") != std::string::npos) {
          if (absorber.empty()) absorber = vol;
        } else if (vol == "CountingGas") {
          // copy numbers one level deeper than from the Converter used by the analysis programs
          const int straw = seg->volumeCopyNumber(2), tube = seg->volumeCopyNumber(4), bank = seg->volumeCopyNumber(6);
          const int base = banks.getBankPixelOffset(bank) + (tube * 7 + straw) * n;
          auto clampToStraw = [&](int id) { return std::min(std::max(id, base), base + n - 1); };
          auto first = seg->firstStep(), last = seg->lastStep();
          const int a = clampToStraw(banks.getPixelId(bank, tube, straw, first->preGlobalX(), first->preGlobalY(), first->preGlobalZ()));
          const int b = clampToStraw(banks.getPixelId(bank, tube, straw, last->postGlobalX(), last->postGlobalY(), last->postGlobalZ()));
          std::cout << ' ' << bank << '/' << tube << '/' << straw << ':' << std::min(a, b) << '-' << std::max(a, b) << (absorber.empty() ? "" : "@" + absorber);
        }
      }
    }
    std::cout << '\n';
  }
  std::cout << "CONVERTER_EXITS " << nConverterSegments << " OUTSIDE_STRAW " << nExitOutsideStraw << '\n';
  return 0;
}
