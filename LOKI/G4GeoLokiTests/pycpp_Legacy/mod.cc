// Test-only module: a FROZEN COPY of the original (pre-BankTransform, commit
// a653830) nominal implementation of the bank placement, the AimHelper pixel
// centres and the PixelatedBanks (x, y) pixel id. The transition tests compare
// the current implementation (BcsBanks::getBankTransform, used by the Geant4
// geometry, AimHelper and the 3D PixelatedBanks::getPixelId) with this copy:
// as long as the geometry parameters are nominal, both must agree exactly.
//
// Do NOT change this file when changing the production code: it is the
// reference. It only reads the (unchanged) nominal parameter tables through the
// public getters of BcsBanks / BcsPack / BcsTube / PixelatedBanks. (Only adapted to API changes:
// the pixels per straw are read from its own PixelatedBanks object, as they are per object; the old
// tube numbering and the Larmor 2022 setup were removed from LOKI (they are in LOKI/Larmor): only the
// LOKI geometry with the current tube numbering is compared.)
//
// Nothing in here is used by production code.

#include "Core/Python.hh"
#include "G4GeoLoki/AimHelper.hh"
#include "G4GeoLoki/PixelatedBanks.hh"
#include "G4GeoLoki/BcsBanks.hh"
#include "G4GeoLoki/BcsPack.hh"
#include "G4GeoLoki/BcsTube.hh"
#include "G4RotationMatrix.hh"
#include <cmath>
#include <array>
#include <algorithm>
#include <stdexcept>

namespace {

  // ---------- verbatim logic of the original AimHelper (a653830) ----------
  void coordinateRotation(double &x, double &y, const double angle) {
    double tempX = std::cos(angle) * x - std::sin(angle) * y;
    double tempY = std::sin(angle) * x + std::cos(angle) * y;
    x = tempX;
    y = tempY;
  }

  int getPackId(const int bankId, const int tubeId, const bool isOldPixelNumbering) {
    const int numberOfPacks = BcsBanks::getNumberOfPacksByBankId(bankId);
    const int normalPackId = isOldPixelNumbering ?
                             (int) tubeId / 8 :
                             (int) (tubeId % (numberOfPacks * 2)) / 2;
    return !BcsBanks::areTubesInverselyNumbered(bankId) ? normalPackId : ((numberOfPacks - 1) - normalPackId);
  }

  int getInPackTubeId(const int bankId, const int tubeId, const bool isOldPixelNumbering) {
    const int numberOfPacks = BcsBanks::getNumberOfPacksByBankId(bankId);
    const int newTubeIdConvertedToOldId = ((tubeId % 2) * 4) + ((int) tubeId / (numberOfPacks * 2));
    if(isOldPixelNumbering){
      return BcsBanks::areTubesInverselyNumbered(bankId) ? (tubeId + 4) % 8 : tubeId % 8;
    }
    else{
      return BcsBanks::areTubesInverselyNumbered(bankId) ? (newTubeIdConvertedToOldId + 4) % 8 : newTubeIdConvertedToOldId % 8;
    }
  }

  struct Legacy {
    PixelatedBanks pb;
    int nbanks;
    Legacy(double rear, int n, int nb) : pb(rear, n, nb), nbanks(nb) {}

    int bankId(int pixelId) const {
      for (int b = 0; b < pb.getNumberOfBanks(); b++)
        if (pixelId < pb.getBankPixelOffset(b+1)) return b;
      throw std::runtime_error("pixel id out of range");
    }
    int tubeId(int pixelId, int b) const {
      return (pixelId - pb.getBankPixelOffset(b)) / (pb.getNumberOfPixelsInStraw(b) * 7);
    }
    int strawId(int pixelId, int b, int t) const {
      const int inTube = pixelId - pb.getBankPixelOffset(b) - t * 7 * pb.getNumberOfPixelsInStraw(b);
      return inTube / pb.getNumberOfPixelsInStraw(b);
    }
    double pixelPositionInStraw(int pixelId, int b) const {
      const int n = pb.getNumberOfPixelsInStraw(b);
      const int loc = pixelId % n;
      const double L = BcsBanks::getStrawLengthByBankId(b);
      const double position = -0.5 * L + (loc + 0.5) * L / n;
      return BcsBanks::areTubesInverselyNumbered(b) ? -position : position;
    }

    // original AimHelper::getPixelCentreCoordinates
    std::array<double,3> pixelCentre(int pixelId) const {
      const bool old = false;
      const int b = bankId(pixelId);
      const int t = tubeId(pixelId, b);
      const int inPackTubeId = getInPackTubeId(b, t, old);
      const int packId = getPackId(b, t, old);
      const int s = strawId(pixelId, b, t);
      double positionZ = pixelPositionInStraw(pixelId, b);
      double positionX = BcsTube::getStrawPositionX(s);
      double positionY = BcsTube::getStrawPositionY(s);
      coordinateRotation(positionX, positionY, BcsPack::getTubeRotation());
      positionX += BcsPack::getHorizontalTubeOffset(inPackTubeId);
      positionY += BcsPack::getVerticalTubeOffset(inPackTubeId);
      coordinateRotation(positionX, positionY, BcsBanks::getPackRotation());
      positionX += BcsBanks::getPackPositionInBank(b, packId, 2);
      positionY += BcsBanks::getPackPositionInBank(b, packId, 1);
      positionZ += BcsBanks::getPackPositionInBank(b, packId, 0);
      coordinateRotation(positionY, positionX, BcsBanks::getBankRotation(b, 2));
      coordinateRotation(positionZ, positionY, BcsBanks::getBankRotation(b, 0));
      coordinateRotation(positionZ, positionX, -BcsBanks::getBankRotation(b, 1));
      positionX += pb.getBankPosition(b, 0);
      positionY += pb.getBankPosition(b, 1);
      positionZ += pb.getBankPosition(b, 2);
      return {positionX, positionY, positionZ};
    }

    // original PixelatedBanks::getPositionPixelId / getPixelId (global x or y)
    int pixelId2D(int b, int t, int s, double x, double y) const {
      const int n = pb.getNumberOfPixelsInStraw(b);
      const double L = BcsBanks::getStrawLengthByBankId(b);
      const double pixelLength = L / n;
      int inStraw;
      if (BcsBanks::isVertical(b)) {
        const double strawBegin = pb.getBankPosition(b, 1) - 0.5 * L;
        inStraw = std::floor((y - strawBegin) / pixelLength);
      } else {
        const double strawBegin = pb.getBankPosition(b, 0) - 0.5 * L;
        inStraw = (n - 1) - (int) std::floor((x - strawBegin) / pixelLength);
      }
      return pb.getBankPixelOffset(b) + (t * 7 + s) * n + inStraw;
    }

    // original Geant4 bank placement (GeoBCSBanks): frame rotation rotateY(a1) rotateX(a0) rotateZ(a2), as rows
    static py::tuple frameRotation(int b) {
      G4RotationMatrix r;
      r.rotateY(BcsBanks::getBankRotation(b, 1));
      r.rotateX(BcsBanks::getBankRotation(b, 0));
      r.rotateZ(BcsBanks::getBankRotation(b, 2));
      return py::make_tuple(py::make_tuple(r.xx(), r.xy(), r.xz()),
                            py::make_tuple(r.yx(), r.yy(), r.yz()),
                            py::make_tuple(r.zx(), r.zy(), r.zz()));
    }
    py::tuple bankPosition(int b) const {
      return py::make_tuple(pb.getBankPosition(b, 0), pb.getBankPosition(b, 1), pb.getBankPosition(b, 2));
    }
    py::tuple pyPixelCentre(int pixelId) const {
      auto c = pixelCentre(pixelId);
      return py::make_tuple(c[0], c[1], c[2]);
    }
    int pyPixelId2D(int b, int t, int s, double x, double y) const { return pixelId2D(b, t, s, x, y); }
  };

  // Compare ALL pixels of the geometry: the current AimHelper centres vs this copy, and the
  // pixel ids of the current 3D PixelatedBanks::getPixelId and of this copy's (x, y) one vs
  // the expected id, at the centre, +-0.49 / +-0.51 pixel pitch along the straw (within the
  // straw) and 3 mm across it (four directions). Returns
  // (pixels, largest centre deviation [mm], pixel of it, pixel id checks, pixel id failures,
  //  first failing pixel or -1).
  py::tuple compare_all(double rear, int n, int nbanks) {
    Legacy legacy(rear, n, nbanks);
    AimHelper aim = (nbanks == 9) ? AimHelper(rear, n) : AimHelper(rear, n, nbanks);
    const int total = aim.getTotalNumberOfPixels();
    double worst = 0.; int worstPixel = -1;
    long nchecks = 0, nfail = 0; int firstFail = -1;
    for (int p = 0; p < total; p++) {
      auto cur = aim.getPixelCentreCoordinates(p);
      auto ref = legacy.pixelCentre(p);
      const double d = std::max({std::fabs(std::get<0>(cur) - ref[0]), std::fabs(std::get<1>(cur) - ref[1]), std::fabs(std::get<2>(cur) - ref[2])});
      if (!(d <= worst)) { worst = d; worstPixel = p; }
      const int b = legacy.bankId(p), t = legacy.tubeId(p, b), s = legacy.strawId(p, b, t);
      const int j = p % n;
      // straw axis (direction of increasing pixel index) and pitch, from the reference centres
      const double L = BcsBanks::getStrawLengthByBankId(b), pitch = L / n;
      std::array<double,3> u;
      if (n > 1) {
        auto a = legacy.pixelCentre(j < n - 1 ? p : p - 1);
        auto c = legacy.pixelCentre(j < n - 1 ? p + 1 : p);
        for (int i = 0; i < 3; i++) u[i] = (c[i] - a[i]) / pitch;
      } else {
        // one pixel: the straw axis is the world x (horizontal) or y (vertical) axis
        u = {0., 0., 0.};
        u[BcsBanks::isVertical(b) ? 1 : 0] = BcsBanks::isVertical(b) ? 1. : -1.;
      }
      // two directions perpendicular to u
      std::array<double,3> e = {0., 0., 1.};
      std::array<double,3> v = {u[1]*e[2]-u[2]*e[1], u[2]*e[0]-u[0]*e[2], u[0]*e[1]-u[1]*e[0]};
      double nv = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); for (auto& x : v) x /= nv;
      std::array<double,3> w = {u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]};
      struct Probe { std::array<double,3> pos; int expected; };
      std::vector<Probe> probes;
      probes.push_back({ref, p});
      for (double f : {0.49, -0.49, 0.51, -0.51}) {
        const int dj = (f > 0.5) ? 1 : (f < -0.5 ? -1 : 0);
        if (j + dj < 0 || j + dj >= n) continue;  // beyond the straw ends: not defined
        probes.push_back({{ref[0] + f*pitch*u[0], ref[1] + f*pitch*u[1], ref[2] + f*pitch*u[2]}, p + dj});
      }
      for (double sgn : {3.0, -3.0}) {
        probes.push_back({{ref[0] + sgn*v[0], ref[1] + sgn*v[1], ref[2] + sgn*v[2]}, p});
        probes.push_back({{ref[0] + sgn*w[0], ref[1] + sgn*w[1], ref[2] + sgn*w[2]}, p});
      }
      for (const auto& pr : probes) {
        nchecks++;
        const int id3 = legacy.pb.getPixelId(b, t, s, pr.pos[0], pr.pos[1], pr.pos[2]);
        const int id2 = legacy.pixelId2D(b, t, s, pr.pos[0], pr.pos[1]);
        if (id3 != pr.expected || id2 != pr.expected) {
          nfail++;
          if (firstFail < 0) firstFail = p;
        }
      }
    }
    return py::make_tuple(total, worst, worstPixel, nchecks, nfail, firstFail);
  }
}

PYTHON_MODULE( mod )
{
  py::class_<Legacy>(mod, "Legacy")
    .def(py::init<double,int,int>())
    .def("pixelCentre", &Legacy::pyPixelCentre)
    .def("pixelId2D", &Legacy::pyPixelId2D)
    .def("bankPosition", &Legacy::bankPosition)
    .def_static("frameRotation", &Legacy::frameRotation)
    ;
  mod.def("compare_all", &compare_all);
}
