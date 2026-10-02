// Test-only helper module for the G4GeoLokiTests package.
//
// Exposes (read-only) internals of the G4GeoLoki library which are not
// available through G4GeoLoki.LokiAimHelper, plus a minimal Geant4 navigation
// utility, so that the test scripts can pin the current behaviour of:
//
//   * the nominal bank/pack/tube/straw placement parameters (BcsBanks, BcsPack, BcsTube),
//   * PixelatedBanks::getPixelId (hit position + copy numbers -> pixel id), which
//     is the function used by all analysis programs,
//   * the constructed Geant4 geometry (GeoBCSBanks), by locating points with a
//     G4Navigator and reporting the touchable history (volume names, copy
//     numbers) and the local coordinates.
//
// Nothing in here is used by production code.

#include "Core/Python.hh"
#include <pybind11/stl.h>
#include "G4Interfaces/GeoConstructBase.hh"
#include "G4GeoLoki/PixelatedBanks.hh"
#include "G4GeoLoki/BcsBanks.hh"
#include "G4GeoLoki/BcsPack.hh"
#include "G4GeoLoki/BcsTube.hh"
#include "LokiMasking/MaskFileCreator.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TouchableHistory.hh"
#include "G4AffineTransform.hh"
#include "G4ThreeVector.hh"
#include <memory>
#include <iostream>
#include <cstdio>

namespace {

  G4VPhysicalVolume* s_world = nullptr;
  std::unique_ptr<G4Navigator> s_nav;

  //Calls Construct() on a geometry module instance (e.g. G4GeoLoki.GeoBCSBanks)
  //and keeps the returned world volume for subsequent locate(..) calls. Note
  //that the GeoConstructBase::place(..) calls perform the standard Geant4
  //overlap check for each placement during construction (printed to G4cout).
  int construct_world( py::object o )
  {
    if ( !py::isinstance<G4Interfaces::GeoConstructBase>( o ) )
      throw std::runtime_error("construct_world: object is not a geometry module");
    auto geo = py::cast<G4Interfaces::GeoConstructBase*>( o );
    s_world = geo->Construct();
    if (!s_world)
      throw std::runtime_error("construct_world: Construct() returned null");
    s_nav = std::make_unique<G4Navigator>();
    s_nav->SetWorldVolume(s_world);
    std::cout.flush();
    std::cerr.flush();
    std::fflush(stdout);
    std::fflush(stderr);
    return (int)s_world->GetLogicalVolume()->GetNoDaughters();
  }

  //Locates the global point (x,y,z) [mm] in the geometry previously built with
  //construct_world. Returns a tuple:
  //   ( [(physvolname, copynumber, logvolname), ...] from deepest volume up to the world,
  //     (local x, local y, local z) of the point in the deepest volume,
  //     (dx,dy,dz) global direction of the local z axis of the deepest volume,
  //     (ox,oy,oz) global position of the origin of the deepest volume )
  py::tuple locate( double x, double y, double z )
  {
    if (!s_nav)
      throw std::runtime_error("locate: call construct_world first");
    G4ThreeVector p(x,y,z);
    s_nav->LocateGlobalPointAndSetup(p, nullptr, false, true);
    std::unique_ptr<G4TouchableHistory> th(s_nav->CreateTouchableHistory());
    py::list hist;
    for (int d = 0; d <= th->GetHistoryDepth(); ++d) {
      auto pv = th->GetVolume(d);
      hist.append(py::make_tuple(std::string(pv->GetName()), th->GetCopyNumber(d),
                                 std::string(pv->GetLogicalVolume()->GetName())));
    }
    G4AffineTransform g2l = s_nav->GetGlobalToLocalTransform();
    G4AffineTransform l2g = s_nav->GetLocalToGlobalTransform();
    G4ThreeVector loc = g2l.TransformPoint(p);
    G4ThreeVector ax = l2g.TransformAxis(G4ThreeVector(0,0,1));
    G4ThreeVector org = l2g.TransformPoint(G4ThreeVector(0,0,0));
    return py::make_tuple( hist,
                           py::make_tuple(loc.x(),loc.y(),loc.z()),
                           py::make_tuple(ax.x(),ax.y(),ax.z()),
                           py::make_tuple(org.x(),org.y(),org.z()) );
  }

  //The volumes placed directly in the world of the geometry built with construct_world:
  //[(physvolname, copynumber, logvolname, frame rotation as 3 rows, translation (x,y,z) [mm]), ...]
  //(the frame rotation is what G4PVPlacement stores: the inverse of the object rotation).
  py::list world_daughters()
  {
    if (!s_world)
      throw std::runtime_error("world_daughters: call construct_world first");
    py::list out;
    auto lv = s_world->GetLogicalVolume();
    for (size_t i = 0; i < lv->GetNoDaughters(); ++i) {
      auto pv = lv->GetDaughter(i);
      const G4RotationMatrix* r = pv->GetRotation();
      G4RotationMatrix identity;
      if (!r) r = &identity;
      const G4ThreeVector t = pv->GetTranslation();
      out.append(py::make_tuple(std::string(pv->GetName()), pv->GetCopyNo(), std::string(pv->GetLogicalVolume()->GetName()),
                                py::make_tuple(py::make_tuple(r->xx(), r->xy(), r->xz()),
                                               py::make_tuple(r->yx(), r->yy(), r->yz()),
                                               py::make_tuple(r->zx(), r->zy(), r->zz())),
                                py::make_tuple(t.x(), t.y(), t.z())));
    }
    return out;
  }

  //Every physical volume of the geometry built with construct_world, depth first:
  //[(depth, physvolname, copynumber, logvolname, global translation (x,y,z) [mm],
  //  global object rotation as 3 rows), ...] (depth 0 = the world).
  void addTree(py::list& out, const G4VPhysicalVolume* pv, int depth, const G4RotationMatrix& Rp, const G4ThreeVector& tp)
  {
    const G4RotationMatrix R = Rp * pv->GetObjectRotationValue();
    const G4ThreeVector t = tp + Rp * pv->GetObjectTranslation();
    out.append(py::make_tuple(depth, std::string(pv->GetName()), pv->GetCopyNo(), std::string(pv->GetLogicalVolume()->GetName()),
                              py::make_tuple(t.x(), t.y(), t.z()),
                              py::make_tuple(py::make_tuple(R.xx(), R.xy(), R.xz()),
                                             py::make_tuple(R.yx(), R.yy(), R.yz()),
                                             py::make_tuple(R.zx(), R.zy(), R.zz()))));
    auto lv = pv->GetLogicalVolume();
    for (size_t i = 0; i < lv->GetNoDaughters(); ++i)
      addTree(out, lv->GetDaughter(i), depth + 1, R, t);
  }
  py::list world_tree()
  {
    if (!s_world)
      throw std::runtime_error("world_tree: call construct_world first");
    py::list out;
    addTree(out, s_world, 0, G4RotationMatrix(), G4ThreeVector());
    return out;
  }

  //Thin wrapper around PixelatedBanks (the class used by the analysis programs
  //to convert hits to pixel ids).
  struct PixelCalc {
    PixelatedBanks pb;
    PixelCalc(double rear, int n, int nbanks) : pb(rear, n, nbanks) {}
    int getPixelId(int bank, int tube, int straw, double x, double y) const { return pb.getPixelId(bank,tube,straw,x,y); }
    int getPixelId3D(int bank, int tube, int straw, double x, double y, double z) const { return pb.getPixelId(bank,tube,straw,x,y,z); }
    //bank transform: (rotation as 3 rows, translation)
    py::tuple getBankTransform(int bank, bool larmor2022) const {
      const BankTransform t = pb.getBankTransform(bank, larmor2022);
      const auto& R = t.rotation;
      return py::make_tuple( py::make_tuple( py::make_tuple(R[0][0],R[0][1],R[0][2]),
                                             py::make_tuple(R[1][0],R[1][1],R[1][2]),
                                             py::make_tuple(R[2][0],R[2][1],R[2][2]) ),
                             py::make_tuple(t.translation[0],t.translation[1],t.translation[2]) );
    }
    double getBankPosition(int bank, int axis) const { return pb.getBankPosition(bank,axis); }
    int getTotalNumberOfPixels() { return pb.getTotalNumberOfPixels(); }
    int getNumberOfBanks() const { return pb.getNumberOfBanks(); }
  };
}

PYTHON_MODULE( mod )
{
  pyextra::pyimport("G4Interfaces");

  mod.def("construct_world", &construct_world);
  mod.def("locate", &locate);
  mod.def("world_daughters", &world_daughters);
  mod.def("world_tree", &world_tree);

  //the mask file writer of the masking analysis (LokiMasking)
  py::class_<MaskFileCreator>(mod, "MaskFileCreator")
    .def(py::init([](const std::string& fileName, int indexOffset, const std::vector<int>& bankPixelLimits, int aimingBankId) {
      return new MaskFileCreator(fileName.c_str(), indexOffset, bankPixelLimits, aimingBankId); }))
    .def("isPixelEntered", &MaskFileCreator::isPixelEntered)
    .def("setPixelEntered", &MaskFileCreator::setPixelEntered)
    .def("isPixelEnteredAimingCheck", &MaskFileCreator::isPixelEnteredAimingCheck)
    .def("setPixelEnteredAimingCheck", &MaskFileCreator::setPixelEnteredAimingCheck)
    .def("createMaskFile", &MaskFileCreator::createMaskFile);

  py::class_<PixelCalc>(mod, "PixelCalc")
    .def(py::init<double,int,int>())
    .def("getPixelId", &PixelCalc::getPixelId)
    .def("getPixelId3D", &PixelCalc::getPixelId3D)
    .def("getBankTransform", &PixelCalc::getBankTransform)
    .def("getBankPosition", &PixelCalc::getBankPosition)
    .def("getTotalNumberOfPixels", &PixelCalc::getTotalNumberOfPixels)
    .def("getNumberOfBanks", &PixelCalc::getNumberOfBanks)
    ;

  //static getters (PixelatedBanks / BcsBanks):
  mod.def("getTubeLayerId", &PixelatedBanks::getTubeLayerId);
  mod.def("getNumberOfPixelsInStraw", &PixelatedBanks::getNumberOfPixelsInStraw);
  mod.def("getBankPixelOffset", &PixelatedBanks::getBankPixelOffset);
  mod.def("getBankRotation", &BcsBanks::getBankRotation);
  mod.def("getBankSize", &BcsBanks::getBankSize);
  mod.def("getStrawLengthByBankId", &BcsBanks::getStrawLengthByBankId);
  mod.def("getNumberOfPacksByBankId", &BcsBanks::getNumberOfPacksByBankId);
  mod.def("getNumberOfTubes", &BcsBanks::getNumberOfTubes);
  mod.def("getPackPositionInBank", &BcsBanks::getPackPositionInBank);
  mod.def("getPackRotation", &BcsBanks::getPackRotation);
  mod.def("getPackPackDistance", &BcsBanks::getPackPackDistance);
  mod.def("detectorSystemFrontDistanceFromBankFront", &BcsBanks::detectorSystemFrontDistanceFromBankFront);
  mod.def("getLarmor2022ExperimentBankPositionY", &BcsBanks::getLarmor2022ExperimentBankPositionY);
  mod.def("isVertical", &BcsBanks::isVertical);
  mod.def("areTubesInverselyNumbered", &BcsBanks::areTubesInverselyNumbered);
  mod.def("getBeamstopSize", &BcsBanks::getBeamstopSize);
  //BcsPack / BcsTube:
  mod.def("getTubeRotation", &BcsPack::getTubeRotation);
  mod.def("getHorizontalTubeOffset", &BcsPack::getHorizontalTubeOffset);
  mod.def("getVerticalTubeOffset", &BcsPack::getVerticalTubeOffset);
  mod.def("getStrawPositionX", &BcsTube::getStrawPositionX);
  mod.def("getStrawPositionY", &BcsTube::getStrawPositionY);
  mod.def("getTubeOuterRadius", &BcsTube::getTubeOuterRadius);
  mod.def("getStrawInnerRadius", &BcsTube::getStrawInnerRadius);
}
