// Test-only helper module for the G4GeoLokiTests package (private module used
// by the G4GeoLokiTests.TestUtils Python module, which collects the results
// passed to the "add" callbacks below into lists).
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

#include "sbldpy/sbldpy.hh"
#include "G4Interfaces/GeoConstructBase.hh"
#include "G4Interfaces/ComponentHandles.hh"
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
#include "G4VSolid.hh"
#include "G4Box.hh"
#include <array>
#include <algorithm>
#include <vector>
#include <memory>
#include <iostream>
#include <cstdio>

namespace {

  using I = std::int64_t;
  using V3 = std::array<double,3>;
  using M3 = std::array<V3,3>;
  V3 v3( const G4ThreeVector& v ) { return { v.x(), v.y(), v.z() }; }
  M3 rows( const G4RotationMatrix& r ) { return { V3{ r.xx(), r.xy(), r.xz() }, V3{ r.yx(), r.yy(), r.yz() }, V3{ r.zx(), r.zy(), r.zz() } }; }

  G4VPhysicalVolume* s_world = nullptr;
  std::unique_ptr<G4Navigator> s_nav;

  //Calls Construct() on a geometry module instance (e.g. G4GeoLoki.GeoBCSBanks,
  //given by its handle) and keeps the returned world volume for subsequent locate(..) calls. Note
  //that the GeoConstructBase::place(..) calls perform the standard Geant4
  //overlap check for each placement during construction (printed to G4cout).
  int construct_world( I geo_handle )
  {
    auto geo = &G4Interfaces::getComponentAs<G4Interfaces::GeoConstructBase>( geo_handle, "geometry module" );
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
  //construct_world. Calls add(physvolname, copynumber, logvolname) from the
  //deepest volume up to the world, and returns a tuple:
  //   ( (local x, local y, local z) of the point in the deepest volume,
  //     (dx,dy,dz) global direction of the local z axis of the deepest volume,
  //     (ox,oy,oz) global position of the origin of the deepest volume )
  std::tuple<V3,V3,V3> locate( double x, double y, double z, sbldpy::Callback<void(std::string,I,std::string)> add )
  {
    if (!s_nav)
      throw std::runtime_error("locate: call construct_world first");
    G4ThreeVector p(x,y,z);
    s_nav->LocateGlobalPointAndSetup(p, nullptr, false, true);
    std::unique_ptr<G4TouchableHistory> th(s_nav->CreateTouchableHistory());
    for (int d = 0; d <= th->GetHistoryDepth(); ++d) {
      auto pv = th->GetVolume(d);
      add(std::string(pv->GetName()), th->GetCopyNumber(d), std::string(pv->GetLogicalVolume()->GetName()));
    }
    G4AffineTransform g2l = s_nav->GetGlobalToLocalTransform();
    G4AffineTransform l2g = s_nav->GetLocalToGlobalTransform();
    G4ThreeVector loc = g2l.TransformPoint(p);
    G4ThreeVector ax = l2g.TransformAxis(G4ThreeVector(0,0,1));
    G4ThreeVector org = l2g.TransformPoint(G4ThreeVector(0,0,0));
    return { v3(loc), v3(ax), v3(org) };
  }

  //The volumes placed directly in the world of the geometry built with construct_world:
  //add(physvolname, copynumber, logvolname, frame rotation as 3 rows, translation (x,y,z) [mm]) for each
  //(the frame rotation is what G4PVPlacement stores: the inverse of the object rotation).
  void world_daughters( sbldpy::Callback<void(std::string,I,std::string,M3,V3)> add )
  {
    if (!s_world)
      throw std::runtime_error("world_daughters: call construct_world first");
    auto lv = s_world->GetLogicalVolume();
    for (size_t i = 0; i < lv->GetNoDaughters(); ++i) {
      auto pv = lv->GetDaughter(i);
      const G4RotationMatrix* r = pv->GetRotation();
      G4RotationMatrix identity;
      if (!r) r = &identity;
      const G4ThreeVector t = pv->GetTranslation();
      add(std::string(pv->GetName()), pv->GetCopyNo(), std::string(pv->GetLogicalVolume()->GetName()), rows(*r), v3(t));
    }
  }

  //Dense overlap check of the volumes placed directly in the world (much denser than the
  //construction-time check, and it also reports the overlap depth): npoints random surface
  //points of each world daughter (G4VSolid::GetPointOnSurface) are tested against every other
  //world daughter. A point counts if it is inside the other volume by more than tolerance
  //(DistanceToOut). add(volume, other volume, number of points, max depth [mm]) for the pairs
  //with any such point (volume = "physvolname:copynumber").
  void dense_world_overlaps(I npoints, double tolerance, sbldpy::Callback<void(std::string,std::string,I,double)> add)
  {
    if (!s_world)
      throw std::runtime_error("dense_world_overlaps: call construct_world first");
    auto lv = s_world->GetLogicalVolume();
    const size_t n = lv->GetNoDaughters();
    std::vector<G4AffineTransform> toWorld(n);
    for (size_t i = 0; i < n; ++i) {
      auto pv = lv->GetDaughter(i);
      // G4AffineTransform(R, t) maps v -> R^-1 v + t, so it takes the frame rotation (as in the Geant4 navigation)
      const G4RotationMatrix* frameRotation = pv->GetRotation();
      toWorld[i] = frameRotation ? G4AffineTransform(frameRotation, pv->GetTranslation()) : G4AffineTransform(pv->GetTranslation());
    }
    auto label = [lv](size_t i) {
      auto pv = lv->GetDaughter(i);
      return std::string(pv->GetName()) + ":" + std::to_string(pv->GetCopyNo());
    };
    for (size_t a = 0; a < n; ++a) {
      const G4VSolid* solidA = lv->GetDaughter(a)->GetLogicalVolume()->GetSolid();
      std::vector<G4ThreeVector> points(npoints);
      for (auto& p : points)
        p = toWorld[a].TransformPoint(solidA->GetPointOnSurface());
      for (size_t b = 0; b < n; ++b) {
        if (a == b)
          continue;
        const G4VSolid* solidB = lv->GetDaughter(b)->GetLogicalVolume()->GetSolid();
        const G4AffineTransform toB = toWorld[b].Inverse();
        int count = 0;
        double maxDepth = 0.;
        for (const auto& p : points) {
          const G4ThreeVector q = toB.TransformPoint(p);
          if (solidB->Inside(q) != kInside)
            continue;
          const double depth = solidB->DistanceToOut(q);
          if (depth <= tolerance)
            continue;
          ++count;
          maxDepth = std::max(maxDepth, depth);
        }
        if (count)
          add(label(a), label(b), count, maxDepth);
      }
    }
  }

  //Dense containment check of the volumes inside the banks (packs, boron masks, calibration masks...):
  //npoints random surface points of each daughter of every bank must be inside its bank solid (which
  //for the banks 6 and 8 has the notch cut out). add(bank copy number, "daughter:copy", number of points
  //outside by more than tolerance, max distance [mm]) for the daughters with such points.
  void dense_bank_containment(I npoints, double tolerance, sbldpy::Callback<void(I,std::string,I,double)> add)
  {
    if (!s_world)
      throw std::runtime_error("dense_bank_containment: call construct_world first");
    auto world = s_world->GetLogicalVolume();
    for (size_t i = 0; i < world->GetNoDaughters(); ++i) {
      auto bankPV = world->GetDaughter(i);
      auto bankLV = bankPV->GetLogicalVolume();
      if (bankLV->GetName() != "Bank")
        continue;
      const G4VSolid* bankSolid = bankLV->GetSolid();
      for (size_t j = 0; j < bankLV->GetNoDaughters(); ++j) {
        auto pv = bankLV->GetDaughter(j);
        const G4RotationMatrix* frameRotation = pv->GetRotation();
        const G4AffineTransform toBank = frameRotation ? G4AffineTransform(frameRotation, pv->GetTranslation()) : G4AffineTransform(pv->GetTranslation());
        const G4VSolid* solid = pv->GetLogicalVolume()->GetSolid();
        int count = 0;
        double maxDistance = 0.;
        for (int k = 0; k < npoints; ++k) {
          const G4ThreeVector q = toBank.TransformPoint(solid->GetPointOnSurface());
          if (bankSolid->Inside(q) != kOutside)
            continue;
          const double distance = bankSolid->DistanceToIn(q);
          if (distance <= tolerance)
            continue;
          ++count;
          maxDistance = std::max(maxDistance, distance);
        }
        if (count)
          add(bankPV->GetCopyNo(), std::string(pv->GetName()) + ":" + std::to_string(pv->GetCopyNo()), count, maxDistance);
      }
    }
  }

  //Smallest distance [mm] of the daughters of a bank (packs, masks) from a box given in the bank frame
  //(e.g. the notch): npoints random surface points per daughter. (distance, "daughter:copy" of the closest).
  std::tuple<double,std::string> bank_daughters_distance_to_box(I bankCopyNo, V3 centre, V3 halfSize, I npoints)
  {
    if (!s_world)
      throw std::runtime_error("bank_daughters_distance_to_box: call construct_world first");
    auto world = s_world->GetLogicalVolume();
    const G4Box box("box", halfSize[0], halfSize[1], halfSize[2]);
    const G4ThreeVector boxCentre(centre[0], centre[1], centre[2]);
    double minDistance = 1e99;
    std::string closest;
    for (size_t i = 0; i < world->GetNoDaughters(); ++i) {
      auto bankPV = world->GetDaughter(i);
      auto bankLV = bankPV->GetLogicalVolume();
      if (bankLV->GetName() != "Bank" || bankPV->GetCopyNo() != bankCopyNo)
        continue;
      for (size_t j = 0; j < bankLV->GetNoDaughters(); ++j) {
        auto pv = bankLV->GetDaughter(j);
        const G4RotationMatrix* frameRotation = pv->GetRotation();
        const G4AffineTransform toBank = frameRotation ? G4AffineTransform(frameRotation, pv->GetTranslation()) : G4AffineTransform(pv->GetTranslation());
        const G4VSolid* solid = pv->GetLogicalVolume()->GetSolid();
        for (int k = 0; k < npoints; ++k) {
          const G4ThreeVector q = toBank.TransformPoint(solid->GetPointOnSurface()) - boxCentre;
          const double distance = box.Inside(q) == kOutside ? box.DistanceToIn(q) : 0.0;
          if (distance < minDistance) {
            minDistance = distance;
            closest = std::string(pv->GetName()) + ":" + std::to_string(pv->GetCopyNo());
          }
        }
      }
    }
    return { minDistance, closest };
  }

  //Every physical volume of the geometry built with construct_world, depth first:
  //add(depth, physvolname, copynumber, logvolname, global translation (x,y,z) [mm],
  //  global object rotation as 3 rows) for each (depth 0 = the world).
  using TreeAdd = sbldpy::Callback<void(I,std::string,I,std::string,V3,M3)>;
  void addTree(TreeAdd& add, const G4VPhysicalVolume* pv, int depth, const G4RotationMatrix& Rp, const G4ThreeVector& tp)
  {
    const G4RotationMatrix R = Rp * pv->GetObjectRotationValue();
    const G4ThreeVector t = tp + Rp * pv->GetObjectTranslation();
    add(depth, std::string(pv->GetName()), pv->GetCopyNo(), std::string(pv->GetLogicalVolume()->GetName()), v3(t), rows(R));
    auto lv = pv->GetLogicalVolume();
    for (size_t i = 0; i < lv->GetNoDaughters(); ++i)
      addTree(add, lv->GetDaughter(i), depth + 1, R, t);
  }
  void world_tree( TreeAdd add )
  {
    if (!s_world)
      throw std::runtime_error("world_tree: call construct_world first");
    addTree(add, s_world, 0, G4RotationMatrix(), G4ThreeVector());
  }

  //Thin wrapper around PixelatedBanks (the class used by the analysis programs
  //to convert hits to pixel ids).
  struct PixelCalc {
    PixelatedBanks pb;
    PixelCalc(double rear, I n, I nbanks, const std::string& bankCalibration) : pb(rear, n, nbanks, bankCalibration) {}
    I getPixelId(I bank, I tube, I straw, double x, double y, double z) const { return pb.getPixelId(bank,tube,straw,x,y,z); }
    //bank transform: (rotation as 3 rows, translation)
    std::pair<M3,V3> getBankTransform(I bank) const {
      const BankTransform t = pb.getBankTransform(bank);
      return { t.rotation, t.translation };
    }
    std::pair<M3,V3> getNominalBankTransform(I bank) const {
      const BankTransform t = pb.getNominalBankTransform(bank);
      return { t.rotation, t.translation };
    }
    //the notch of bank 6 / 8: (centre, half sizes, isNominal)
    std::tuple<V3,V3,bool> getBankNotch(I bank) const {
      const auto notch = pb.getBankNotch(bank);
      return { notch.centre, notch.halfSize, notch.isNominal };
    }
    //the points spanning the overlap region of two bank boxes, in the frame of
    //the first bank (as a flat array, x0,y0,z0,x1,...)
    sbldpy::NumpyCopy<double> getBankBoxOverlapPoints(I bank, I otherBank) const {
      std::vector<double> v;
      for (const auto& p : pb.getBankBoxOverlapPoints(bank, otherBank))
        v.insert(v.end(), p.begin(), p.end());
      return sbldpy::NumpyCopy<double>( std::move( v ) );
    }
    bool isBankCalibrated(I bank) const { return pb.isBankCalibrated(bank); }
    std::string getBankCalibrationName() const { return pb.getBankCalibration().name(); }
    std::string getBankCalibrationFile() const { return pb.getBankCalibration().fileName(); }
    I getTotalNumberOfPixels() const { return pb.getTotalNumberOfPixels(); }
    I getNumberOfBanks() const { return pb.getNumberOfBanks(); }
  };

  //The mask file writer of the masking analysis (LokiMasking), with the bank
  //pixel limits given as a numpy array:
  struct MaskFileCreatorImpl {
    MaskFileCreator mfc;
    MaskFileCreatorImpl( const std::string& fileName, I indexOffset, sbldpy::NumpyView<const std::int64_t> bankPixelLimits, I aimingBankId )
      : mfc( fileName.c_str(), indexOffset, std::vector<int>( bankPixelLimits.begin(), bankPixelLimits.end() ), aimingBankId ) {}
  };
}

SBLDPY_MODULE( _TestUtils )
{
  m.def("construct_world", &construct_world, { "geo_handle" });
  m.def("locate", &locate, { "x", "y", "z", "add" });
  m.def("world_daughters", &world_daughters, { "add" });
  m.def("world_tree", &world_tree, { "add" });
  m.def("dense_world_overlaps", &dense_world_overlaps, { "npoints", "tolerance", "add" });
  m.def("dense_bank_containment", &dense_bank_containment, { "npoints", "tolerance", "add" });
  m.def("bank_daughters_distance_to_box", &bank_daughters_distance_to_box, { "bankCopyNo", "centre", "halfSize", "npoints" });

  using MFC = MaskFileCreatorImpl;
  m.add_class<MFC>("MaskFileCreator")
    .def_init<std::string, I, sbldpy::NumpyView<const std::int64_t>, I>( { "fileName", "indexOffset", "bankPixelLimits", "aimingBankId" } )
    .def("isPixelEntered", []( const MFC& o, I pixel ) { return o.mfc.isPixelEntered( pixel ); }, { "pixelNumber" })
    .def("setPixelEntered", []( MFC& o, I pixel ) { o.mfc.setPixelEntered( pixel ); }, { "pixelNumber" })
    .def("isPixelEnteredAimingCheck", []( const MFC& o, I pixel ) { return o.mfc.isPixelEnteredAimingCheck( pixel ); }, { "pixelNumber" })
    .def("setPixelEnteredAimingCheck", []( MFC& o, I pixel ) { o.mfc.setPixelEnteredAimingCheck( pixel ); }, { "pixelNumber" })
    .def("createMaskFile", []( const MFC& o ) { o.mfc.createMaskFile(); });

  m.add_class<PixelCalc>("PixelCalc")
    .def_init<double,I,I,std::string>( { "rear", "n", "nbanks", "bankCalibration" } )
    .def("getPixelId", &PixelCalc::getPixelId, { "bank", "tube", "straw", "x", "y", "z" })
    .def("getBankTransform", &PixelCalc::getBankTransform, { "bank" })
    .def("getNominalBankTransform", &PixelCalc::getNominalBankTransform, { "bank" })
    .def("getBankNotch", &PixelCalc::getBankNotch, { "bank" })
    .def("getBankBoxOverlapPoints", &PixelCalc::getBankBoxOverlapPoints, { "bank", "otherBank" })
    .def("isBankCalibrated", &PixelCalc::isBankCalibrated, { "bank" })
    .def("getBankCalibrationName", &PixelCalc::getBankCalibrationName)
    .def("getBankCalibrationFile", &PixelCalc::getBankCalibrationFile)
    .def("getTotalNumberOfPixels", &PixelCalc::getTotalNumberOfPixels)
    .def("getNumberOfBanks", &PixelCalc::getNumberOfBanks);

  //static getters (PixelatedBanks / BcsBanks):
  m.def("getTubeLayerId", []( I bank, I tube ) { return BcsBanks::getTubeLayerId( bank, tube ); }, { "bankId", "tubeId" });
  m.def("getTubeIdInBank", []( I bank, I pack, I inPackTube ) { return BcsBanks::getTubeIdInBank( bank, pack, inPackTube ); },
        { "bankId", "packId", "inPackTubeId" });
  m.def("getPackId", []( I bank, I tube ) { return BcsBanks::getPackId( bank, tube ); }, { "bankId", "tubeId" });
  m.def("getInPackTubeId", []( I bank, I tube ) { return BcsBanks::getInPackTubeId( bank, tube ); }, { "bankId", "tubeId" });
  m.def("getBankSize", []( I bank, I axis ) { return BcsBanks::getBankSize( bank, axis ); }, { "bankId", "axisIndex" });
  m.def("getStrawLengthByBankId", []( I bank ) { return BcsBanks::getStrawLengthByBankId( bank ); }, { "bankId" });
  m.def("getNumberOfPacksByBankId", []( I bank ) { return BcsBanks::getNumberOfPacksByBankId( bank ); }, { "bankId" });
  m.def("getNumberOfTubes", []( I bank ) { return BcsBanks::getNumberOfTubes( bank ); }, { "bankId" });
  m.def("getPackPositionInBank", []( I bank, I pack, I axis ) { return BcsBanks::getPackPositionInBank( bank, pack, axis ); },
        { "bankId", "packNumber", "axisIndex" });
  m.def("getPackRotation", &BcsBanks::getPackRotation);
  m.def("getPackPackDistance", &BcsBanks::getPackPackDistance);
  m.def("detectorSystemFrontDistanceFromBankFront", []( I bank ) { return BcsBanks::detectorSystemFrontDistanceFromBankFront( bank ); },
        { "bankId" });
  m.def("isVertical", []( I bank ) { return BcsBanks::isVertical( bank ); }, { "bankId" });
  m.def("areTubesInverselyNumbered", []( I bank ) { return BcsBanks::areTubesInverselyNumbered( bank ); }, { "bankId" });
  m.def("getBeamstopSize", []( I beamstop, I axis ) { return BcsBanks::getBeamstopSize( beamstop, axis ); }, { "beamstopId", "axisIndex" });
  m.def("getFrontFaceCentreInBank", []( I bank ) { return BcsBanks::getFrontFaceCentreInBank( bank ); }, { "bankId" });
  //BcsPack / BcsTube:
  m.def("getTubeRotation", &BcsPack::getTubeRotation);
  m.def("getHorizontalTubeOffset", []( I tube ) { return BcsPack::getHorizontalTubeOffset( tube ); }, { "inPackTubeId" });
  m.def("getVerticalTubeOffset", []( I tube ) { return BcsPack::getVerticalTubeOffset( tube ); }, { "inPackTubeId" });
  m.def("getStrawPositionX", []( I straw ) { return BcsTube::getStrawPositionX( straw ); }, { "strawId" });
  m.def("getStrawPositionY", []( I straw ) { return BcsTube::getStrawPositionY( straw ); }, { "strawId" });
  m.def("getTubeOuterRadius", &BcsTube::getTubeOuterRadius);
  m.def("getStrawInnerRadius", &BcsTube::getStrawInnerRadius);
}
