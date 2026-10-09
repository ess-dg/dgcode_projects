#include "G4GeoLoki/BcsBanks.hh"
#include "Units/Units.hh"
#include <cmath>
#include <array>
#include <cassert>
#include <algorithm>
#include <vector>

// The notch of the banks 6 and 8, and the overlap region of two bank boxes it is computed from (see BcsBanks.hh).

/// bank notch ///
const double BcsBanks::nominalNotchCentre[3] = {-145.0 *Units::mm, -42.0 *Units::mm, -425.0 *Units::mm};
const double BcsBanks::nominalNotchHalfSize[3] = {15.0 *Units::mm, 70.0 *Units::mm, 60.0 *Units::mm};
const double BcsBanks::notchMargin = 1.0 *Units::mm;

bool BcsBanks::hasBankNotch(const int bankId) {
  return bankId == 6 || bankId == 8;
}

int BcsBanks::getBankNotchNeighbour(const int bankId) {
  assert(hasBankNotch(bankId));
  return bankId - 1;
}

namespace {
  struct Box { // oriented box: p_world = rotation * p_box + centre, |p_box[k]| <= halfSize[k]
    BankTransform transform;
    std::array<double,3> halfSize;
  };

  // The part of the segment p0-p1 inside the box (Liang-Barsky clipping); false if it is outside.
  bool clipSegment(const std::array<double,3>& p0, const std::array<double,3>& p1, const Box& box,
                   std::array<double,3>& q0, std::array<double,3>& q1) {
    const auto a = box.transform.toLocal(p0);
    const auto b = box.transform.toLocal(p1);
    double tMin = 0.0, tMax = 1.0;
    for (int k = 0; k < 3; k++) {
      const double d = b[k] - a[k];
      if (std::abs(d) < 1e-12) {
        if (std::abs(a[k]) > box.halfSize[k])
          return false;
        continue;
      }
      double t1 = (-box.halfSize[k] - a[k]) / d;
      double t2 = (box.halfSize[k] - a[k]) / d;
      if (t1 > t2)
        std::swap(t1, t2);
      tMin = std::max(tMin, t1);
      tMax = std::min(tMax, t2);
      if (tMin > tMax)
        return false;
    }
    for (int i = 0; i < 3; i++) {
      q0[i] = p0[i] + tMin * (p1[i] - p0[i]);
      q1[i] = p0[i] + tMax * (p1[i] - p0[i]);
    }
    return true;
  }

  // The parts of the 12 edges of box a that are inside box b
  void addEdgePointsInside(const Box& a, const Box& b, std::vector<std::array<double,3>>& points) {
    for (int corner = 0; corner < 8; corner++) {
      for (int axis = 0; axis < 3; axis++) {
        if (corner & (1 << axis))
          continue; // each edge once: from the corner with the lower coordinate along the axis
        std::array<double,3> c0, c1;
        for (int k = 0; k < 3; k++) {
          c0[k] = ((corner >> k) & 1) ? a.halfSize[k] : -a.halfSize[k];
          c1[k] = c0[k];
        }
        c1[axis] = a.halfSize[axis];
        std::array<double,3> q0, q1;
        if (clipSegment(a.transform.toGlobal(c0), a.transform.toGlobal(c1), b, q0, q1)) {
          points.push_back(q0);
          points.push_back(q1);
        }
      }
    }
  }
}

std::vector<std::array<double,3>> BcsBanks::getBankBoxOverlapPoints(const int bankId, const int otherBankId) const {
  // The overlap region of two boxes is a convex polyhedron; its vertices are on the edges of the two boxes, so
  // it is spanned by the parts of the edges of each box that are inside the other box.
  auto bankBox = [this](const int id) {
    return Box{getBankTransform(id), {0.5 * getBankSize(id, 2), 0.5 * getBankSize(id, 1), 0.5 * getBankSize(id, 0)}};
  };
  const Box a = bankBox(bankId);
  const Box b = bankBox(otherBankId);
  std::vector<std::array<double,3>> points;
  addEdgePointsInside(a, b, points);
  addEdgePointsInside(b, a, points);
  for (auto& point : points)
    point = a.transform.toLocal(point);
  return points;
}

BcsBanks::BankNotch BcsBanks::getBankNotch(const int bankId) const {
  assert(hasBankNotch(bankId));
  BankNotch notch;
  notch.isNominal = true;
  for (int k = 0; k < 3; k++) {
    notch.centre[k] = nominalNotchCentre[k];
    notch.halfSize[k] = nominalNotchHalfSize[k];
  }
  const int neighbour = getBankNotchNeighbour(bankId);
  if (!isBankCalibrated(bankId) && !isBankCalibrated(neighbour))
    return notch; // nominal placements
  const auto points = getBankBoxOverlapPoints(bankId, neighbour);
  bool insideNominalNotch = true;
  for (const auto& point : points)
    for (int k = 0; k < 3; k++)
      if (std::abs(point[k] - nominalNotchCentre[k]) > nominalNotchHalfSize[k])
        insideNominalNotch = false;
  if (insideNominalNotch)
    return notch; // (also if the boxes don't overlap)
  for (int k = 0; k < 3; k++) {
    double low = nominalNotchCentre[k] - nominalNotchHalfSize[k];
    double high = nominalNotchCentre[k] + nominalNotchHalfSize[k];
    for (const auto& point : points) {
      low = std::min(low, point[k] - notchMargin);
      high = std::max(high, point[k] + notchMargin);
    }
    notch.centre[k] = 0.5 * (low + high);
    notch.halfSize[k] = 0.5 * (high - low);
  }
  notch.isNominal = false;
  return notch;
}
