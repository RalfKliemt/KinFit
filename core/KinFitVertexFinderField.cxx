//****************************************************************************
//*                    This file is part of KinFit.                          *
//*                                                                          *
//*             KinFit is distributed under the terms of the                 *
//*              GNU General Public License (GPL) version 3,                 *
//*                 copied verbatim in the file "LICENSE".                   *
//*                                                                          *
//*                             Copyright 2024                               *
//*             GSI Helmholtzzentrum für Schwerionenforschung                *
//*          This software is distributed under the terms of the             *
//*          GNU General Public Licence version 3 (GPL Version 3)            *
//*                                                                          *
//*     The copyright holders are listed in the file "COPYRIGHTHOLDERS".     *
//*               The authors are listed in the file "AUTHORS".              *
//****************************************************************************

#include "KinFitVertexFinderField.h"
#include "KinFitTrajectory.h"

#include <TMatrixDSym.h>
#include <TVectorD.h>

#include <cmath>

namespace {
/** Least-squares point closest to lines (r_i, unit t_i):
 * sum_i (I - t t^T) (v - r_i) = 0. */
bool solveLines(const std::vector<TVector3> &r, const std::vector<TVector3> &t,
                TVector3 &v) {
  TMatrixDSym A(3);
  TVectorD b(3);
  for (size_t i = 0; i < r.size(); ++i) {
    const double tt[3] = {t[i].X(), t[i].Y(), t[i].Z()};
    const double rr[3] = {r[i].X(), r[i].Y(), r[i].Z()};
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        const double m = (j == k ? 1. : 0.) - tt[j] * tt[k];
        A(j, k) += m;
        b(j) += m * rr[k];
      }
  }
  const double det = A.Determinant();
  const double scale = (A(0, 0) + A(1, 1) + A(2, 2)) / 3.;
  if (!(std::abs(det) > 1e-12 * scale * scale * scale))
    return false;
  A.Invert();
  const TVectorD x = A * b;
  v.SetXYZ(x(0), x(1), x(2));
  return true;
}
} // namespace

KinFitVertexFinderField::KinFitVertexFinderField(
    const std::vector<KinFitParticle> &cands, const KinFitField *field,
    int maxIterations, double tolerance) {
  std::vector<KinFitTrajectory> traj;
  for (const KinFitParticle &c : cands) {
    if (c.isNeutralHit())
      continue;
    traj.push_back(KinFitTrajectory::fromTrackParameters(
        1. / c.P(), c.Theta(), c.Phi(), c.getR(), c.getZ(), c.getCharge(),
        field));
  }
  fNTracks = traj.size();
  if (fNTracks < 2)
    return;

  // Start: tangent lines at the reference states
  std::vector<TVector3> r(fNTracks), t(fNTracks);
  for (int i = 0; i < fNTracks; ++i) {
    r[i] = traj[i].getReferencePosition();
    t[i] = traj[i].getReferenceDirection();
  }
  TVector3 v;
  if (!solveLines(r, t, v)) {
    v.SetXYZ(0, 0, 0);
    for (const TVector3 &ri : r)
      v += ri;
    v *= 1. / fNTracks;
  }

  for (fIterations = 1; fIterations <= maxIterations; ++fIterations) {
    for (int i = 0; i < fNTracks; ++i)
      traj[i].closestApproach(v, r[i], t[i]);
    TVector3 vNew;
    if (!solveLines(r, t, vNew))
      return;
    const double shift = (vNew - v).Mag();
    v = vNew;
    if (shift < tolerance) {
      fVertex = v;
      fValid = true;
      return;
    }
  }
  // not converged: keep the last estimate but flag it
  fVertex = v;
}
