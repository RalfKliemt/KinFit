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

#include "KinFitTrajectory.h"
#include "KinFitField.h"

#include <algorithm>
#include <cmath>

KinFitTrajectory::KinFitTrajectory(const TVector3 &pos, const TVector3 &mom,
                                   double charge, const KinFitField *field)
    : fPos(pos), fDir(mom.Unit()), fP(mom.Mag()), fCharge(charge),
      fField(field) {
  fScale = (fP > 0) ? kKappa * fCharge / fP : 0.;
  fUniform = !fField || fField->isUniform() || fCharge == 0;
  if (fUniform) {
    const TVector3 B = (fField && fCharge != 0) ? fField->getField(fPos)
                                                : TVector3(0, 0, 0);
    const double b = B.Mag();
    if (b > 0 && fScale != 0) {
      fAxis = B * (1. / b);
      fOmega = fScale * b;
    }
  }
}

KinFitTrajectory KinFitTrajectory::fromTrackParameters(
    double invP, double theta, double phi, double R, double Z, double charge,
    const KinFitField *field) {
  const double sp = std::sin(phi), cp = std::cos(phi);
  const double st = std::sin(theta), ct = std::cos(theta);
  const double p = 1. / invP;
  return KinFitTrajectory(TVector3(-R * sp, R * cp, Z),
                          TVector3(p * st * cp, p * st * sp, p * ct), charge,
                          field);
}

void KinFitTrajectory::propagate(double s, TVector3 &pos, TVector3 &dir,
                                 TVector3 *dDirDs) const {
  if (fUniform)
    propagateHelix(s, pos, dir, dDirDs);
  else
    propagateRK(s, pos, dir, dDirDs);
}

void KinFitTrajectory::propagateHelix(double s, TVector3 &pos, TVector3 &dir,
                                      TVector3 *dDirDs) const {
  if (fOmega == 0) {
    pos = fPos + s * fDir;
    dir = fDir;
    if (dDirDs)
      dDirDs->SetXYZ(0, 0, 0);
    return;
  }
  // t(s) = t_par b + t_perp cos(ws) + (t_perp x b) sin(ws)
  // r(s) = r0 + t_par b s + t_perp sin(ws)/w + (t_perp x b)(1 - cos(ws))/w
  const double tPar = fDir.Dot(fAxis);
  const TVector3 tPerp = fDir - tPar * fAxis;
  const TVector3 tCross = tPerp.Cross(fAxis);
  const double ws = fOmega * s;
  const double c = std::cos(ws), sn = std::sin(ws);
  const double sinOverW = (std::abs(ws) < 1e-8) ? s : sn / fOmega;
  const double oneMinusCosOverW =
      2. * std::sin(0.5 * ws) * std::sin(0.5 * ws) / fOmega;
  pos = fPos + (tPar * s) * fAxis + sinOverW * tPerp +
        oneMinusCosOverW * tCross;
  dir = tPar * fAxis + c * tPerp + sn * tCross;
  if (dDirDs)
    *dDirDs = fOmega * (-sn * tPerp + c * tCross);
}

void KinFitTrajectory::propagateRK(double s, TVector3 &pos, TVector3 &dir,
                                   TVector3 *dDirDs) const {
  // Step size: at most fMaxStep and 1/20 of the radius at the start. The
  // number of steps depends only on |s| (smooth in the track parameters).
  double h = fMaxStep;
  const double curv = std::abs(fScale) * fField->getField(fPos).Mag();
  if (curv > 0)
    h = std::min(h, 0.05 / curv);
  const int n = std::max(1, int(std::ceil(std::abs(s) / h)));
  const double ds = s / n;

  auto accel = [&](const TVector3 &r, const TVector3 &t) {
    return fScale * t.Cross(fField->getField(r));
  };
  TVector3 r = fPos, t = fDir;
  for (int i = 0; i < n; ++i) {
    const TVector3 k1r = t, k1t = accel(r, t);
    const TVector3 r2 = r + 0.5 * ds * k1r, t2 = t + 0.5 * ds * k1t;
    const TVector3 k2r = t2, k2t = accel(r2, t2);
    const TVector3 r3 = r + 0.5 * ds * k2r, t3 = t + 0.5 * ds * k2t;
    const TVector3 k3r = t3, k3t = accel(r3, t3);
    const TVector3 r4 = r + ds * k3r, t4 = t + ds * k3t;
    const TVector3 k4r = t4, k4t = accel(r4, t4);
    r += (ds / 6.) * (k1r + 2. * k2r + 2. * k3r + k4r);
    t += (ds / 6.) * (k1t + 2. * k2t + 2. * k3t + k4t);
    t = t.Unit();
  }
  pos = r;
  dir = t;
  if (dDirDs)
    *dDirDs = accel(r, t);
}

double KinFitTrajectory::closestApproach(const TVector3 &point, TVector3 &pos,
                                         TVector3 &dir) const {
  // Newton on g(s) = (r(s) - P) . t(s), g'(s) = 1 + (r - P) . dt/ds, from
  // the straight-line estimate. Steps are limited to a fraction of the
  // radius (stay on the same turn) and only taken where g' > 0 (minimum).
  double s = (point - fPos).Dot(fDir);
  double maxStep = 1e30;
  if (fUniform && fOmega != 0)
    maxStep = 0.5 / std::abs(fOmega);
  else if (!fUniform) {
    const double curv = std::abs(fScale) * fField->getField(point).Mag();
    if (curv > 0)
      maxStep = 0.5 / curv;
  }
  TVector3 dt;
  for (int iter = 0; iter < 50; ++iter) {
    propagate(s, pos, dir, &dt);
    const TVector3 diff = pos - point;
    const double g = diff.Dot(dir);
    double gp = 1. + diff.Dot(dt);
    if (gp < 0.1)
      gp = 1.; // not convex here: plain gradient step
    const double step = std::clamp(g / gp, -maxStep, maxStep);
    s -= step;
    if (std::abs(step) < 1e-10 * (1. + std::abs(s)))
      break;
  }
  propagate(s, pos, dir);
  return s;
}

bool KinFitTrajectory::trackParameters(double &invP, double &theta,
                                       double &phi, double &R,
                                       double &Z) const {
  // Newton on h(s) = r_T(s) . t_T(s), h' = |t_T|^2 + r_T . (dt/ds)_T
  TVector3 pos, dir, dt;
  double s = 0;
  const double tT2 = fDir.Perp2();
  if (tT2 > 0)
    s = -(fPos.X() * fDir.X() + fPos.Y() * fDir.Y()) / tT2;
  bool ok = tT2 == 0;
  for (int iter = 0; iter < 50 && !ok; ++iter) {
    propagate(s, pos, dir, &dt);
    const double h = pos.X() * dir.X() + pos.Y() * dir.Y();
    double hp = dir.Perp2() + pos.X() * dt.X() + pos.Y() * dt.Y();
    if (!(hp > 0))
      hp = dir.Perp2();
    if (!(hp > 0))
      break;
    const double step = h / hp;
    s -= step;
    ok = std::abs(step) < 1e-10 * (1. + std::abs(s));
  }
  propagate(s, pos, dir);
  invP = 1. / fP;
  theta = dir.Theta();
  phi = dir.Phi();
  R = -pos.X() * std::sin(phi) + pos.Y() * std::cos(phi);
  Z = pos.Z();
  return ok;
}
