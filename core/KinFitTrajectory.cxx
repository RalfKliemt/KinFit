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

namespace {
using Vec = KinFitTrajectory::Vec;
inline Vec add(const Vec &a, const Vec &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
inline Vec sub(const Vec &a, const Vec &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
inline Vec mul(double f, const Vec &a) { return {f * a[0], f * a[1], f * a[2]}; }
inline double dot(const Vec &a, const Vec &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline Vec cross(const Vec &a, const Vec &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}
inline Vec unit(const Vec &a) {
  const double n = std::sqrt(dot(a, a));
  return n > 0 ? mul(1. / n, a) : a;
}
inline Vec fromT(const TVector3 &v) { return {v.X(), v.Y(), v.Z()}; }
inline void toT(const Vec &a, TVector3 &v) { v.SetXYZ(a[0], a[1], a[2]); }
} // namespace

KinFitTrajectory::KinFitTrajectory(const TVector3 &pos, const TVector3 &mom,
                                   double charge, const KinFitField *field)
    : fPos(fromT(pos)), fDir(unit(fromT(mom))), fP(mom.Mag()),
      fCharge(charge), fField(field) {
  init();
}

KinFitTrajectory KinFitTrajectory::fromTrackParameters(
    double invP, double theta, double phi, double R, double Z, double charge,
    const KinFitField *field) {
  const double sp = std::sin(phi), cp = std::cos(phi);
  const double st = std::sin(theta), ct = std::cos(theta);
  KinFitTrajectory t;
  t.fPos = {-R * sp, R * cp, Z};
  t.fDir = {st * cp, st * sp, ct};
  t.fP = 1. / invP;
  t.fCharge = charge;
  t.fField = field;
  t.init();
  return t;
}

KinFitTrajectory::Vec KinFitTrajectory::field(const Vec &r) const {
  return fromT(fField->getField(r[0], r[1], r[2]));
}

void KinFitTrajectory::init() {
  fScale = (fP > 0) ? kKappa * fCharge / fP : 0.;
  fUniform = !fField || fField->isUniform() || fCharge == 0;
  fOmega = 0;
  if (fUniform && fField && fCharge != 0) {
    const Vec B = field(fPos);
    const double b = std::sqrt(dot(B, B));
    if (b > 0 && fScale != 0) {
      fAxis = mul(1. / b, B);
      fOmega = fScale * b;
      fTPar = dot(fDir, fAxis);
      fTPerp = sub(fDir, mul(fTPar, fAxis));
      fTCross = cross(fTPerp, fAxis);
    }
  }
}

void KinFitTrajectory::propagateRaw(double s, Vec &pos, Vec &dir,
                                    Vec *dDirDs) const {
  if (fUniform) {
    if (fOmega == 0) {
      pos = add(fPos, mul(s, fDir));
      dir = fDir;
      if (dDirDs)
        *dDirDs = {0., 0., 0.};
      return;
    }
    // t(s) = t_par b + t_perp cos(ws) + (t_perp x b) sin(ws)
    // r(s) = r0 + t_par b s + t_perp sin(ws)/w + (t_perp x b)(1 - cos(ws))/w
    const double ws = fOmega * s;
    const double c = std::cos(ws), sn = std::sin(ws);
    const double sinOverW = (std::abs(ws) < 1e-8) ? s : sn / fOmega;
    const double sh = std::sin(0.5 * ws);
    const double oneMinusCosOverW = 2. * sh * sh / fOmega;
    for (int i = 0; i < 3; ++i) {
      pos[i] = fPos[i] + fTPar * s * fAxis[i] + sinOverW * fTPerp[i] +
               oneMinusCosOverW * fTCross[i];
      dir[i] = fTPar * fAxis[i] + c * fTPerp[i] + sn * fTCross[i];
    }
    if (dDirDs)
      for (int i = 0; i < 3; ++i)
        (*dDirDs)[i] = fOmega * (-sn * fTPerp[i] + c * fTCross[i]);
    return;
  }

  // Runge-Kutta. Step size: at most fMaxStep and 1/20 of the radius at the
  // start; the number of steps depends only on |s| (smooth in the track
  // parameters).
  double h = fMaxStep;
  const Vec B0 = field(fPos);
  const double curv = std::abs(fScale) * std::sqrt(dot(B0, B0));
  if (curv > 0)
    h = std::min(h, 0.05 / curv);
  const int n = std::max(1, int(std::ceil(std::abs(s) / h)));
  const double ds = s / n;
  auto accel = [&](const Vec &r, const Vec &t) {
    return mul(fScale, cross(t, field(r)));
  };
  Vec r = fPos, t = fDir;
  for (int i = 0; i < n; ++i) {
    const Vec k1r = t, k1t = accel(r, t);
    const Vec r2 = add(r, mul(0.5 * ds, k1r)), t2 = add(t, mul(0.5 * ds, k1t));
    const Vec k2r = t2, k2t = accel(r2, t2);
    const Vec r3 = add(r, mul(0.5 * ds, k2r)), t3 = add(t, mul(0.5 * ds, k2t));
    const Vec k3r = t3, k3t = accel(r3, t3);
    const Vec r4 = add(r, mul(ds, k3r)), t4 = add(t, mul(ds, k3t));
    const Vec k4r = t4, k4t = accel(r4, t4);
    for (int j = 0; j < 3; ++j) {
      r[j] += ds / 6. * (k1r[j] + 2. * k2r[j] + 2. * k3r[j] + k4r[j]);
      t[j] += ds / 6. * (k1t[j] + 2. * k2t[j] + 2. * k3t[j] + k4t[j]);
    }
    t = unit(t);
  }
  pos = r;
  dir = t;
  if (dDirDs)
    *dDirDs = accel(r, t);
}

void KinFitTrajectory::propagate(double s, TVector3 &pos, TVector3 &dir,
                                 TVector3 *dDirDs) const {
  Vec p, d, dt;
  propagateRaw(s, p, d, dDirDs ? &dt : nullptr);
  toT(p, pos);
  toT(d, dir);
  if (dDirDs)
    toT(dt, *dDirDs);
}

double KinFitTrajectory::closestApproach(const TVector3 &point, TVector3 &pos,
                                         TVector3 &dir) const {
  return closestApproach(point, pos, dir,
                         std::numeric_limits<double>::quiet_NaN());
}

double KinFitTrajectory::closestApproach(const TVector3 &point, TVector3 &pos,
                                         TVector3 &dir, double sStart) const {
  Vec p, d;
  const double s = closestApproach(fromT(point), p, d, sStart);
  toT(p, pos);
  toT(d, dir);
  return s;
}

double KinFitTrajectory::closestApproach(const Vec &point, Vec &pos, Vec &dir,
                                         double sStart) const {
  // Newton on g(s) = (r(s) - P) . t(s), g'(s) = 1 + (r - P) . dt/ds. Steps
  // are limited to a fraction of the radius (stay on the same turn) and
  // only taken where g' > 0 (minimum).
  double s = std::isfinite(sStart) ? sStart : dot(sub(point, fPos), fDir);
  double maxStep = 1e30;
  if (fUniform && fOmega != 0)
    maxStep = 0.5 / std::abs(fOmega);
  else if (!fUniform) {
    const Vec B = field(point);
    const double curv = std::abs(fScale) * std::sqrt(dot(B, B));
    if (curv > 0)
      maxStep = 0.5 / curv;
  }
  Vec dt;
  for (int iter = 0; iter < 50; ++iter) {
    propagateRaw(s, pos, dir, &dt);
    const Vec diff = sub(pos, point);
    const double g = dot(diff, dir);
    double gp = 1. + dot(diff, dt);
    if (gp < 0.1)
      gp = 1.; // not convex here: plain gradient step
    const double step = std::clamp(g / gp, -maxStep, maxStep);
    s -= step;
    if (std::abs(step) < 1e-10 * (1. + std::abs(s)))
      break;
  }
  propagateRaw(s, pos, dir, nullptr);
  return s;
}

bool KinFitTrajectory::trackParameters(double &invP, double &theta,
                                       double &phi, double &R,
                                       double &Z) const {
  // Newton on h(s) = r_T(s) . t_T(s), h' = |t_T|^2 + r_T . (dt/ds)_T
  Vec pos, dir, dt;
  double s = 0;
  const double tT2 = fDir[0] * fDir[0] + fDir[1] * fDir[1];
  if (tT2 > 0)
    s = -(fPos[0] * fDir[0] + fPos[1] * fDir[1]) / tT2;
  bool ok = tT2 == 0;
  for (int iter = 0; iter < 50 && !ok; ++iter) {
    propagateRaw(s, pos, dir, &dt);
    const double h = pos[0] * dir[0] + pos[1] * dir[1];
    const double dT2 = dir[0] * dir[0] + dir[1] * dir[1];
    double hp = dT2 + pos[0] * dt[0] + pos[1] * dt[1];
    if (!(hp > 0))
      hp = dT2;
    if (!(hp > 0))
      break;
    const double step = h / hp;
    s -= step;
    ok = std::abs(step) < 1e-10 * (1. + std::abs(s));
  }
  propagateRaw(s, pos, dir, nullptr);
  invP = 1. / fP;
  theta = std::atan2(std::sqrt(dir[0] * dir[0] + dir[1] * dir[1]), dir[2]);
  phi = std::atan2(dir[1], dir[0]);
  R = -pos[0] * std::sin(phi) + pos[1] * std::cos(phi);
  Z = pos[2];
  return ok;
}
