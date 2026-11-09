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

/**
 * ToyTracks.h
 *
 * Toy Monte Carlo helpers for the QA suite: tracks with known truth in the
 * KinFit parametrisation (1/p, theta, phi, R, Z), Gaussian smearing, a
 * Lambda -> p pi- generator and simple statistics.
 *
 * Units: MeV/c, MeV/c^2, mm, rad, T.
 */

#pragma once

#include "KinFitParticle.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace toy {

constexpr double kMassLambda = 1115.683;
constexpr double kMassProton = 938.272;
constexpr double kMassPion = 139.570;
constexpr double kMassKaon = 493.677;

/** Track resolution. sigma(1/p) is relative to 1/p, the others absolute. */
struct Resolution {
  double relInvP = 0.01;
  double theta = 0.002;
  double phi = 0.002;
  double R = 0.5;
  double Z = 1.0;
};

/** (R, Z) of the straight track through `v` with momentum `p` at its point
 *  of closest approach (PCA) to the z axis; PCA = (-R sin phi, R cos phi, Z). */
inline void straightRZ(const TVector3 &v, const TVector3 &p, double &R,
                       double &Z) {
  const double t = -(v.X() * p.X() + v.Y() * p.Y()) / p.Perp2();
  const TVector3 pca = v + t * p;
  R = -pca.X() * std::sin(p.Phi()) + pca.Y() * std::cos(p.Phi());
  Z = pca.Z();
}

/** Diagonal 5x5 covariance for a track with momentum magnitude p. */
inline TMatrixD covariance(double p, const Resolution &res) {
  TMatrixD cov(5, 5);
  const double sInvP = res.relInvP / p;
  cov(0, 0) = sInvP * sInvP;
  cov(1, 1) = res.theta * res.theta;
  cov(2, 2) = res.phi * res.phi;
  cov(3, 3) = res.R * res.R;
  cov(4, 4) = res.Z * res.Z;
  return cov;
}

/** Exact (unsmeared) straight track from vertex `v`. */
inline KinFitParticle exactTrack(const TLorentzVector &truth,
                                 const TVector3 &v, const Resolution &res) {
  double R, Z;
  straightRZ(v, truth.Vect(), R, Z);
  KinFitParticle cand(truth, R, Z);
  cand.setCovariance(covariance(truth.P(), res));
  return cand;
}

/** Straight track from vertex `v`, smeared in (1/p, theta, phi, R, Z)
 *  according to `res`. The covariance is that of the truth momentum. */
inline KinFitParticle smearedTrack(const TLorentzVector &truth,
                                   const TVector3 &v, const Resolution &res,
                                   TRandom3 &rnd) {
  double R, Z;
  straightRZ(v, truth.Vect(), R, Z);
  const double invP = 1.0 / truth.P();
  const double mInvP = rnd.Gaus(invP, res.relInvP * invP);
  const double mTheta = rnd.Gaus(truth.Theta(), res.theta);
  const double mPhi = rnd.Gaus(truth.Phi(), res.phi);
  const double p = 1.0 / mInvP;
  TLorentzVector lv;
  lv.SetXYZM(p * std::sin(mTheta) * std::cos(mPhi),
             p * std::sin(mTheta) * std::sin(mPhi), p * std::cos(mTheta),
             truth.M());
  KinFitParticle cand(lv, rnd.Gaus(R, res.R), rnd.Gaus(Z, res.Z));
  cand.setCovariance(covariance(truth.P(), res));
  return cand;
}

/** Two-body decay of `mother` into masses m1, m2, isotropic in the rest
 *  frame. */
inline void twoBodyDecay(const TLorentzVector &mother, double m1, double m2,
                         TRandom3 &rnd, TLorentzVector &d1,
                         TLorentzVector &d2) {
  const double M = mother.M();
  const double pStar =
      std::sqrt((M * M - (m1 + m2) * (m1 + m2)) * (M * M - (m1 - m2) * (m1 - m2))) /
      (2 * M);
  TVector3 dir;
  dir.SetMagThetaPhi(1.0, std::acos(rnd.Uniform(-1, 1)),
                     rnd.Uniform(-TMath::Pi(), TMath::Pi()));
  d1.SetVectM(pStar * dir, m1);
  d2.SetVectM(-pStar * dir, m2);
  d1.Boost(mother.BoostVector());
  d2.Boost(mother.BoostVector());
}

/** A Lambda with p = 800..1500 MeV/c, theta = 20..60 deg, produced at a
 *  primary vertex near the origin (sigma_xy = 1 mm, |z| < 5 mm) and
 *  decaying after 20..80 mm. */
struct LambdaEvent {
  TLorentzVector lambda, proton, pion;
  TVector3 primaryVertex, decayVertex;
};

inline LambdaEvent lambdaEvent(TRandom3 &rnd) {
  LambdaEvent ev;
  TVector3 pL;
  pL.SetMagThetaPhi(rnd.Uniform(800, 1500),
                    rnd.Uniform(20, 60) * TMath::DegToRad(),
                    rnd.Uniform(-TMath::Pi(), TMath::Pi()));
  ev.lambda.SetVectM(pL, kMassLambda);
  ev.primaryVertex.SetXYZ(rnd.Gaus(0, 1), rnd.Gaus(0, 1), rnd.Uniform(-5, 5));
  ev.decayVertex = ev.primaryVertex + pL.Unit() * rnd.Uniform(20, 80);
  twoBodyDecay(ev.lambda, kMassProton, kMassPion, rnd, ev.proton, ev.pion);
  return ev;
}

/** Helix track in constant Bz [T] from vertex `v`, returned at its PCA to
 *  the z axis (theta, phi, p, R, Z at the PCA). Positive charges turn
 *  clockwise for Bz > 0. */
inline KinFitParticle helixTrackAtPCA(const TVector3 &v, const TVector3 &p,
                                      int q, double bz, double mass) {
  const double pt = p.Perp();
  const double rho = pt / (0.299792458 * std::abs(q * bz));
  const double h = (q * bz > 0) ? -1.0 : 1.0;
  const double phi0 = p.Phi();
  const TVector2 c(v.X() - h * rho * std::sin(phi0),
                   v.Y() + h * rho * std::cos(phi0));
  const TVector2 pca = c - rho * c.Unit();
  const double aV = std::atan2(v.Y() - c.Y(), v.X() - c.X());
  const double aP = std::atan2(pca.Y() - c.Y(), pca.X() - c.X());
  const double dA = TVector2::Phi_mpi_pi(aP - aV);
  const double s = h * dA * rho;
  const double phiPCA = phi0 + dA;
  const double zPCA = v.Z() + s * p.Z() / pt;
  const double R = -pca.X() * std::sin(phiPCA) + pca.Y() * std::cos(phiPCA);
  TLorentzVector lv;
  lv.SetXYZM(pt * std::cos(phiPCA), pt * std::sin(phiPCA), p.Z(), mass);
  KinFitParticle cand(lv, R, zPCA);
  cand.setCharge(q);
  return cand;
}

/** Parameter vector (1/p, theta, phi, R, Z) of a set of particles, laid out
 *  as KinFitter stores it internally. */
inline TMatrixD parameterVector(const std::vector<KinFitParticle> &cands) {
  TMatrixD y(5 * cands.size(), 1);
  for (size_t i = 0; i < cands.size(); ++i) {
    y(5 * i + 0, 0) = 1.0 / cands[i].P();
    y(5 * i + 1, 0) = cands[i].Theta();
    y(5 * i + 2, 0) = cands[i].Phi();
    y(5 * i + 3, 0) = cands[i].getR();
    y(5 * i + 4, 0) = cands[i].getZ();
  }
  return y;
}

/** Running mean / RMS. */
struct Stats {
  int n = 0;
  double sum = 0, sum2 = 0;
  void add(double x) {
    ++n;
    sum += x;
    sum2 += x * x;
  }
  double mean() const { return n ? sum / n : 0; }
  double rms() const {
    return n ? std::sqrt(std::max(0.0, sum2 / n - mean() * mean())) : 0;
  }
  /** Standard error of the mean. */
  double meanError() const { return n > 1 ? rms() / std::sqrt(n) : 0; }
};

} // namespace toy
