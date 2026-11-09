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
 * toy_lambda_fit.cxx
 *
 * Self-contained example, no input file needed:
 *   - generates Lambda -> p pi- decays at a displaced vertex (no B field,
 *     straight tracks),
 *   - converts the tracks to the KinFit parametrisation (1/p, theta, phi, R, Z),
 *   - smears them according to a diagonal covariance,
 *   - runs a mass-constraint fit and a vertex-constraint fit,
 *   - reports chi2/probability and pull widths.
 *
 * Units: MeV/c, MeV/c^2, mm, rad.
 *
 * Build: cmake -DKINFIT_BUILD_EXAMPLES=ON .. && make
 *   or : root -l -e 'gSystem->Load("libKinFit")' toy_lambda_fit.cxx
 */

// In an interpreted ROOT session the KinFit classes are autoloaded from the
// library dictionary; including the headers again would clash with it.
#ifndef __CLING__
#include "KFitVertexFinder.h"
#include "KinFitParticle.h"
#include "KinFitter.h"
#endif

#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

constexpr double kMassLambda = 1115.683; // MeV/c^2
constexpr double kMassProton = 938.272;  // MeV/c^2
constexpr double kMassPion = 139.570;    // MeV/c^2

// Track resolution: sigma(1/p) is relative, the others absolute.
constexpr double kRelSigmaInvP = 0.02; // 2 % on 1/p
constexpr double kSigmaTheta = 0.002;  // rad
constexpr double kSigmaPhi = 0.002;    // rad
constexpr double kSigmaR = 0.5;        // mm
constexpr double kSigmaZ = 1.0;        // mm

/** Straight track from vertex `v` with momentum `p`: compute (R, Z) at the
 *  point of closest approach to the z axis, using the KinFit convention
 *  PCA = (-R sin(phi), R cos(phi), Z). */
void straightTrackRZ(const TVector3 &v, const TVector3 &p, double &R,
                     double &Z) {
  const double t = -(v.X() * p.X() + v.Y() * p.Y()) / p.Perp2();
  const TVector3 pca = v + t * p;
  R = -pca.X() * std::sin(p.Phi()) + pca.Y() * std::cos(p.Phi());
  Z = pca.Z();
}

/** Build a smeared KinFitParticle with diagonal covariance. */
KinFitParticle makeMeasuredTrack(const TLorentzVector &truth,
                                 const TVector3 &vertex, TRandom3 &rnd) {
  double R, Z;
  straightTrackRZ(vertex, truth.Vect(), R, Z);

  const double invP = 1.0 / truth.P();
  const double sInvP = kRelSigmaInvP * invP;

  const double mInvP = rnd.Gaus(invP, sInvP);
  const double mTheta = rnd.Gaus(truth.Theta(), kSigmaTheta);
  const double mPhi = rnd.Gaus(truth.Phi(), kSigmaPhi);
  const double mR = rnd.Gaus(R, kSigmaR);
  const double mZ = rnd.Gaus(Z, kSigmaZ);

  TLorentzVector lv;
  const double p = 1.0 / mInvP;
  lv.SetXYZM(p * std::sin(mTheta) * std::cos(mPhi),
             p * std::sin(mTheta) * std::sin(mPhi), p * std::cos(mTheta),
             truth.M());

  KinFitParticle cand(lv, mR, mZ);
  TMatrixD cov(5, 5);
  cov(0, 0) = sInvP * sInvP;
  cov(1, 1) = kSigmaTheta * kSigmaTheta;
  cov(2, 2) = kSigmaPhi * kSigmaPhi;
  cov(3, 3) = kSigmaR * kSigmaR;
  cov(4, 4) = kSigmaZ * kSigmaZ;
  cand.setCovariance(cov);
  return cand;
}

/** Generate one Lambda -> p pi- decay; returns daughters and decay vertex. */
void generateLambda(TRandom3 &rnd, TLorentzVector &proton,
                    TLorentzVector &pion, TVector3 &vertex) {
  // Lambda: p = 800..1500 MeV/c, theta = 20..60 deg, flight 20..80 mm.
  TVector3 pL;
  pL.SetMagThetaPhi(rnd.Uniform(800, 1500),
                    rnd.Uniform(20, 60) * TMath::DegToRad(),
                    rnd.Uniform(-TMath::Pi(), TMath::Pi()));
  TLorentzVector lambda;
  lambda.SetVectM(pL, kMassLambda);
  vertex = pL.Unit() * rnd.Uniform(20, 80);

  // Two-body decay, isotropic in the Lambda rest frame.
  const double pStar =
      std::sqrt((kMassLambda * kMassLambda -
                 std::pow(kMassProton + kMassPion, 2)) *
                (kMassLambda * kMassLambda -
                 std::pow(kMassProton - kMassPion, 2))) /
      (2 * kMassLambda);
  TVector3 dir;
  dir.SetMagThetaPhi(1.0, std::acos(rnd.Uniform(-1, 1)),
                     rnd.Uniform(-TMath::Pi(), TMath::Pi()));
  proton.SetVectM(pStar * dir, kMassProton);
  pion.SetVectM(-pStar * dir, kMassPion);
  proton.Boost(lambda.BoostVector());
  pion.Boost(lambda.BoostVector());
}

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
};

} // namespace

int toy_lambda_fit(int nEvents = 2000) {
  TRandom3 rnd(4711);

  Stats massBefore, massAfter, massProb, vtxProb, vtxDist;
  Stats massPull[5], vtxPull[5]; // pulls of the proton parameters
  int massOk = 0, vtxOk = 0;

  for (int ev = 0; ev < nEvents; ++ev) {
    TLorentzVector pTrue, piTrue;
    TVector3 vTrue;
    generateLambda(rnd, pTrue, piTrue, vTrue);

    std::vector<KinFitParticle> cands{makeMeasuredTrack(pTrue, vTrue, rnd),
                                      makeMeasuredTrack(piTrue, vTrue, rnd)};
    massBefore.add((cands[0] + cands[1]).M());

    // --- Mass constraint: M(p pi-) = M(Lambda) -----------------------------
    KinFitter massFit(cands);
    massFit.addMassConstraint(kMassLambda);
    if (massFit.fit()) {
      ++massOk;
      massAfter.add((massFit.getDaughter(0) + massFit.getDaughter(1)).M());
      massProb.add(massFit.getProb());
      for (int k = 0; k < 5; ++k)
        massPull[k].add(massFit.getPull(k));
    }

    // --- Vertex constraint: both tracks from one point ---------------------
    KinFitter vtxFit(cands);
    vtxFit.addVertexConstraint();
    if (vtxFit.fit()) {
      ++vtxOk;
      vtxProb.add(vtxFit.getProb());
      for (int k = 0; k < 5; ++k)
        vtxPull[k].add(vtxFit.getPull(k));

      // Vertex position from the fitted tracks.
      std::vector<KinFitParticle> fitted;
      vtxFit.getDaughters(fitted);
      KFitVertexFinder finder(fitted);
      vtxDist.add((finder.getVertex() - vTrue).Mag());
    }
  }

  const char *names[5] = {"1/p", "theta", "phi", "R", "Z"};

  std::printf("Lambda -> p pi-  toy, %d events\n\n", nEvents);
  std::printf("Mass fit:   converged %d/%d\n", massOk, nEvents);
  std::printf("  M(p pi-) before: mean %8.2f  rms %6.2f MeV/c^2\n",
              massBefore.mean(), massBefore.rms());
  std::printf("  M(p pi-) after : mean %8.2f  rms %6.2f MeV/c^2\n",
              massAfter.mean(), massAfter.rms());
  std::printf("  <prob> = %.3f (expect 0.5 for a flat distribution)\n",
              massProb.mean());
  std::printf("  proton pulls (mean / width):");
  for (int k = 0; k < 3; ++k)
    std::printf("  %s %+.2f/%.2f", names[k], massPull[k].mean(),
                massPull[k].rms());
  std::printf("\n\n");

  std::printf("Vertex fit: converged %d/%d\n", vtxOk, nEvents);
  std::printf("  <prob> = %.3f\n", vtxProb.mean());
  // 1/p does not enter the vertex constraint: its pull is undefined
  // (KinFitter reports -10000), so start at theta.
  std::printf("  proton pulls (mean / width):");
  for (int k = 1; k < 5; ++k)
    std::printf("  %s %+.2f/%.2f", names[k], vtxPull[k].mean(),
                vtxPull[k].rms());
  std::printf("\n  |vertex - truth| from fitted tracks: mean %.2f mm\n",
              vtxDist.mean());
  return 0;
}

#ifndef __CLING__
int main(int argc, char **argv) {
  return toy_lambda_fit(argc > 1 ? std::atoi(argv[1]) : 2000);
}
#endif
