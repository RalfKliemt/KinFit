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
 * vertex_field_fit.cxx
 *
 * Self-contained example, no input file needed: common-vertex fit in a
 * dipole-like magnetic field with a neutral particle.
 *
 *   - field: By(z) from a TH1 (Gaussian profile, 0.8 T at z = 600 mm);
 *     any KinFitField works (KinFitFieldUniform::Bz/By, KinFitFieldByXZ,
 *     KinFitFieldFunction, ...)
 *   - event: a mother (1300 MeV/c^2, known four-momentum) decays at a
 *     displaced vertex into p pi- gamma
 *   - charged tracks: propagated in the field with KinFitTrajectory and
 *     stored at their point of closest approach to the z axis (the usual
 *     KinFit parametrisation), then smeared
 *   - photon: measured only by its energy and a hit 1.5 m away
 *     (KinFitParticle::fromNeutralHit); its direction is hit - vertex
 *   - fit: 4C (known mother) + common vertex; the vertex is fitted, the
 *     charged momenta are taken at the vertex
 *
 * Reports fit probability, vertex resolution (start value from
 * KinFitVertexFinderField vs. fit) and pulls.
 *
 * Units: MeV/c, MeV/c^2, mm, rad, Tesla.
 *
 * Build: cmake -DKINFIT_BUILD_EXAMPLES=ON .. && make
 *   or : root -l -e 'gSystem->Load("libKinFit")' vertex_field_fit.cxx
 */

// In an interpreted ROOT session the KinFit classes are autoloaded from the
// library dictionary; including the headers again would clash with it.
#ifndef __CLING__
#include "KinFitFieldByZ.h"
#include "KinFitParticle.h"
#include "KinFitTrajectory.h"
#include "KinFitVertexFinderField.h"
#include "KinFitter.h"
#endif

#include <TGenPhaseSpace.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

const double kMassProton = 938.272, kMassPion = 139.570, kMassMother = 1300.;

// Track resolution: sigma(1/p)/(1/p), theta, phi [rad], R, Z [mm]
const double kRes[5] = {0.01, 0.002, 0.002, 0.5, 1.0};
// Photon: sigma(E)/E, hit x/y and z [mm]
const double kResE = 0.03, kResHit = 3.0;

struct Stats {
  int n = 0;
  double sum = 0, sum2 = 0;
  void add(double x) { ++n, sum += x, sum2 += x * x; }
  double mean() const { return n ? sum / n : 0; }
  double rms() const {
    return n ? std::sqrt(std::max(0., sum2 / n - mean() * mean())) : 0;
  }
};

/** 68 % quantile of |x|: a resolution measure that is robust against the
 *  few events with a poorly determined vertex (nearly collinear tracks). */
double quantile68(std::vector<double> x) {
  if (x.empty())
    return 0;
  for (double &v : x)
    v = std::abs(v);
  std::sort(x.begin(), x.end());
  return x[x.size() * 68 / 100];
}

/** Charged track from vertex v in the field, at its transverse PCA to the z
 *  axis, smeared with kRes. */
KinFitParticle chargedTrack(const TLorentzVector &truth, const TVector3 &v,
                            int q, const KinFitField &field, TRandom3 &rnd) {
  double invP, theta, phi, R, Z;
  KinFitTrajectory(v, truth.Vect(), q, &field)
      .trackParameters(invP, theta, phi, R, Z);
  const double truthInvP = invP;
  invP = rnd.Gaus(invP, kRes[0] * invP);
  theta = rnd.Gaus(theta, kRes[1]);
  phi = rnd.Gaus(phi, kRes[2]);
  R = rnd.Gaus(R, kRes[3]);
  Z = rnd.Gaus(Z, kRes[4]);

  TLorentzVector lv;
  lv.SetXYZM(std::sin(theta) * std::cos(phi) / invP,
             std::sin(theta) * std::sin(phi) / invP, std::cos(theta) / invP,
             truth.M());
  KinFitParticle track(lv, R, Z);
  track.setCharge(q); // needed for the track model in the field
  TMatrixD cov(5, 5);
  cov(0, 0) = std::pow(kRes[0] * truthInvP, 2);
  for (int i = 1; i < 5; ++i)
    cov(i, i) = kRes[i] * kRes[i];
  track.setCovariance(cov);
  return track;
}

/** Photon from vertex v hitting a detector 1.5 m away (smeared). */
KinFitParticle photon(const TLorentzVector &truth, const TVector3 &v,
                      TRandom3 &rnd) {
  const TVector3 hitTrue = v + 1500. * truth.Vect().Unit();
  const TVector3 hit(rnd.Gaus(hitTrue.X(), kResHit),
                     rnd.Gaus(hitTrue.Y(), kResHit),
                     rnd.Gaus(hitTrue.Z(), kResHit));
  const double E = rnd.Gaus(truth.E(), kResE * truth.E());
  TMatrixD cov(4, 4); // (E, X, Y, Z)
  cov(0, 0) = std::pow(kResE * truth.E(), 2);
  for (int i = 1; i < 4; ++i)
    cov(i, i) = kResHit * kResHit;
  // origin: only used without a vertex fit
  return KinFitParticle::fromNeutralHit(E, 0., hit, cov, TVector3(0, 0, 0));
}

} // namespace

int vertex_field_fit(int nEvents = 500) {
  // Dipole-like field By(z)
  TH1D hBy("hBy", "B_{y}(z);z [mm];B_{y} [T]", 200, -1000., 3000.);
  for (int i = 1; i <= hBy.GetNbinsX(); ++i) {
    const double z = hBy.GetBinCenter(i);
    hBy.SetBinContent(i, 0.8 * std::exp(-0.5 * std::pow((z - 600.) / 400., 2)));
  }
  const KinFitFieldByZ field(hBy); // copies the histogram

  TRandom3 rnd(1);
  int nOk = 0;
  Stats prob, vtxPull[3], pullE, pullTheta, massCheck;
  std::vector<double> startRes[3], fitRes[3]; // (estimate - truth) per axis
  int ndf = 0;

  for (int ev = 0; ev < nEvents; ++ev) {
    // Mother at a displaced vertex, decaying to p pi- gamma
    const TVector3 vTrue(rnd.Gaus(0, 5), rnd.Gaus(0, 5), rnd.Uniform(0, 60));
    TVector3 pm;
    pm.SetMagThetaPhi(rnd.Uniform(800, 1500), rnd.Uniform(0.1, 0.3),
                      rnd.Uniform(-TMath::Pi(), TMath::Pi()));
    TLorentzVector mother;
    mother.SetVectM(pm, kMassMother);
    TGenPhaseSpace decay;
    const double masses[3] = {kMassProton, kMassPion, 0.};
    decay.SetDecay(mother, 3, masses);
    decay.Generate();

    const std::vector<KinFitParticle> cands = {
        chargedTrack(*decay.GetDecay(0), vTrue, +1, field, rnd),
        chargedTrack(*decay.GetDecay(1), vTrue, -1, field, rnd),
        photon(*decay.GetDecay(2), vTrue, rnd)};

    // Unweighted vertex of the charged tracks (also the fit's start value)
    const KinFitVertexFinderField finder(cands, &field);

    KinFitter fitter(cands);
    fitter.setField(field);
    fitter.add4Constraint(mother);
    fitter.addCommonVertexConstraint();
    if (!fitter.fit())
      continue;
    ++nOk;
    ndf = fitter.getNDF();
    prob.add(fitter.getProb());

    const TVector3 &v = fitter.getVertex();
    const TMatrixD &cv = fitter.getVertexCovariance();
    for (int j = 0; j < 3; ++j) {
      startRes[j].push_back(finder.getVertex()[j] - vTrue[j]);
      fitRes[j].push_back(v[j] - vTrue[j]);
      vtxPull[j].add((v[j] - vTrue[j]) / std::sqrt(cv(j, j)));
    }
    pullE.add(fitter.getPull(10));    // photon E (particle 2, slot 0)
    pullTheta.add(fitter.getPull(1)); // proton theta
    // daughters (charged momenta at the vertex) add up to the mother
    const TLorentzVector sum = fitter.getDaughter(0) + fitter.getDaughter(1) +
                               fitter.getDaughter(2);
    massCheck.add(sum.M() - kMassMother);
  }

  std::printf("p pi- gamma in By(z), 4C + common vertex, %d events\n\n",
              nEvents);
  std::printf("converged %d/%d, NDF = %d (4 + 2 tracks x 2 - 3)\n", nOk,
              nEvents, ndf);
  std::printf("<prob> = %.3f (expect 0.5 for a flat distribution)\n",
              prob.mean());
  std::printf("M(sum of fitted daughters) - M(mother): mean %.1e, rms %.1e "
              "MeV/c^2\n\n",
              massCheck.mean(), massCheck.rms());
  const char *xyz[3] = {"x", "y", "z"};
  std::printf("vertex resolution [mm, 68%% of |estimate - truth|]\n");
  std::printf("        start (finder)      fit     fit pull mean/width\n");
  for (int j = 0; j < 3; ++j)
    std::printf("  %s        %8.3f       %8.3f      %+.2f/%.2f\n", xyz[j],
                quantile68(startRes[j]), quantile68(fitRes[j]),
                vtxPull[j].mean(), vtxPull[j].rms());
  std::printf("\npulls (mean/width): photon E %+.2f/%.2f, proton theta "
              "%+.2f/%.2f\n",
              pullE.mean(), pullE.rms(), pullTheta.mean(), pullTheta.rms());
  return 0;
}

#ifndef __CLING__
int main(int argc, char **argv) {
  return vertex_field_fit(argc > 1 ? std::atoi(argv[1]) : 500);
}
#endif
