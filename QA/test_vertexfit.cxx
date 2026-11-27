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

// Tests for the common-vertex constraint of KinFitter: all charged tracks,
// track model in every field type (B = 0, Bz, By, By(z), By(x, z), user
// function), neutral particles measured by energy and hit position.
//
//   exact_*    : unsmeared input -> chi2 ~ 0, vertex = truth
//   pulls_*    : smeared input   -> flat probability, pulls ~ N(0, 1) for
//                track parameters and the fitted vertex
//   jacobian_* : Feta_eval / Fxi_eval == finite differences of f_eval

#include "JacobianCheck.h"
#include "KinFitTest.h"
#include "ToyTracks.h"

#include "KinFitFieldByXZ.h"
#include "KinFitFieldByZ.h"
#include "KinFitFieldFunction.h"
#include "KinFitFieldUniform.h"
#include "KinFitter.h"

#include <TGenPhaseSpace.h>
#include <TH1D.h>
#include <TH2D.h>

#include <memory>

namespace {

enum class F { Zero, Bz, By, ByZ, ByXZ, Function };

const char *name(F f) {
  switch (f) {
  case F::Zero: return "B=0";
  case F::Bz: return "Bz";
  case F::By: return "By";
  case F::ByZ: return "By(z)";
  case F::ByXZ: return "By(x,z)";
  case F::Function: return "function";
  }
  return "?";
}

/** Dipole-like profile: 0.8 T peak at z = 600 mm, width 400 mm. */
double profile(double z) { return 0.8 * std::exp(-0.5 * std::pow((z - 600.) / 400., 2)); }

std::unique_ptr<KinFitField> makeField(F f) {
  switch (f) {
  case F::Zero:
    return std::make_unique<KinFitFieldUniform>();
  case F::Bz:
    return std::make_unique<KinFitFieldUniform>(KinFitFieldUniform::Bz(1.0));
  case F::By:
    return std::make_unique<KinFitFieldUniform>(KinFitFieldUniform::By(0.7));
  case F::ByZ: {
    TH1D h("hByZ", "", 120, -1000., 3000.);
    for (int i = 1; i <= 120; ++i)
      h.SetBinContent(i, profile(h.GetBinCenter(i)));
    return std::make_unique<KinFitFieldByZ>(h);
  }
  case F::ByXZ: {
    TH2D h("hByXZ", "", 40, -2000., 2000., 120, -1000., 3000.);
    for (int ix = 1; ix <= 40; ++ix)
      for (int iz = 1; iz <= 120; ++iz) {
        const double x = h.GetXaxis()->GetBinCenter(ix);
        h.SetBinContent(ix, iz,
                        profile(h.GetYaxis()->GetBinCenter(iz)) *
                            (1. - 1e-7 * x * x));
      }
    return std::make_unique<KinFitFieldByXZ>(h);
  }
  case F::Function:
    return std::make_unique<KinFitFieldFunction>(
        [](double x, double, double z) {
          return TVector3(0.05, profile(z) * (1. + 1e-4 * x), 0.3);
        });
  }
  return nullptr;
}

const toy::Resolution kRes;
const toy::NeutralResolution kNeutralRes;

TVector3 randomVertex(TRandom3 &rnd) {
  return TVector3(rnd.Gaus(0, 5), rnd.Gaus(0, 5), rnd.Uniform(0, 60));
}

/** 4 charged tracks (charges alternating) from vertex v. */
std::vector<KinFitParticle> chargedTracks(const TVector3 &v,
                                          const KinFitField *field,
                                          TRandom3 &rnd, bool smear) {
  std::vector<KinFitParticle> cands;
  for (int k = 0; k < 4; ++k) {
    TVector3 p;
    p.SetMagThetaPhi(rnd.Uniform(300, 1500), rnd.Uniform(0.15, 0.6),
                     rnd.Uniform(-TMath::Pi(), TMath::Pi()));
    TLorentzVector lv;
    lv.SetVectM(p, toy::kMassPion);
    cands.push_back(toy::fieldTrack(lv, v, k % 2 ? 1 : -1, field, kRes,
                                    smear ? &rnd : nullptr));
  }
  return cands;
}

/** Mother (1300 MeV, known four-momentum) at vertex v -> p pi- gamma;
 *  the photon is measured by energy and hit 1500 mm away. */
struct NeutralEvent {
  std::vector<KinFitParticle> cands; // p, pi-, gamma
  TLorentzVector mother;
  TVector3 vertex;
};

NeutralEvent neutralEvent(const KinFitField *field, TRandom3 &rnd,
                          bool smear) {
  NeutralEvent ev;
  ev.vertex = randomVertex(rnd);
  TVector3 pm;
  pm.SetMagThetaPhi(rnd.Uniform(800, 1500), rnd.Uniform(0.1, 0.3),
                    rnd.Uniform(-TMath::Pi(), TMath::Pi()));
  ev.mother.SetVectM(pm, 1300.);
  TGenPhaseSpace ps;
  const double masses[3] = {toy::kMassProton, toy::kMassPion, 0.};
  ps.SetDecay(ev.mother, 3, masses);
  ps.Generate(); // unweighted is fine for a toy
  TRandom3 *r = smear ? &rnd : nullptr;
  ev.cands.push_back(toy::fieldTrack(*ps.GetDecay(0), ev.vertex, +1, field,
                                     kRes, r));
  ev.cands.push_back(toy::fieldTrack(*ps.GetDecay(1), ev.vertex, -1, field,
                                     kRes, r));
  ev.cands.push_back(
      toy::neutralHit(*ps.GetDecay(2), ev.vertex, 1500., kNeutralRes, r));
  return ev;
}

void runExactCharged(F fld) {
  const auto field = makeField(fld);
  TRandom3 rnd(10);
  double worstChi2 = 0, worstV = 0;
  for (int ev = 0; ev < 20; ++ev) {
    const TVector3 v = randomVertex(rnd);
    KinFitter f(chargedTracks(v, field.get(), rnd, false));
    f.setField(*field);
    f.addCommonVertexConstraint();
    KF_CHECK(f.fit());
    KF_CHECK(f.getNDF() == 5);
    worstChi2 = std::max(worstChi2, f.getChi2());
    worstV = std::max(worstV, (f.getVertex() - v).Mag());
  }
  kftest::info("%s exact: worst chi2 %g, worst |v - truth| %g mm", name(fld),
               worstChi2, worstV);
  KF_CHECK(worstChi2 < 1e-6);
  KF_CHECK(worstV < 1e-4);
}

/** Pull statistics; pull values -10000 (undefined) are skipped. */
struct PullSet {
  std::vector<toy::Stats> s;
  explicit PullSet(int n) : s(n) {}
  void add(int i, double v) {
    if (v > -9999)
      s[i].add(v);
  }
};

void checkPull(const char *what, const toy::Stats &st) {
  if (st.n < 50)
    return;
  kftest::info("  %-12s mean %+.3f width %.3f", what, st.mean(), st.rms());
  KF_CHECK_NEAR(st.mean(), 0, 0.15);
  KF_CHECK_NEAR(st.rms(), 1, 0.12);
}

void checkProb(const char *what, int nConv, int nEvents,
               const toy::Stats &prob) {
  kftest::info("%s: converged %d/%d, <prob> = %.3f +- %.3f", what, nConv,
               nEvents, prob.mean(), prob.meanError());
  KF_CHECK(nConv >= 0.98 * nEvents);
  KF_CHECK_NEAR(prob.mean(), 0.5, 0.05);
}

void runPullsCharged(F fld) {
  const auto field = makeField(fld);
  TRandom3 rnd(20);
  const int nEvents = 400;
  int nConv = 0;
  toy::Stats prob;
  PullSet par(20), vtx(3);
  for (int ev = 0; ev < nEvents; ++ev) {
    const TVector3 v = randomVertex(rnd);
    KinFitter f(chargedTracks(v, field.get(), rnd, true));
    f.setField(*field);
    f.addCommonVertexConstraint();
    if (!f.fit())
      continue;
    ++nConv;
    prob.add(f.getProb());
    for (int i = 0; i < 20; ++i)
      par.add(i, f.getPull(i));
    const TMatrixD &c = f.getVertexCovariance();
    for (int j = 0; j < 3; ++j)
      vtx.add(j, (f.getVertex()[j] - v[j]) / std::sqrt(c(j, j)));
  }
  checkProb(name(fld), nConv, nEvents, prob);
  const char *pn[5] = {"1/p", "theta", "phi", "R", "Z"};
  for (int i = 0; i < 20; ++i)
    checkPull((std::string("trk") + std::to_string(i / 5) + " " + pn[i % 5])
                  .c_str(),
              par.s[i]);
  const char *vn[3] = {"vertex x", "vertex y", "vertex z"};
  for (int j = 0; j < 3; ++j)
    checkPull(vn[j], vtx.s[j]);
}

void runJacobianCharged(F fld) {
  const auto field = makeField(fld);
  TRandom3 rnd(30);
  for (int ev = 0; ev < 3; ++ev) {
    KinFitter f(chargedTracks(randomVertex(rnd), field.get(), rnd, true));
    f.setField(*field);
    f.addCommonVertexConstraint();
    TMatrixD y, xi;
    KF_CHECK(f.getStartValues(y, xi));
    compareJacobian("Feta", f.Feta_eval(y, xi),
                    [&](const TMatrixD &yy) { return f.f_eval(yy, xi); }, y,
                    1e-2);
    compareJacobian("Fxi", f.Fxi_eval(y, xi),
                    [&](const TMatrixD &xx) { return f.f_eval(y, xx); }, xi,
                    1.);
  }
}

void runNeutral4C(F fld) {
  const auto field = makeField(fld);
  TRandom3 rnd(40);
  const int nEvents = 400;
  int nConv = 0;
  toy::Stats prob;
  PullSet par(15), vtx(3);
  for (int ev = 0; ev < nEvents; ++ev) {
    const NeutralEvent e = neutralEvent(field.get(), rnd, true);
    KinFitter f(e.cands);
    f.setField(*field);
    f.add4Constraint(e.mother);
    f.addCommonVertexConstraint();
    if (!f.fit())
      continue;
    ++nConv;
    KF_CHECK(f.getNDF() == 5);
    prob.add(f.getProb());
    for (int i = 0; i < 15; ++i)
      par.add(i, f.getPull(i));
    const TMatrixD &c = f.getVertexCovariance();
    for (int j = 0; j < 3; ++j)
      vtx.add(j, (f.getVertex()[j] - e.vertex[j]) / std::sqrt(c(j, j)));
    if (nConv <= 10) { // constraint fulfilled, photon points from the vertex
      const TLorentzVector sum =
          f.getDaughter(0) + f.getDaughter(1) + f.getDaughter(2);
      KF_CHECK_NEAR((sum - e.mother).Vect().Mag(), 0, 1e-2);
      KF_CHECK_NEAR(sum.E() - e.mother.E(), 0, 1e-2);
      const KinFitParticle &g = f.getDaughter(2);
      KF_CHECK_NEAR(g.Vect().Unit().Dot((g.getHit() - f.getVertex()).Unit()),
                    1, 1e-12);
    }
  }
  checkProb((std::string("p pi- gamma 4C+vertex ") + name(fld)).c_str(), nConv,
            nEvents, prob);
  // nearly collinear p pi- pairs used to fail or end in a wrong minimum
  // (second start value in fit())
  KF_CHECK(nConv == nEvents);
  const char *gn[4] = {"gamma E", "gamma X", "gamma Y", "gamma Z"};
  for (int i = 0; i < 4; ++i)
    checkPull(gn[i], par.s[10 + i]);
  const char *pn[5] = {"1/p", "theta", "phi", "R", "Z"};
  for (int i = 0; i < 10; ++i)
    checkPull((std::string("trk") + std::to_string(i / 5) + " " + pn[i % 5])
                  .c_str(),
              par.s[i]);
  const char *vn[3] = {"vertex x", "vertex y", "vertex z"};
  for (int j = 0; j < 3; ++j)
    checkPull(vn[j], vtx.s[j]);
}

} // namespace

KF_TEST(exact_B0) { runExactCharged(F::Zero); }
KF_TEST(exact_Bz) { runExactCharged(F::Bz); }
KF_TEST(exact_By) { runExactCharged(F::By); }
KF_TEST(exact_ByZ) { runExactCharged(F::ByZ); }
KF_TEST(exact_ByXZ) { runExactCharged(F::ByXZ); }
KF_TEST(exact_function) { runExactCharged(F::Function); }

KF_TEST(pulls_B0) { runPullsCharged(F::Zero); }
KF_TEST(pulls_Bz) { runPullsCharged(F::Bz); }
KF_TEST(pulls_By) { runPullsCharged(F::By); }
KF_TEST(pulls_ByZ) { runPullsCharged(F::ByZ); }
KF_TEST(pulls_ByXZ) { runPullsCharged(F::ByXZ); }
KF_TEST(pulls_function) { runPullsCharged(F::Function); }

KF_TEST(jacobian_B0) { runJacobianCharged(F::Zero); }
KF_TEST(jacobian_Bz) { runJacobianCharged(F::Bz); }
KF_TEST(jacobian_ByZ) { runJacobianCharged(F::ByZ); }
KF_TEST(jacobian_function) { runJacobianCharged(F::Function); }

// Neutral photon + two tracks, 4C with known mother + common vertex.
KF_TEST(neutral_4C_vertex_Bz) { runNeutral4C(F::Bz); }
KF_TEST(neutral_4C_vertex_ByZ) { runNeutral4C(F::ByZ); }

// Jacobians with a neutral hit: its momentum depends on the vertex (Fxi).
KF_TEST(jacobian_neutral_4C_vertex) {
  const auto field = makeField(F::ByZ);
  TRandom3 rnd(50);
  for (int ev = 0; ev < 3; ++ev) {
    const NeutralEvent e = neutralEvent(field.get(), rnd, true);
    KinFitter f(e.cands);
    f.setField(*field);
    f.add4Constraint(e.mother);
    f.addCommonVertexConstraint();
    f.addMassConstraint(1000.); // mass row with the neutral too
    TMatrixD y, xi;
    KF_CHECK(f.getStartValues(y, xi));
    // larger steps: rows of O(1e6 MeV^2) need h >> rounding noise
    compareJacobian("Feta", f.Feta_eval(y, xi),
                    [&](const TMatrixD &yy) { return f.f_eval(yy, xi); }, y,
                    1e-2);
    compareJacobian("Fxi", f.Fxi_eval(y, xi),
                    [&](const TMatrixD &xx) { return f.f_eval(y, xx); }, xi,
                    100.);
  }
}

// pi0 -> gamma gamma from a known production point (no vertex fit): the
// photon directions are hit - origin; mass constraint.
KF_TEST(neutral_mass_fixed_origin) {
  TRandom3 rnd(60);
  const double mPi0 = 134.977;
  const int nEvents = 1000;
  int nConv = 0;
  toy::Stats prob;
  PullSet par(10);
  for (int ev = 0; ev < nEvents; ++ev) {
    const TVector3 v = randomVertex(rnd);
    TVector3 p;
    p.SetMagThetaPhi(rnd.Uniform(300, 1500), rnd.Uniform(0.1, 0.5),
                     rnd.Uniform(-TMath::Pi(), TMath::Pi()));
    TLorentzVector pi0, g1, g2;
    pi0.SetVectM(p, mPi0);
    toy::twoBodyDecay(pi0, 0., 0., rnd, g1, g2);
    std::vector<KinFitParticle> cands = {
        toy::neutralHit(g1, v, 1500., kNeutralRes, &rnd, v),
        toy::neutralHit(g2, v, 1500., kNeutralRes, &rnd, v)};
    KinFitter f(cands);
    f.addMassConstraint(mPi0);
    if (!f.fit())
      continue;
    ++nConv;
    KF_CHECK(f.getNDF() == 1);
    prob.add(f.getProb());
    for (int i = 0; i < 10; ++i)
      par.add(i, f.getPull(i));
    if (nConv <= 10)
      KF_CHECK_NEAR((f.getDaughter(0) + f.getDaughter(1)).M(), mPi0, 1e-3);
  }
  checkProb("pi0 -> gamma gamma", nConv, nEvents, prob);
  const char *gn[4] = {"E", "X", "Y", "Z"};
  for (int i = 0; i < 10; ++i)
    if (i % 5 < 4)
      checkPull((std::string("gamma") + std::to_string(i / 5) + " " +
                 gn[i % 5]).c_str(),
                par.s[i]);
}

// Lambda -> p pi- in a By(z) field: invariant mass at the fitted vertex
// (momenta transported from the reference points to the vertex).
KF_TEST(lambda_mass_vertex_ByZ) {
  const auto field = makeField(F::ByZ);
  TRandom3 rnd(80);
  const int nEvents = 400;
  int nConv = 0;
  toy::Stats prob;
  PullSet vtx(3);
  for (int ev = 0; ev < nEvents; ++ev) {
    const toy::LambdaEvent l = toy::lambdaEvent(rnd);
    std::vector<KinFitParticle> cands = {
        toy::fieldTrack(l.proton, l.decayVertex, +1, field.get(), kRes, &rnd),
        toy::fieldTrack(l.pion, l.decayVertex, -1, field.get(), kRes, &rnd)};
    KinFitter f(cands);
    f.setField(*field);
    f.addMassVtxConstraint(toy::kMassLambda);
    if (!f.fit())
      continue;
    ++nConv;
    KF_CHECK(f.getNDF() == 2);
    prob.add(f.getProb());
    const TMatrixD &c = f.getVertexCovariance();
    for (int j = 0; j < 3; ++j)
      vtx.add(j, (f.getVertex()[j] - l.decayVertex[j]) / std::sqrt(c(j, j)));
    if (nConv <= 10) // daughters (at the vertex) give the Lambda mass
      KF_CHECK_NEAR((f.getDaughter(0) + f.getDaughter(1)).M(),
                    toy::kMassLambda, 1e-2);
  }
  checkProb("Lambda mass+vertex By(z)", nConv, nEvents, prob);
  const char *vn[3] = {"vertex x", "vertex y", "vertex z"};
  for (int j = 0; j < 3; ++j)
    checkPull(vn[j], vtx.s[j]);
}

// The legacy vertex constraint now uses all charged tracks.
KF_TEST(legacy_vertex_uses_all_tracks) {
  const KinFitFieldUniform zero;
  TRandom3 rnd(70);
  const TVector3 v = randomVertex(rnd);
  KinFitter f(chargedTracks(v, &zero, rnd, true));
  f.addVertexConstraint();
  KF_CHECK(f.fit());
  KF_CHECK(f.getNDF() == 5); // 2 * 4 - 3 (was 1: first two tracks)
}

KF_TEST_MAIN()
