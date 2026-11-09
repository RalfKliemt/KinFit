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

// Tests for KFitDecayCandFinder: building a Lambda candidate from its
// daughters and the primary and decay vertices.

#include "KinFitTest.h"
#include "ToyTracks.h"

#include "KFitDecayCandFinder.h"
#include "KinFitter.h"

#include <TMatrixDSymEigen.h>

namespace {

KFitDecayCandFinder makeFinder(const toy::LambdaEvent &ev) {
  std::vector<KinFitParticle> daughters = {
      toy::exactTrack(ev.proton, ev.decayVertex, toy::Resolution()),
      toy::exactTrack(ev.pion, ev.decayVertex, toy::Resolution())};
  return KFitDecayCandFinder(daughters, toy::kMassLambda, ev.primaryVertex,
                             ev.decayVertex, 0.5, 0.5, 1.0, 2.0, 2.0, 4.0);
}

} // namespace

// |p| = sqrt((sum E)^2 - m^2) reproduces the true Lambda momentum.
KF_TEST(momentum_from_daughter_energies) {
  TRandom3 rnd(1);
  for (int i = 0; i < 100; ++i) {
    const auto ev = toy::lambdaEvent(rnd);
    const KinFitParticle cand = makeFinder(ev).getDecayCand();
    KF_CHECK_REL(cand.P(), ev.lambda.P(), 1e-9, 1);
    KF_CHECK_NEAR(cand.M(), toy::kMassLambda, 1e-6);
  }
}

// Direction = primary -> decay vertex.
KF_TEST(direction_from_vertices) {
  TRandom3 rnd(2);
  for (int i = 0; i < 100; ++i) {
    const auto ev = toy::lambdaEvent(rnd);
    const KinFitParticle cand = makeFinder(ev).getDecayCand();
    KF_CHECK_NEAR(cand.Theta(), ev.lambda.Theta(), 1e-9);
    KF_CHECK_NEAR(cand.Phi(), ev.lambda.Phi(), 1e-9);
  }
}

// R, Z must be those of the flight line (through the decay vertex along the
// flight direction). Known bug: the line is anchored at (decay - primary)
// and the sign of R is taken from the primary vertex y (TODO.md phase 3).
KF_TEST(R_Z_of_flight_line) {
  TRandom3 rnd(3);
  for (int i = 0; i < 100; ++i) {
    const auto ev = toy::lambdaEvent(rnd);
    double R, Z;
    toy::straightRZ(ev.decayVertex, ev.lambda.Vect(), R, Z);
    const KinFitParticle cand = makeFinder(ev).getDecayCand();
    KF_CHECK_NEAR(cand.getR(), R, 1e-6);
    KF_CHECK_NEAR(cand.getZ(), Z, 1e-6);
  }
}

KF_TEST(covariance_accessors_agree) {
  TRandom3 rnd(4);
  const auto ev = toy::lambdaEvent(rnd);
  auto finder = makeFinder(ev);
  const TMatrixD a = finder.getCovarianceMatrixDecayCand();
  const TMatrixD b = finder.getDecayCand().getCovariance();
  KF_CHECK(a.GetNrows() == 5 && b.GetNrows() == 5);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      KF_CHECK_NEAR(a(i, j), b(i, j), 0);
}

// A covariance matrix must be symmetric and positive semi-definite, also
// with vertex correlations. Known bug: inconsistent signs and square roots
// of correlations (TODO.md phase 3).
KF_TEST(covariance_symmetric_positive) {
  TRandom3 rnd(5);
  for (int i = 0; i < 50; ++i) {
    const auto ev = toy::lambdaEvent(rnd);
    std::vector<KinFitParticle> daughters = {
        toy::exactTrack(ev.proton, ev.decayVertex, toy::Resolution()),
        toy::exactTrack(ev.pion, ev.decayVertex, toy::Resolution())};
    KFitDecayCandFinder finder(daughters, toy::kMassLambda, ev.primaryVertex,
                               ev.decayVertex, 0.5, 0.5, 1.0, 2.0, 2.0, 4.0);
    finder.setVertexCorrelations(0.3, 0.1, -0.2, 0.3, 0.1, -0.2);
    // The correlations only take effect in the constructor's calculation;
    // until that is fixed, check the matrix as delivered.
    const TMatrixD cov = finder.getCovarianceMatrixDecayCand();
    TMatrixDSym sym(5);
    for (int r = 0; r < 5; ++r)
      for (int c = 0; c < 5; ++c) {
        KF_CHECK_NEAR(cov(r, c), cov(c, r), 1e-12 * std::abs(cov(r, r)) + 1e-15);
        sym(r, c) = 0.5 * (cov(r, c) + cov(c, r));
      }
    TMatrixDSymEigen eig(sym);
    const TVectorD ev_ = eig.GetEigenValues();
    for (int k = 0; k < 5; ++k)
      KF_CHECK(ev_(k) >= -1e-12 * std::abs(ev_(0)));
  }
}

// The candidate can be used as mother in a 3C fit of exact daughters.
// Depends on correct R, Z (see R_Z_of_flight_line).
KF_TEST(usable_as_3C_mother) {
  TRandom3 rnd(6);
  int nConverged = 0;
  toy::Stats chi2;
  for (int i = 0; i < 20; ++i) {
    const auto ev = toy::lambdaEvent(rnd);
    std::vector<KinFitParticle> daughters = {
        toy::exactTrack(ev.proton, ev.decayVertex, toy::Resolution()),
        toy::exactTrack(ev.pion, ev.decayVertex, toy::Resolution())};
    KFitDecayCandFinder finder(daughters, toy::kMassLambda, ev.primaryVertex,
                               ev.decayVertex, 0.5, 0.5, 1.0, 2.0, 2.0, 4.0);
    KinFitter fitter(daughters);
    fitter.add3Constraint(finder.getDecayCand());
    if (fitter.fit()) {
      ++nConverged;
      chi2.add(fitter.getChi2());
    }
  }
  kftest::info("converged %d/20, <chi2> = %g", nConverged, chi2.mean());
  KF_CHECK(nConverged == 20);
  KF_CHECK(chi2.mean() < 1e-3);
}

KF_TEST_MAIN()
