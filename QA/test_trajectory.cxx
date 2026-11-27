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

// Tests for KinFitTrajectory: helix vs Runge-Kutta, closest approach,
// track-parameter round trip.

#include "KinFitTest.h"
#include "ToyTracks.h"

#include "KinFitFieldByZ.h"
#include "KinFitFieldFunction.h"
#include "KinFitFieldUniform.h"
#include "KinFitTrajectory.h"

#include <TH1D.h>
#include <TRandom3.h>

namespace {
/** Same constant field, but declared non-uniform: forces Runge-Kutta. */
KinFitFieldFunction asNumeric(const TVector3 &B) {
  return KinFitFieldFunction([B](double, double, double) { return B; });
}

TVector3 randomMomentum(TRandom3 &rnd) {
  TVector3 p;
  p.SetMagThetaPhi(rnd.Uniform(200, 1500), rnd.Uniform(0.2, 1.4),
                   rnd.Uniform(-TMath::Pi(), TMath::Pi()));
  return p;
}
} // namespace

// p = 300 MeV/c perpendicular to Bz = 1 T: rho = 300 / 0.29979 mm; a
// positive charge turns clockwise (towards -y for p along +x).
KF_TEST(radius_and_turning_direction) {
  const auto B = KinFitFieldUniform::Bz(1.);
  const KinFitTrajectory t(TVector3(0, 0, 0), TVector3(300, 0, 0), +1, &B);
  const double rho = 300. / KinFitTrajectory::kKappa;
  TVector3 pos, dir;
  t.propagate(TMath::Pi() * rho, pos, dir); // half turn
  KF_CHECK_NEAR(pos.X(), 0, 1e-9);
  KF_CHECK_NEAR(pos.Y(), -2 * rho, 1e-9);
  KF_CHECK_NEAR(dir.X(), -1, 1e-12);
  const KinFitTrajectory n(TVector3(0, 0, 0), TVector3(300, 0, 0), -1, &B);
  n.propagate(TMath::Pi() * rho, pos, dir);
  KF_CHECK_NEAR(pos.Y(), 2 * rho, 1e-9);
}

// B = 0 and neutral particles: straight lines.
KF_TEST(straight_lines) {
  const KinFitFieldUniform zero;
  const auto Bz = KinFitFieldUniform::Bz(1.);
  const TVector3 p(100, 200, 300), r0(1, 2, 3);
  for (const KinFitTrajectory &t :
       {KinFitTrajectory(r0, p, 1, &zero), KinFitTrajectory(r0, p, 0, &Bz),
        KinFitTrajectory(r0, p, 1, nullptr)}) {
    TVector3 pos, dir;
    t.propagate(250., pos, dir);
    KF_CHECK_NEAR((pos - (r0 + 250. * p.Unit())).Mag(), 0, 1e-12);
  }
}

// Exact helix == Runge-Kutta for the same constant field (Bz, By, oblique).
KF_TEST(helix_equals_runge_kutta) {
  TRandom3 rnd(1);
  int nBad = 0;
  for (const TVector3 &B : {TVector3(0, 0, 1.), TVector3(0, -0.8, 0),
                            TVector3(0.3, 0.5, 0.7)}) {
    const KinFitFieldUniform helix(B.X(), B.Y(), B.Z());
    const KinFitFieldFunction rk = asNumeric(B);
    for (int i = 0; i < 20; ++i) {
      const TVector3 p = randomMomentum(rnd);
      const TVector3 r0(rnd.Gaus(0, 5), rnd.Gaus(0, 5), rnd.Gaus(0, 20));
      const double q = rnd.Rndm() < 0.5 ? -1 : 1;
      const KinFitTrajectory th(r0, p, q, &helix), tr(r0, p, q, &rk);
      for (double s : {-300., 50., 1000.}) {
        TVector3 ph, dh, pr, dr, dth, dtr;
        th.propagate(s, ph, dh, &dth);
        tr.propagate(s, pr, dr, &dtr);
        if ((ph - pr).Mag() > 1e-5 || (dh - dr).Mag() > 1e-8 ||
            (dth - dtr).Mag() > 1e-8) {
          if (++nBad <= 5)
            kftest::info("B (%g,%g,%g) s %g: |dpos| %g |ddir| %g", B.X(),
                         B.Y(), B.Z(), s, (ph - pr).Mag(), (dh - dr).Mag());
        }
        KF_CHECK_NEAR(dh.Mag(), 1, 1e-12);
      }
    }
  }
  KF_CHECK(nBad == 0);
}

// Closest approach: residual perpendicular to the direction and a local
// minimum of the distance, for a uniform and a By(z) field.
KF_TEST(closest_approach) {
  TH1D h("hBy", "", 50, -500., 2500.);
  for (int i = 1; i <= 50; ++i) {
    const double z = h.GetBinCenter(i);
    h.SetBinContent(i, 1.0 * std::exp(-0.5 * std::pow((z - 800.) / 400., 2)));
  }
  const KinFitFieldByZ byz(h);
  const auto bz = KinFitFieldUniform::Bz(0.5);
  TRandom3 rnd(2);
  for (const KinFitField *f : {static_cast<const KinFitField *>(&bz),
                               static_cast<const KinFitField *>(&byz)}) {
    for (int i = 0; i < 30; ++i) {
      const KinFitTrajectory t(TVector3(rnd.Gaus(0, 3), rnd.Gaus(0, 3), 0),
                               randomMomentum(rnd), rnd.Rndm() < .5 ? -1 : 1,
                               f);
      TVector3 onTrack, dir;
      t.propagate(rnd.Uniform(-100, 800), onTrack, dir);
      const TVector3 point = onTrack + TVector3(rnd.Gaus(0, 2), rnd.Gaus(0, 2),
                                                rnd.Gaus(0, 2));
      TVector3 pca, pdir;
      const double s = t.closestApproach(point, pca, pdir);
      KF_CHECK_NEAR((pca - point).Dot(pdir), 0, 1e-8);
      TVector3 a, b, da, db;
      t.propagate(s - 0.5, a, da);
      t.propagate(s + 0.5, b, db);
      KF_CHECK((a - point).Mag() >= (pca - point).Mag());
      KF_CHECK((b - point).Mag() >= (pca - point).Mag());
    }
  }
}

// trackParameters(): transverse PCA to the z axis; fromTrackParameters()
// of the result describes the same trajectory.
KF_TEST(track_parameter_round_trip) {
  TRandom3 rnd(3);
  for (const TVector3 &B : {TVector3(0, 0, 1.), TVector3(0, 0.7, 0)}) {
    const KinFitFieldUniform f(B.X(), B.Y(), B.Z());
    for (int i = 0; i < 30; ++i) {
      const double q = rnd.Rndm() < .5 ? -1 : 1;
      const TVector3 v(rnd.Gaus(0, 5), rnd.Gaus(0, 5), rnd.Gaus(20, 10));
      const KinFitTrajectory t(v, randomMomentum(rnd), q, &f);
      double invP, theta, phi, R, Z;
      KF_CHECK(t.trackParameters(invP, theta, phi, R, Z));
      const KinFitTrajectory u =
          KinFitTrajectory::fromTrackParameters(invP, theta, phi, R, Z, q, &f);
      const TVector3 &r0 = u.getReferencePosition();
      const TVector3 &t0 = u.getReferenceDirection();
      KF_CHECK_NEAR(r0.X() * t0.X() + r0.Y() * t0.Y(), 0, 1e-8);
      // u passes through the original vertex v
      TVector3 pca, dir;
      u.closestApproach(v, pca, dir);
      KF_CHECK_NEAR((pca - v).Mag(), 0, 1e-6);
    }
  }
}

// Independent check: trackParameters() in constant Bz agrees with the
// closed-form helix PCA of the toy generator.
KF_TEST(track_parameters_match_toy_helix) {
  TRandom3 rnd(4);
  const double bz = 0.8;
  const auto f = KinFitFieldUniform::Bz(bz);
  for (int i = 0; i < 50; ++i) {
    const int q = rnd.Rndm() < .5 ? -1 : 1;
    const TVector3 v(rnd.Gaus(0, 5), rnd.Gaus(0, 5), rnd.Gaus(20, 10));
    const TVector3 p = randomMomentum(rnd);
    const KinFitParticle ref = toy::helixTrackAtPCA(v, p, q, bz, 139.57);
    double invP, theta, phi, R, Z;
    KinFitTrajectory(v, p, q, &f).trackParameters(invP, theta, phi, R, Z);
    KF_CHECK_NEAR(theta, ref.Theta(), 1e-9);
    KF_CHECK_NEAR(TVector2::Phi_mpi_pi(phi - ref.Phi()), 0, 1e-9);
    KF_CHECK_NEAR(R, ref.getR(), 1e-6);
    KF_CHECK_NEAR(Z, ref.getZ(), 1e-6);
  }
}

KF_TEST_MAIN()
