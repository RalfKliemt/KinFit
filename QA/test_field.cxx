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

// Tests for the magnetic field models.

#include "KinFitTest.h"

#include "KinFitFieldByXZ.h"
#include "KinFitFieldByZ.h"
#include "KinFitFieldFunction.h"
#include "KinFitFieldUniform.h"

#include <TH1D.h>
#include <TH2D.h>

#include <memory>

KF_TEST(uniform) {
  const KinFitFieldUniform zero;
  KF_CHECK(zero.isUniform());
  KF_CHECK_NEAR(zero.getField(1, 2, 3).Mag(), 0, 0);
  const auto bz = KinFitFieldUniform::Bz(0.5);
  KF_CHECK_NEAR(bz.getField(10, -20, 300).Z(), 0.5, 0);
  const auto by = KinFitFieldUniform::By(-1.2);
  KF_CHECK_NEAR(by.getField(TVector3(0, 0, 0)).Y(), -1.2, 0);
  std::unique_ptr<KinFitField> c(by.clone());
  KF_CHECK_NEAR(c->getField(5, 5, 5).Y(), -1.2, 0);
}

// By(z): linear interpolation between bin centres, constant outside, and
// independent of the original histogram (copied).
KF_TEST(by_of_z_histogram) {
  auto h = std::make_unique<TH1D>("hBy", "", 10, 0., 1000.); // centres 50..950
  for (int i = 1; i <= 10; ++i)
    h->SetBinContent(i, 0.1 * i); // By = 0.1 at z = 50, 1.0 at z = 950
  const KinFitFieldByZ f(*h, 2.);
  h.reset();
  KF_CHECK(!f.isUniform());
  KF_CHECK_NEAR(f.getField(0, 0, 50.).Y(), 0.2, 1e-12);
  KF_CHECK_NEAR(f.getField(0, 0, 100.).Y(), 2 * 0.15, 1e-12);
  KF_CHECK_NEAR(f.getField(0, 0, -500.).Y(), 0.2, 1e-12);
  KF_CHECK_NEAR(f.getField(0, 0, 5000.).Y(), 2.0, 1e-12);
  KF_CHECK_NEAR(f.getField(7, 8, 100.).X(), 0, 0);
  KF_CHECK_NEAR(f.getField(7, 8, 100.).Z(), 0, 0);
  std::unique_ptr<KinFitField> c(f.clone());
  KF_CHECK_NEAR(c->getField(0, 0, 100.).Y(), 0.3, 1e-12);
}

// By(x, z): bilinear interpolation, clamped outside.
KF_TEST(by_of_xz_histogram) {
  TH2D h("hByXZ", "", 4, -200., 200., 5, 0., 1000.);
  for (int ix = 1; ix <= 4; ++ix)
    for (int iz = 1; iz <= 5; ++iz) {
      const double x = h.GetXaxis()->GetBinCenter(ix);
      const double z = h.GetYaxis()->GetBinCenter(iz);
      h.SetBinContent(ix, iz, 0.001 * x + 0.0005 * z); // bilinear itself
    }
  const KinFitFieldByXZ f(h);
  KF_CHECK_NEAR(f.getField(12., 0, 345.).Y(), 0.001 * 12 + 0.0005 * 345, 1e-12);
  // clamped to the outer bin centres x = +150, z = 900
  KF_CHECK_NEAR(f.getField(1000., 0, 5000.).Y(), 0.001 * 150 + 0.0005 * 900,
                1e-12);
}

KF_TEST(function) {
  const KinFitFieldFunction f(
      [](double x, double, double z) { return TVector3(0, 1e-3 * x, z); });
  KF_CHECK(!f.isUniform());
  KF_CHECK_NEAR(f.getField(100, 0, 7).Y(), 0.1, 1e-15);
  KF_CHECK_NEAR(f.getField(100, 0, 7).Z(), 7, 0);
  std::unique_ptr<KinFitField> c(f.clone());
  KF_CHECK_NEAR(c->getField(200, 0, 0).Y(), 0.2, 1e-15);
}

KF_TEST_MAIN()
