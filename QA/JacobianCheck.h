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

// Comparison of an analytic Jacobian with central finite differences.

#pragma once

#include "KinFitTest.h"

#include <TMatrixD.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>

/** Compare an analytic Jacobian with central finite differences of f;
 *  step 1e-6 * max(|x|, minScale). */

inline void compareJacobian(const char *what, const TMatrixD &analytic,
                            const std::function<TMatrixD(const TMatrixD &)> &f,
                            const TMatrixD &x0, double minScale = 1e-3) {
  const TMatrixD f0 = f(x0);
  KF_CHECK(analytic.GetNrows() == f0.GetNrows());
  KF_CHECK(analytic.GetNcols() == x0.GetNrows());
  if (analytic.GetNrows() != f0.GetNrows() ||
      analytic.GetNcols() != x0.GetNrows())
    return;

  int nBad = 0;
  for (int k = 0; k < x0.GetNrows(); ++k) {
    const double h = 1e-6 * std::max(std::abs(x0(k, 0)), minScale);
    TMatrixD xp(x0), xm(x0);
    xp(k, 0) += h;
    xm(k, 0) -= h;
    const TMatrixD fp = f(xp), fm = f(xm);
    for (int r = 0; r < f0.GetNrows(); ++r) {
      // Scale: largest analytic entry of this row (derivatives of one
      // constraint share its units).
      double scale = 0;
      for (int cc = 0; cc < analytic.GetNcols(); ++cc)
        scale = std::max(scale, std::abs(analytic(r, cc)));
      const double numeric = (fp(r, 0) - fm(r, 0)) / (2 * h);
      const double tol = 1e-5 * std::max(scale, 1e-12) +
                         1e-6 * std::abs(numeric);
      if (std::abs(numeric - analytic(r, k)) > tol) {
        if (++nBad <= 5)
          kftest::info("%s: d f[%d] / d x[%d]: analytic %.6g, numeric %.6g",
                       what, r, k, analytic(r, k), numeric);
      }
    }
  }
  if (nBad)
    kftest::reportFailure(__FILE__, __LINE__,
                          std::string(what) + ": " + std::to_string(nBad) +
                              " Jacobian entries differ");
}

