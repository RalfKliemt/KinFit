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

#include "KinFitFieldByXZ.h"

#include <algorithm>

KinFitFieldByXZ::KinFitFieldByXZ(const TH2 &byOfXZ, double scale)
    : fScale(scale) {
  TH2 *h = static_cast<TH2 *>(byOfXZ.Clone());
  h->SetDirectory(nullptr);
  fHist.reset(h);
}

TVector3 KinFitFieldByXZ::getField(double x, double, double z) const {
  const TAxis *ax = fHist->GetXaxis();
  const TAxis *az = fHist->GetYaxis();
  // TH2::Interpolate refuses points beyond the outer bin centres: clamp.
  const double xc = std::min(std::max(x, ax->GetBinCenter(1)),
                             ax->GetBinCenter(ax->GetNbins()));
  const double zc = std::min(std::max(z, az->GetBinCenter(1)),
                             az->GetBinCenter(az->GetNbins()));
  return TVector3(0., fScale * fHist->Interpolate(xc, zc), 0.);
}
