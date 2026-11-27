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

#include "KinFitFieldByZ.h"

#include <algorithm>

KinFitFieldByZ::KinFitFieldByZ(const TH1 &byOfZ, double scale)
    : fScale(scale) {
  TH1 *h = static_cast<TH1 *>(byOfZ.Clone());
  h->SetDirectory(nullptr);
  fHist.reset(h);
}

TVector3 KinFitFieldByZ::getField(double, double, double z) const {
  const TAxis *ax = fHist->GetXaxis();
  // Clamp to the range of bin centres (constant extrapolation).
  const double lo = ax->GetBinCenter(1), hi = ax->GetBinCenter(ax->GetNbins());
  const double zc = std::min(std::max(z, lo), hi);
  return TVector3(0., fScale * fHist->Interpolate(zc), 0.);
}
