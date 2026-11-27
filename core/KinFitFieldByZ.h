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
 * KinFitFieldByZ.h
 *
 * Dipole-like field B = (0, By(z), 0) with By(z) taken from a TH1 (x axis:
 * z in mm, content: By in Tesla), linearly interpolated between bin
 * centres and constant beyond the first/last bin centre.
 */

#pragma once

#include "KinFitField.h"

#include <TH1.h>

#include <memory>

class KinFitFieldByZ : public KinFitField {
private:
  std::shared_ptr<const TH1> fHist; // private copy, not attached to a file
  double fScale = 1.;

public:
  /** @param byOfZ histogram By(z) [T] over z [mm]; copied
   *  @param scale factor applied to the histogram content (e.g. polarity) */
  explicit KinFitFieldByZ(const TH1 &byOfZ, double scale = 1.);

  TVector3 getField(double x, double y, double z) const override;
  using KinFitField::getField;
  KinFitField *clone() const override { return new KinFitFieldByZ(*this); }
};
