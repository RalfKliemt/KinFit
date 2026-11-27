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
 * KinFitFieldByXZ.h
 *
 * Dipole-like field B = (0, By(x, z), 0) with By(x, z) taken from a TH2
 * (x axis: x in mm, y axis: z in mm, content: By in Tesla), bilinearly
 * interpolated and constant beyond the outermost bin centres.
 */

#pragma once

#include "KinFitField.h"

#include <TH2.h>

#include <memory>

class KinFitFieldByXZ : public KinFitField {
private:
  std::shared_ptr<const TH2> fHist; // private copy, not attached to a file
  double fScale = 1.;

public:
  /** @param byOfXZ histogram By(x, z) [T]; x axis = x [mm], y axis = z [mm];
   *  copied
   *  @param scale factor applied to the histogram content (e.g. polarity) */
  explicit KinFitFieldByXZ(const TH2 &byOfXZ, double scale = 1.);

  TVector3 getField(double x, double y, double z) const override;
  using KinFitField::getField;
  KinFitField *clone() const override { return new KinFitFieldByXZ(*this); }
};
