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
 * KinFitFourMomentum.h
 *
 * Internal helper of KinFitter (not installed, not in the dictionary):
 * a four-momentum as plain doubles, used to sum tracks in the constraint
 * equations without TLorentzVector overhead.
 */

#pragma once

#include <TLorentzVector.h>

class KinFitTrackKinematics;

class KinFitFourMomentum {
public:
  double px = 0, py = 0, pz = 0, E = 0;

  KinFitFourMomentum() = default;
  explicit KinFitFourMomentum(const TLorentzVector &lv)
      : px(lv.Px()), py(lv.Py()), pz(lv.Pz()), E(lv.E()) {}

  /** Add (sign = +1) or subtract (sign = -1) a track. */
  void add(const KinFitTrackKinematics &k, double sign = 1.);

  /** Squared invariant mass E^2 - p^2. */
  double m2() const { return E * E - px * px - py * py - pz * pz; }
};
