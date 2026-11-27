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

#include "KinFitTrackKinematics.h"
#include "KinFitFourMomentum.h"

#include <cmath>

KinFitTrackKinematics::KinFitTrackKinematics(const TMatrixD &m, int offset,
                                             double mass) {
  const double th = m(offset + 1, 0), ph = m(offset + 2, 0);
  p = 1. / m(offset, 0);
  E = std::sqrt(p * p + mass * mass);
  st = std::sin(th);
  ct = std::cos(th);
  sp = std::sin(ph);
  cp = std::cos(ph);
}

void KinFitTrackKinematics::massDerivatives(const KinFitFourMomentum &P,
                                            double sign, double &dInvP,
                                            double &dTheta,
                                            double &dPhi) const {
  // dE/d(1/p) = -p^3/E, dp/d(1/p) = -p^2 u
  dInvP = sign * (-2. * P.E * p * p * p / E +
                  2. * p * p * (P.px * st * cp + P.py * st * sp + P.pz * ct));
  dTheta = sign * (-2. * p * (P.px * ct * cp + P.py * ct * sp - P.pz * st));
  dPhi = sign * (-2. * p * st * (-P.px * sp + P.py * cp));
}

void KinFitTrackKinematics::momentumJacobian(TMatrixD &H, int col0) const {
  const double p2 = p * p;
  H(0, col0) = -p2 * st * cp;
  H(1, col0) = -p2 * st * sp;
  H(2, col0) = -p2 * ct;
  H(3, col0) = -p2 * p / E;
  H(0, col0 + 1) = p * ct * cp;
  H(1, col0 + 1) = p * ct * sp;
  H(2, col0 + 1) = -p * st;
  H(0, col0 + 2) = -p * st * sp;
  H(1, col0 + 2) = p * st * cp;
}
