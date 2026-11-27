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
 * KinFitTrackKinematics.h
 *
 * Internal helper of KinFitter (not installed, not in the dictionary):
 * kinematics of one track from its (1/p, theta, phi) fit parameters, with
 * the trigonometric functions evaluated once, and the derivatives of the
 * momentum and mass constraints with respect to these parameters.
 */

#pragma once

#include <TMatrixD.h>

class KinFitFourMomentum;

class KinFitTrackKinematics {
public:
  double p = 0;  ///< |p|
  double E = 0;  ///< energy for the track's mass
  double st = 0; ///< sin(theta)
  double ct = 1; ///< cos(theta)
  double sp = 0; ///< sin(phi)
  double cp = 1; ///< cos(phi)

  KinFitTrackKinematics() = default;
  /** @param m column vector of fit parameters
   *  @param offset row of 1/p of this track (theta, phi follow)
   *  @param mass track mass */
  KinFitTrackKinematics(const TMatrixD &m, int offset, double mass);

  /**
   * d(E^2 - p^2)/d(1/p, theta, phi) of this track entering the total P with
   * `sign` (+1: invariant mass of a sum, -1: missing mass P_init - sum).
   */
  void massDerivatives(const KinFitFourMomentum &P, double sign,
                       double &dInvP, double &dTheta, double &dPhi) const;

  /** d(px, py, pz, E)/d(1/p, theta, phi) into rows 0..3, columns
   * col0..col0+2 of H. */
  void momentumJacobian(TMatrixD &H, int col0) const;
};
