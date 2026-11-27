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
 * KinFitTrackLine.h
 *
 * Internal helper of KinFitter (not installed, not in the dictionary):
 * straight track line B + s d with B = (-R sin(phi), R cos(phi), Z), and the
 * vertex constraint (d1 x d2).(B1 - B2) of two lines with its derivatives.
 */

#pragma once

#include <TMatrixD.h>
#include <TVector3.h>

class KinFitTrackKinematics;

class KinFitTrackLine {
public:
  TVector3 b;         ///< base point (point of closest approach to z)
  TVector3 d;         ///< unit direction
  TVector3 dd_dTheta; ///< d(d)/d(theta)
  TVector3 dd_dPhi;   ///< d(d)/d(phi)
  TVector3 db_dPhi;   ///< d(B)/d(phi)
  TVector3 db_dR;     ///< d(B)/d(R)

  /** @param m column vector of fit parameters
   *  @param offset row of 1/p of this track (theta, phi, R, Z follow)
   *  @param k kinematics of the same track */
  KinFitTrackLine(const TMatrixD &m, int offset,
                  const KinFitTrackKinematics &k);

  /** (d1 x d2).(B1 - B2): zero if the lines intersect. */
  static double vertexConstraint(const KinFitTrackLine &l1,
                                 const KinFitTrackLine &l2);

  /** Derivatives of vertexConstraint() into row `row` of H; the parameters
   * (1/p, theta, phi, R, Z) of l1 start at column col1, those of l2 at col2.
   */
  static void vertexJacobian(TMatrixD &H, int row, const KinFitTrackLine &l1,
                             int col1, const KinFitTrackLine &l2, int col2);
};
