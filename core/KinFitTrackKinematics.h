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
 * four-momentum of one particle from its fit parameters, with the
 * derivatives with respect to its five parameter slots and, for neutral
 * hits, with respect to the vertex.
 *
 * Parameter slots: tracks (1/p, theta, phi, R, Z), neutral hits
 * (E, X, Y, Z, unused) with direction (hit - vertex).
 */

#pragma once

#include <TMatrixD.h>
#include <TVector3.h>

class KinFitFourMomentum;

class KinFitTrackKinematics {
public:
  double px = 0, py = 0, pz = 0, E = 0;
  double dEta[4][5] = {}; ///< d(px, py, pz, E) / d(parameter slot)
  double dVtx[4][3] = {}; ///< d(px, py, pz, E) / d(vertex), neutral hits only

  KinFitTrackKinematics() = default;

  /** @brief Track: (1/p, theta, phi, R, Z) at rows offset.. of m. */
  KinFitTrackKinematics(const TMatrixD &m, int offset, double mass);

  /** @brief Neutral hit: (E, X, Y, Z) at rows offset.. of m, produced at
   * `vertex`. */
  static KinFitTrackKinematics neutralHit(const TMatrixD &m, int offset,
                                          double mass, const TVector3 &vertex);

  /**
   * d(E_tot^2 - p_tot^2) of the total P to which this particle contributes
   * with `sign` (+1: invariant mass of a sum, -1: missing mass P_init - sum),
   * with respect to its parameter slots and to the vertex.
   */
  void massDerivatives(const KinFitFourMomentum &P, double sign,
                       double dEtaOut[5], double dVtxOut[3]) const;
};
