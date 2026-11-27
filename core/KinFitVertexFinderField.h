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
 * KinFitVertexFinderField.h
 *
 * Vertex finder for all tracks in any field model (KinFitField): each
 * trajectory is linearised by its tangent at the point of closest approach
 * to the current vertex estimate, and the 3x3 least-squares system of the
 * tangent lines is solved; repeated until the vertex moves less than the
 * tolerance. Neutral hits are skipped (their direction depends on the
 * vertex). Unweighted: meant as start value (KinFitter's common-vertex
 * constraint gives the weighted vertex and its covariance).
 */

#pragma once

#include "KinFitParticle.h"

#include <TVector3.h>

#include <vector>

class KinFitField;

class KinFitVertexFinderField {
private:
  TVector3 fVertex{-1000., -1000., -1000.};
  bool fValid = false;
  bool fHasEstimate = false;
  int fNTracks = 0;
  int fIterations = 0;

public:
  /** @param cands tracks; charged ones follow the field, charge 0 straight
   *         lines, neutral hits are ignored
   *  @param field magnetic field, nullptr = B = 0
   *  @param maxIterations maximum number of re-linearisations
   *  @param tolerance convergence threshold on the vertex shift [mm] */
  KinFitVertexFinderField(const std::vector<KinFitParticle> &cands,
                          const KinFitField *field = nullptr,
                          int maxIterations = 20, double tolerance = 1e-4);

  /** @brief Vertex [mm]; the last estimate if not converged,
   * (-1000, -1000, -1000) if there is none. */
  TVector3 getVertex() const { return fVertex; }
  /** @brief True if converged (shift below the tolerance). */
  bool isValid() const { return fValid; }
  /** @brief True if a vertex estimate exists, converged or not (false for
   * fewer than two usable tracks or a singular system). */
  bool hasEstimate() const { return fHasEstimate; }
  int getNumberOfTracks() const { return fNTracks; }
  int getIterations() const { return fIterations; }
};
