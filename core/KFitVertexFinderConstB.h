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
 * KFitVertexFinderConstB.h
 *
 * Vertex finder for charged tracks in a constant magnetic field Bz.
 */

#pragma once

#include "KinFitParticle.h"

#include <TMatrixD.h>
#include <TVector3.h>

#include <vector>

class KFitVertexFinderConstB {
private:
  /** @brief Per-track state used by the local helix linearization. */
  struct HelixTrack {
    TVector3 fRefPoint;    ///< Reference point on trajectory in mm.
    double fPhi;           ///< Azimuthal direction at reference point in rad.
    double fCotTheta;      ///< cot(theta), controls z advance per arc length.
    double fSignedRho;     ///< Signed curvature radius in mm.
    bool fUseStraightLine; ///< True for zero-curvature fallback.
  };

  int fVerbose = 0;              ///< Verbosity level.
  int fMaxIterations = 10;       ///< Max outer fixed-point iterations.
  int fMaxHelixNewtonSteps = 12; ///< Max Newton steps for point-to-helix projection.
  double fConvergenceTolMm = 1e-3; ///< Outer-loop convergence threshold in mm.
  double fBzTesla = 0.0;           ///< Constant z-component magnetic field in Tesla.

  std::vector<KinFitParticle> fCands; ///< Input candidate tracks.
  std::vector<HelixTrack> fTracks;    ///< Internal trajectory model.

  TVector3 fVertex; ///< Fitted vertex in mm.

  // Workspace matrices used by line-based least-squares update.
  TMatrixD fM;   ///< Per-track contribution matrix.
  TMatrixD fSys; ///< Accumulated normal-equation system matrix.
  TVector3 fB;   ///< Accumulated right-hand side.

  /** @brief Build helix or straight-line model for each input track. */
  void buildTrackModel();
  /** @brief Convert (R,Z,phi) track parametrization to a 3D reference point. */
  TVector3 lineBaseFromCand(const KinFitParticle &cand) const;
  /** @brief Local straight-line direction from (theta,phi). */
  TVector3 lineDirectionFromCand(const KinFitParticle &cand) const;

  /** @brief Position on helix at arc length s (mm). */
  TVector3 helixPoint(const HelixTrack &track, double sMm) const;
  /** @brief First derivative d(r)/d(s) at arc length s (unit tangent up to normalization). */
  TVector3 helixDerivative(const HelixTrack &track, double sMm) const;
  /** @brief Second derivative d2(r)/d(s)2 used by Newton projection. */
  TVector3 helixSecondDerivative(const HelixTrack &track, double sMm) const;

  /** @brief Estimate arc length s that minimizes distance to a space point. */
  double estimateArcLengthToPoint(const HelixTrack &track,
                                  const TVector3 &point) const;

  /** @brief Reset linear-system accumulators for one outer iteration. */
  void resetLinearSystem();
  /** @brief Add one linearized track constraint to the normal equations. */
  void addLineToVertexSystem(const TVector3 &r, const TVector3 &alpha,
                             double w = 1.0);
  /** @brief Solve current 3x3 normal equations and return fitted vertex. */
  TVector3 solveVertexFromSystem() const;
  /** @brief Run iterative fit loop with helix re-linearization. */
  void fitVertex();
  /** @brief True if v is the (-1000,-1000,-1000) failure marker. */
  static bool isInvalidVertex(const TVector3 &v);

public:
  /** @brief Constructor. Builds the track model and runs the fit once
   * with default settings.
   * @param cands Input tracks in KinFitParticle representation.
   * @param bzTesla Constant magnetic field Bz in Tesla.
   */
  KFitVertexFinderConstB(std::vector<KinFitParticle> &cands, double bzTesla);

  ~KFitVertexFinderConstB() {}

  /** @brief Rebuild the track model and rerun the fit.
   *
   * Call this after changing any setting (field, iterations, tolerance,
   * verbosity) for it to take effect.
   */
  void fit();

  /** @brief Set verbosity level (takes effect on next fit()). */
  void setVerbosity(int val) { fVerbose = val; }
  /** @brief Set max outer re-linearization iterations (next fit()). */
  void setMaxIterations(int val) { fMaxIterations = val; }
  /** @brief Set convergence tolerance on vertex shift in mm (next fit()). */
  void setConvergenceTolerance(double val) { fConvergenceTolMm = val; }
  /** @brief Set magnetic field Bz in Tesla (takes effect on next fit()). */
  void setMagneticFieldBz(double bzTesla) { fBzTesla = bzTesla; }

  /** @brief Return fitted vertex position in mm.
   * @return TVector3(x,y,z) of the fitted vertex, or (-1000,-1000,-1000)
   *         if the fit failed.
   */
  TVector3 getVertex() const { return fVertex; }
};

