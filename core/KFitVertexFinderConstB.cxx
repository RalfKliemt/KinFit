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

#include "KFitVertexFinderConstB.h"

#include <TMath.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <iostream>

/**
 * @file KFitVertexFinderConstB.cxx
 * @brief Helix-aware vertex fit in a constant magnetic field Bz.
 *
 * The algorithm alternates between:
 * 1) projecting the current vertex guess onto each track helix,
 * 2) linearizing each helix by its local tangent,
 * 3) solving a 3x3 least-squares vertex system.
 */

/** @brief Build and run the fit immediately from input tracks and field. */
KFitVertexFinderConstB::KFitVertexFinderConstB(std::vector<KinFitParticle> &cands,
                                               double bzTesla)
    : fBzTesla(bzTesla), fCands(cands) {
  fM.ResizeTo(3, 3);
  fSys.ResizeTo(3, 3);

  fit();
}

/** @brief Rebuild the track model with current settings and refit. */
void KFitVertexFinderConstB::fit() {
  buildTrackModel();
  fitVertex();
}

bool KFitVertexFinderConstB::isInvalidVertex(const TVector3 &v) {
  return v.X() < -999.0 && v.Y() < -999.0 && v.Z() < -999.0;
}

/** @brief Build a 3D point on the track from the (R,phi,Z) representation. */
TVector3 KFitVertexFinderConstB::lineBaseFromCand(const KinFitParticle &cand) const {
  return TVector3(cand.getR() * std::cos(cand.Phi() + TMath::PiOver2()),
                  cand.getR() * std::sin(cand.Phi() + TMath::PiOver2()),
                  cand.getZ());
}

/** @brief Build local track direction from spherical angles. */
TVector3 KFitVertexFinderConstB::lineDirectionFromCand(const KinFitParticle &cand) const {
  return TVector3(std::sin(cand.Theta()) * std::cos(cand.Phi()),
                  std::sin(cand.Theta()) * std::sin(cand.Phi()),
                  std::cos(cand.Theta()));
}

/**
 * @brief Convert input candidates into either helix or straight-line models.
 *
 * Straight-line fallback is used for neutral tracks, tiny Bz, or degenerate pT.
 */
void KFitVertexFinderConstB::buildTrackModel() {
  fTracks.clear();
  fTracks.reserve(fCands.size());

  for (size_t i = 0; i < fCands.size(); ++i) {
    const KinFitParticle &cand = fCands[i];
    HelixTrack track;

    track.fRefPoint = lineBaseFromCand(cand);
    track.fPhi = cand.Phi();

    const double sinTheta = std::sin(cand.Theta());
    if (std::abs(sinTheta) < 1e-9) {
      track.fCotTheta = 0.0;
      track.fUseStraightLine = true;
      track.fSignedRho = 0.0;
      fTracks.push_back(track);
      continue;
    }

    track.fCotTheta = std::cos(cand.Theta()) / sinTheta;

    const double q = cand.getCharge();
    const double ptMeV = cand.P() * sinTheta;
    const double absB = std::abs(fBzTesla);

    if (absB < 1e-12 || std::abs(q) < 1e-12 || ptMeV < 1e-9) {
      track.fUseStraightLine = true;
      track.fSignedRho = 0.0;
      fTracks.push_back(track);
      continue;
    }

    // rho[mm] = pT[MeV/c] / (0.299792458 * |q| * B[T]).
    const double rhoAbsMm = ptMeV / (0.299792458 * std::abs(q) * absB);
    const double chargeFieldSign = (q * fBzTesla >= 0.0) ? 1.0 : -1.0;

    // Signed radius defines left/right turning sense at increasing arc length.
    track.fSignedRho = -chargeFieldSign * rhoAbsMm;
    track.fUseStraightLine = false;
    fTracks.push_back(track);
  }
}

/** @brief Evaluate helix position at arc length s from the reference point. */
TVector3 KFitVertexFinderConstB::helixPoint(const HelixTrack &track, double sMm) const {
  if (track.fUseStraightLine || std::abs(track.fSignedRho) < 1e-12) {
    const TVector3 dir(std::cos(track.fPhi), std::sin(track.fPhi), track.fCotTheta);
    return track.fRefPoint + sMm * dir;
  }

  const double phiS = track.fPhi + sMm / track.fSignedRho;
  const double sinPhi0 = std::sin(track.fPhi);
  const double cosPhi0 = std::cos(track.fPhi);
  const double signedR = track.fSignedRho;

  const double cx = track.fRefPoint.X() + signedR * (-sinPhi0);
  const double cy = track.fRefPoint.Y() + signedR * (cosPhi0);

  return TVector3(cx + signedR * std::sin(phiS),
                  cy - signedR * std::cos(phiS),
                  track.fRefPoint.Z() + sMm * track.fCotTheta);
}

/** @brief Evaluate first derivative dr/ds at arc length s. */
TVector3 KFitVertexFinderConstB::helixDerivative(const HelixTrack &track,
                                                 double sMm) const {
  if (track.fUseStraightLine || std::abs(track.fSignedRho) < 1e-12) {
    return TVector3(std::cos(track.fPhi), std::sin(track.fPhi), track.fCotTheta);
  }

  const double phiS = track.fPhi + sMm / track.fSignedRho;
  return TVector3(std::cos(phiS), std::sin(phiS), track.fCotTheta);
}

/** @brief Evaluate second derivative d2r/ds2 at arc length s. */
TVector3 KFitVertexFinderConstB::helixSecondDerivative(const HelixTrack &track,
                                                       double sMm) const {
  if (track.fUseStraightLine || std::abs(track.fSignedRho) < 1e-12) {
    return TVector3(0.0, 0.0, 0.0);
  }

  const double phiS = track.fPhi + sMm / track.fSignedRho;
  return TVector3(-std::sin(phiS) / track.fSignedRho,
                  std::cos(phiS) / track.fSignedRho,
                  0.0);
}

/**
 * @brief Project a point onto a track: closest approach in arc length.
 *
 * Straight lines are projected exactly. For helices the start value comes
 * from the exact closest point on the transverse circle, choosing the turn
 * whose z is closest to the point (handles loopers). Newton iterations on
 * the full 3D distance then refine it; steps are only taken where the
 * distance has positive curvature (to avoid converging onto a maximum)
 * and are limited to a fraction of the radius to stay on the same branch.
 * @return Estimated arc length s in mm.
 */
double KFitVertexFinderConstB::estimateArcLengthToPoint(const HelixTrack &track,
                                                         const TVector3 &point) const {
  const TVector3 delta0 = point - track.fRefPoint;

  if (track.fUseStraightLine || std::abs(track.fSignedRho) < 1e-12) {
    const TVector3 dir = helixDerivative(track, 0.0);
    return delta0.Dot(dir) / dir.Mag2();
  }

  const double rho = track.fSignedRho;
  const double absRho = std::abs(rho);
  const double sinPhi0 = std::sin(track.fPhi);
  const double cosPhi0 = std::cos(track.fPhi);
  const double cx = track.fRefPoint.X() - rho * sinPhi0;
  const double cy = track.fRefPoint.Y() + rho * cosPhi0;

  // Angles of reference point and target point seen from the circle centre.
  const double aRef = std::atan2(track.fRefPoint.Y() - cy, track.fRefPoint.X() - cx);
  double s = 0.0;
  if (std::hypot(point.X() - cx, point.Y() - cy) > 1e-9) {
    const double aPt = std::atan2(point.Y() - cy, point.X() - cx);
    // Position angle around the centre advances by s/rho, like the direction.
    const double dA = TVector2::Phi_mpi_pi(aPt - aRef);
    s = rho * dA;
  }

  // Pick the helix turn closest in z.
  const double zPerTurn = TMath::TwoPi() * absRho * track.fCotTheta;
  if (std::abs(zPerTurn) > 1e-9) {
    const double zRes = point.Z() - (track.fRefPoint.Z() + s * track.fCotTheta);
    s += std::round(zRes / zPerTurn) * TMath::TwoPi() * absRho;
  }

  const double maxStep = 0.5 * absRho;
  for (int iter = 0; iter < fMaxHelixNewtonSteps; ++iter) {
    const TVector3 r = helixPoint(track, s);
    const TVector3 dr = helixDerivative(track, s);
    const TVector3 d2r = helixSecondDerivative(track, s);
    const TVector3 diff = r - point;

    const double g = diff.Dot(dr);
    const double gp = dr.Mag2() + diff.Dot(d2r);
    if (gp <= 1e-12) {
      break; // Not a minimum locally; keep the geometric start value.
    }

    const double step = std::clamp(g / gp, -maxStep, maxStep);
    s -= step;

    if (std::abs(step) < 1e-4) {
      break;
    }
  }

  return s;
}

/** @brief Clear normal-equation accumulators before building one fit step. */
void KFitVertexFinderConstB::resetLinearSystem() {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      fM(i, j) = 0.0;
      fSys(i, j) = 0.0;
    }
  }
  fB.SetXYZ(0.0, 0.0, 0.0);
}

/** @brief Add one local tangent-line contribution to the least-squares system. */
void KFitVertexFinderConstB::addLineToVertexSystem(const TVector3 &r,
                                                    const TVector3 &alpha,
                                                    double w) {
  fM(0, 0) = w * (alpha.Y() * alpha.Y() + alpha.Z() * alpha.Z());
  fM(0, 1) = -w * (alpha.X() * alpha.Y());
  fM(0, 2) = -w * (alpha.X() * alpha.Z());
  fM(1, 1) = w * (alpha.X() * alpha.X() + alpha.Z() * alpha.Z());
  fM(1, 2) = -w * (alpha.Y() * alpha.Z());
  fM(2, 2) = w * (alpha.Y() * alpha.Y() + alpha.X() * alpha.X());

  fSys(0, 0) += fM(0, 0);
  fSys(0, 1) += fM(0, 1);
  fSys(0, 2) += fM(0, 2);
  fSys(1, 1) += fM(1, 1);
  fSys(1, 2) += fM(1, 2);
  fSys(2, 2) += fM(2, 2);

  fB.SetX(fB.X() + fM(0, 0) * r.X() + fM(0, 1) * r.Y() + fM(0, 2) * r.Z());
  fB.SetY(fB.Y() + fM(0, 1) * r.X() + fM(1, 1) * r.Y() + fM(1, 2) * r.Z());
  fB.SetZ(fB.Z() + fM(0, 2) * r.X() + fM(1, 2) * r.Y() + fM(2, 2) * r.Z());
}

/**
 * @brief Solve the symmetric 3x3 normal-equation system.
 * @return Fitted vertex; returns (-1000,-1000,-1000) for singular systems.
 */
TVector3 KFitVertexFinderConstB::solveVertexFromSystem() const {
  TMatrixD m(3, 3);

  const double det =
      fSys(0, 0) * fSys(1, 1) * fSys(2, 2) +
      fSys(0, 1) * fSys(1, 2) * fSys(0, 2) +
      fSys(0, 1) * fSys(1, 2) * fSys(0, 2) -
      fSys(0, 2) * fSys(1, 1) * fSys(0, 2) -
      fSys(0, 1) * fSys(0, 1) * fSys(2, 2) -
      fSys(1, 2) * fSys(1, 2) * fSys(0, 0);

  if (std::abs(det) < 1e-18) {
    return TVector3(-1000.0, -1000.0, -1000.0);
  }

  m(0, 0) = fSys(1, 1) * fSys(2, 2) - fSys(1, 2) * fSys(1, 2);
  m(0, 1) = fSys(1, 2) * fSys(0, 2) - fSys(0, 1) * fSys(2, 2);
  m(0, 2) = fSys(0, 1) * fSys(1, 2) - fSys(1, 1) * fSys(0, 2);
  m(1, 1) = fSys(0, 0) * fSys(2, 2) - fSys(0, 2) * fSys(0, 2);
  m(1, 2) = fSys(0, 1) * fSys(0, 2) - fSys(0, 0) * fSys(1, 2);
  m(2, 2) = fSys(0, 0) * fSys(1, 1) - fSys(0, 1) * fSys(0, 1);
  m(1, 0) = m(0, 1);
  m(2, 0) = m(0, 2);
  m(2, 1) = m(1, 2);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      m(i, j) /= det;
    }
  }

  return m * fB;
}

/**
 * @brief Iterative helix-linearized vertex fit.
 *
 * Starts from a straight-line estimate and then re-linearizes each helix
 * around the current vertex guess until convergence.
 */
void KFitVertexFinderConstB::fitVertex() {
  if (fTracks.empty()) {
    fVertex.SetXYZ(-1000.0, -1000.0, -1000.0);
    return;
  }

  // Initial estimate from straight-line approximation around reference points.
  resetLinearSystem();
  TVector3 refCentroid(0.0, 0.0, 0.0);
  for (size_t i = 0; i < fTracks.size(); ++i) {
    const TVector3 base = fTracks[i].fRefPoint;
    refCentroid += base;
    TVector3 dir = helixDerivative(fTracks[i], 0.0);
    if (dir.Mag2() < 1e-18) {
      continue;
    }
    dir *= 1.0 / std::sqrt(dir.Mag2());
    addLineToVertexSystem(base, dir);
  }
  refCentroid *= 1.0 / fTracks.size();

  TVector3 vertexGuess = solveVertexFromSystem();
  if (isInvalidVertex(vertexGuess)) {
    // Degenerate straight-line system (e.g. parallel tracks): start from the
    // centroid of the reference points instead of the failure marker.
    vertexGuess = refCentroid;
  }

  if (fVerbose > 0) {
    std::cout << "KFitVertexFinderConstB::fitVertex(): " << fTracks.size()
              << " tracks, Bz = " << fBzTesla << " T, start ("
              << vertexGuess.X() << ", " << vertexGuess.Y() << ", "
              << vertexGuess.Z() << ") mm" << std::endl;
  }

  for (int iter = 0; iter < std::max(1, fMaxIterations); ++iter) {
    resetLinearSystem();

    for (size_t i = 0; i < fTracks.size(); ++i) {
      const HelixTrack &track = fTracks[i];
      // Closest point on current helix to current vertex guess.
      const double s = estimateArcLengthToPoint(track, vertexGuess);
      const TVector3 point = helixPoint(track, s);

      // Linearize helix locally by tangent and add to normal equations.
      TVector3 tangent = helixDerivative(track, s);
      if (tangent.Mag2() < 1e-18) {
        continue;
      }
      tangent *= 1.0 / std::sqrt(tangent.Mag2());
      addLineToVertexSystem(point, tangent);
    }

    const TVector3 newVertex = solveVertexFromSystem();
    if (isInvalidVertex(newVertex)) {
      if (fVerbose > 0) {
        std::cout << "KFitVertexFinderConstB::fitVertex(): singular system in "
                  << "iteration " << iter << ", fit failed" << std::endl;
      }
      fVertex = newVertex;
      return;
    }

    const double shift = (newVertex - vertexGuess).Mag();
    if (fVerbose > 1) {
      std::cout << "  iter " << iter << ": (" << newVertex.X() << ", "
                << newVertex.Y() << ", " << newVertex.Z()
                << ") mm, shift " << shift << " mm" << std::endl;
    }

    vertexGuess = newVertex;
    if (shift < fConvergenceTolMm) {
      if (fVerbose > 0) {
        std::cout << "KFitVertexFinderConstB::fitVertex(): converged after "
                  << iter + 1 << " iterations" << std::endl;
      }
      fVertex = vertexGuess;
      return;
    }
  }

  if (fVerbose > 0) {
    std::cout << "KFitVertexFinderConstB::fitVertex(): not converged after "
              << std::max(1, fMaxIterations) << " iterations" << std::endl;
  }
  fVertex = vertexGuess;
}
