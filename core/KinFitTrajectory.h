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
 * KinFitTrajectory.h
 *
 * Trajectory of a charged (or neutral) particle in a magnetic field,
 * parametrised by the signed arc length s from a reference state
 * (position, momentum, charge). Uniform fields (including B = 0) use the
 * exact helix around the field direction, other fields 4th-order
 * Runge-Kutta integration of dt/ds = kappa q / p (t x B(r)).
 *
 * Units: mm, MeV/c, Tesla, charge in units of e.
 */

#pragma once

#include <TVector3.h>

class KinFitField;

class KinFitTrajectory {
public:
  /** kappa = 0.299792458 (MeV/c) / (T mm): curvature 1/rho = kappa q B / p */
  static constexpr double kKappa = 0.299792458;

  /** @param pos reference position [mm]
   *  @param mom momentum at pos [MeV/c]
   *  @param charge charge [e]; 0 gives a straight line
   *  @param field magnetic field, not owned; nullptr means B = 0 */
  KinFitTrajectory(const TVector3 &pos, const TVector3 &mom, double charge,
                   const KinFitField *field);

  /** @brief State from KinFit track parameters: position
   * (-R sin(phi), R cos(phi), Z) with direction (theta, phi) and |p| = 1/invP
   * there. */
  static KinFitTrajectory fromTrackParameters(double invP, double theta,
                                              double phi, double R, double Z,
                                              double charge,
                                              const KinFitField *field);

  /** @brief Position, unit direction and its derivative d(dir)/ds at arc
   * length s. */
  void propagate(double s, TVector3 &pos, TVector3 &dir,
                 TVector3 *dDirDs = nullptr) const;

  /** @brief Point of closest approach to `point`, searched from the
   * reference state (local minimum of the distance).
   * @return arc length s of that point; pos/dir are filled. */
  double closestApproach(const TVector3 &point, TVector3 &pos,
                         TVector3 &dir) const;

  /** @brief KinFit track parameters at the point of closest approach to the
   * z axis in the transverse plane (r_T . t_T = 0) nearest to the reference
   * state. For a field along z this is the usual (1/p, theta, phi, R, Z) of a
   * helix; for other fields it is the state at that point.
   * @return false if the search did not converge. */
  bool trackParameters(double &invP, double &theta, double &phi, double &R,
                       double &Z) const;

  /** @brief Maximum Runge-Kutta step [mm] (non-uniform fields), default 5. */
  void setMaxStep(double val) { fMaxStep = val; }

  const TVector3 &getReferencePosition() const { return fPos; }
  const TVector3 &getReferenceDirection() const { return fDir; }
  double getMomentum() const { return fP; }
  double getCharge() const { return fCharge; }

private:
  TVector3 fPos, fDir;       // reference state, |fDir| = 1
  double fP = 0;             // |p| [MeV/c]
  double fCharge = 0;        // [e]
  const KinFitField *fField; // not owned, may be nullptr
  double fMaxStep = 5.;      // [mm]

  // Uniform field: helix around fAxis with curvature fOmega (signed) [1/mm]
  bool fUniform = true;
  TVector3 fAxis;     // unit field direction
  double fOmega = 0;  // kappa q |B| / p, 0 for straight lines
  double fScale = 0;  // kappa q / p, for the general equation of motion

  void propagateHelix(double s, TVector3 &pos, TVector3 &dir,
                      TVector3 *dDirDs) const;
  void propagateRK(double s, TVector3 &pos, TVector3 &dir,
                   TVector3 *dDirDs) const;
};
