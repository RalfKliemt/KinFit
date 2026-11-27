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

#include "KinFitParticle.h"

#include <cmath>
#include <iostream>

/**
 * @file KinFitParticle.cxx
 * @brief Implementation of particle parameter and covariance container.
 */

KinFitParticle::KinFitParticle(TLorentzVector cand, double R, double Z)
    : TLorentzVector(cand), cand(cand), fMomentum(cand.P()),
      fTheta(cand.Theta()), fPhi(cand.Phi()), fR(R), fZ(Z) {
  fCov.ResizeTo(5, 5);
}

KinFitParticle::KinFitParticle(TLorentzVector cand, double X, double Y,
                               double Z)
    : TLorentzVector(cand), cand(cand), fMomentum(cand.P()),
      fTheta(cand.Theta()), fPhi(cand.Phi()) {
  // The covariance was left 0x0 here, so KinFitter read out of bounds.
  fCov.ResizeTo(5, 5);

  // Point of closest approach (PCA) of the line through (X, Y, Z) to the z
  // axis, PCA = (-R sin(phi), R cos(phi), Z). Moving along the line keeps
  // R = -X sin(phi) + Y cos(phi) fixed and reaches the PCA after a transverse
  // step of -(X cos(phi) + Y sin(phi)). (Same result as the former
  // dot-product construction, which divided by |d x z| and gave NaN for a
  // track along z.)
  const double sp = std::sin(fPhi), cp = std::cos(fPhi);
  const double st = std::sin(fTheta), ct = std::cos(fTheta);
  fR = -X * sp + Y * cp;
  fZ = Z;
  if (std::abs(st) > 1e-12) // parallel to z: every point is a PCA
    fZ -= ct / st * (X * cp + Y * sp);
}

KinFitParticle::KinFitParticle() : TLorentzVector() {
  cand.SetXYZM(0, 0, 0, 0);
  fCov.ResizeTo(5, 5);
}

void KinFitParticle::setCovariance(const TMatrixD &cov) {
  // HADES default track parametrization
  // 0 = 1/p
  // 1 = theta
  // 2 = phi
  // 3 = R
  // 4 = z

  // An alternative parametrization
  // 0 = px
  // 1 = py
  // 2 = pz
  // 3 = x (vertex x position)
  // 4 = y (vertex y position)
  // 5 = z (vertex z position)

  if (cov.GetNoElements() == 5 * 5) {
    fCov.ResizeTo(5, 5);
    fCov = cov;
  }

  else if (cov.GetNoElements() == 6 * 6) {
    // TO-DO: implement error-propagation calculation for 6x6 inputs.
  }
}

void KinFitParticle::Print(Option_t *option) const {
  // Print TLorentzVector state first.
  TLorentzVector::Print(option);

  std::cout << "KinFitParticle:"
            << " p=" << fMomentum
            << " theta=" << fTheta
            << " phi=" << fPhi
            << " R=" << fR
            << " Z=" << fZ
            << " q=" << fCharge
            << " pid=" << fPid
            << " trackId=" << fTrackId
            << std::endl;
}

void KinFitParticle::PrintCovariance() const {
  std::cout << "KinFitParticle covariance matrix (parameter order: 1/p, theta, phi, R, Z)"
            << std::endl;
  fCov.Print();
}
