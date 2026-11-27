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

#include "KinFitter.h"
#include "KinFitFourMomentum.h"
#include "KinFitTrackKinematics.h"
#include "KinFitTrackLine.h"
#include "RtypesCore.h"

#include <TVector2.h>

#include <cmath>
#include <numeric>

/**
 * @file KinFitter.cxx
 * @brief Implementation of iterative constrained kinematic fitting.
 */

const size_t cov_dim = 5;

namespace {
/** Resize dst to the shape of src and copy (TMatrixD::operator= requires
 * equal shapes). */
void assign(TMatrixD &dst, const TMatrixD &src) {
  dst.ResizeTo(src);
  dst = src;
}

/**
 * Invert a symmetric positive-definite matrix in place (Cholesky, uses the
 * lower triangle; the result is exactly symmetric). Fails if the matrix is
 * not positive definite or near-singular: det / prod(diag) is the
 * determinant of the corresponding correlation matrix, in (0, 1] and
 * independent of the units of the rows. (Before, TMatrixD::Invert() and a
 * test for det == 0 only, so near-singular systems gave garbage.)
 */
bool invertSPD(TMatrixD &m, double minCorrDet = 1e-12) {
  const int n = m.GetNrows();
  std::vector<double> L(n * n, 0.), Li(n * n, 0.);
  double logCorrDet = 0;
  for (int j = 0; j < n; ++j) {
    const double ajj = m(j, j);
    if (!(ajj > 0) || !std::isfinite(ajj))
      return false;
    double sum = ajj;
    for (int k = 0; k < j; ++k)
      sum -= L[j * n + k] * L[j * n + k];
    if (!(sum > 0))
      return false;
    const double ljj = std::sqrt(sum);
    L[j * n + j] = ljj;
    logCorrDet += std::log(sum / ajj);
    for (int i = j + 1; i < n; ++i) {
      double t = m(i, j);
      for (int k = 0; k < j; ++k)
        t -= L[i * n + k] * L[j * n + k];
      L[i * n + j] = t / ljj;
    }
  }
  if (logCorrDet < std::log(minCorrDet))
    return false;
  // L^-1 (lower triangular), then m^-1 = L^-T L^-1
  for (int i = 0; i < n; ++i) {
    Li[i * n + i] = 1. / L[i * n + i];
    for (int j = 0; j < i; ++j) {
      double t = 0;
      for (int k = j; k < i; ++k)
        t -= L[i * n + k] * Li[k * n + j];
      Li[i * n + j] = t / L[i * n + i];
    }
  }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j <= i; ++j) {
      double t = 0;
      for (int k = i; k < n; ++k)
        t += Li[k * n + i] * Li[k * n + j];
      m(i, j) = m(j, i) = t;
    }
  return true;
}

/** Bring (theta, phi) of a fitted track back to [0, pi] x (-pi, pi]. A
 * negative theta describes the same line as (-theta, phi + pi) with R -> -R;
 * the covariance rows of theta and R change sign accordingly. Indices refer
 * to the (1/p, theta, phi, R, Z) ordering of cov. */
void normaliseAngles(double &theta, double &phi, double &R, TMatrixD &cov) {
  if (theta < 0 || theta > TMath::Pi()) {
    theta = theta < 0 ? -theta : TMath::TwoPi() - theta;
    phi += TMath::Pi();
    R = -R;
    const double sign[5] = {1, -1, 1, -1, 1};
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 5; ++j)
        cov(i, j) *= sign[i] * sign[j];
  }
  phi = TVector2::Phi_mpi_pi(phi);
}
} // namespace

// Main KinFitter constructor.
KinFitter::KinFitter(const std::vector<KinFitParticle> &cands) : fCands(cands) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::KinFitter() -----------------"
              << std::endl;
  }

  // fN is the number of daughters, e.g. (Lambda -> p pi-) gives n=2.
  fN = fCands.size();
  fyDim = fN * cov_dim; // Dimension of full covariance matrix (number of
                        // measured variables x cov_dim)
  fMassVec.clear();

  y.ResizeTo(fyDim, 1);
  x.ResizeTo(1, 1);
  V.ResizeTo(fyDim, fyDim);
  Vx.ResizeTo(1, 1);

  y.Zero();
  x.Zero();
  V.Zero();
  Vx.Zero();

  // set 'y=alpha' measurements
  // and the covariance
  for (int ix = 0; ix < fN; ix++) {
    KinFitParticle cand = fCands[ix];
    y(0 + ix * cov_dim, 0) = 1. / cand.P();
    y(1 + ix * cov_dim, 0) = cand.Theta();
    y(2 + ix * cov_dim, 0) = cand.Phi();
    y(3 + ix * cov_dim, 0) = cand.getR();
    y(4 + ix * cov_dim, 0) = cand.getZ();
    fM.push_back(cand.M());

    TMatrixD covariance = cand.getCovariance();
    V(0 + ix * cov_dim, 0 + ix * cov_dim) = covariance(0, 0);
    V(1 + ix * cov_dim, 0 + ix * cov_dim) = covariance(1, 0);
    V(2 + ix * cov_dim, 0 + ix * cov_dim) = covariance(2, 0);
    V(3 + ix * cov_dim, 0 + ix * cov_dim) = covariance(3, 0);
    V(4 + ix * cov_dim, 0 + ix * cov_dim) = covariance(4, 0);
    V(0 + ix * cov_dim, 1 + ix * cov_dim) = covariance(0, 1);
    V(1 + ix * cov_dim, 1 + ix * cov_dim) = covariance(1, 1);
    V(2 + ix * cov_dim, 1 + ix * cov_dim) = covariance(2, 1);
    V(3 + ix * cov_dim, 1 + ix * cov_dim) = covariance(3, 1);
    V(4 + ix * cov_dim, 1 + ix * cov_dim) = covariance(4, 1);
    V(0 + ix * cov_dim, 2 + ix * cov_dim) = covariance(0, 2);
    V(1 + ix * cov_dim, 2 + ix * cov_dim) = covariance(1, 2);
    V(2 + ix * cov_dim, 2 + ix * cov_dim) = covariance(2, 2);
    V(3 + ix * cov_dim, 2 + ix * cov_dim) = covariance(3, 2);
    V(4 + ix * cov_dim, 2 + ix * cov_dim) = covariance(4, 2);
    V(0 + ix * cov_dim, 3 + ix * cov_dim) = covariance(0, 3);
    V(1 + ix * cov_dim, 3 + ix * cov_dim) = covariance(1, 3);
    V(2 + ix * cov_dim, 3 + ix * cov_dim) = covariance(2, 3);
    V(3 + ix * cov_dim, 3 + ix * cov_dim) = covariance(3, 3);
    V(4 + ix * cov_dim, 3 + ix * cov_dim) = covariance(4, 3);
    V(0 + ix * cov_dim, 4 + ix * cov_dim) = covariance(0, 4);
    V(1 + ix * cov_dim, 4 + ix * cov_dim) = covariance(1, 4);
    V(2 + ix * cov_dim, 4 + ix * cov_dim) = covariance(2, 4);
    V(3 + ix * cov_dim, 4 + ix * cov_dim) = covariance(3, 4);
    V(4 + ix * cov_dim, 4 + ix * cov_dim) = covariance(4, 4);
  }
}

void KinFitter::addMassConstraint(double mass) {
  if (fVerbose > 0) {
    std::cout
        << "--------------- KinFitter::addMassConstraint() -----------------"
        << std::endl;
  }

  // One equation per call; the particles selected with setUseCandInFit()
  // since the previous call belong to this constraint. fMass is left alone:
  // it belongs to MassVtx/MissingMass, which may be combined with this.
  fMassFitPair.push_back(fFlexiParticlesInFit);
  fFlexiParticlesInFit.clear();
  fMassVec.push_back(mass);
  fMassConstraint = true;
}

void KinFitter::addMissingMassConstraint(TLorentzVector lv, double mass) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::addMissingMassConstraint() "
                 "-----------------"
              << std::endl;
  }

  fMass = mass;
  fInit = lv;
  fMissingMassConstraint = true;
}

void KinFitter::addMassVtxConstraint(double mass) {
  if (fVerbose > 0) {
    std::cout
        << "--------------- KinFitter::addMassVtxConstraint() -----------------"
        << std::endl;
  }

  fMass = mass;
  fMassVtxConstraint = true;
}

void KinFitter::add4Constraint(TLorentzVector lv) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::add4Constraint() -----------------"
              << std::endl;
  }

  fInit = lv;

  f4Constraint = true;
}

void KinFitter::add3Constraint(KinFitParticle mother) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::add3Constraint() -----------------"
              << std::endl;
  }

  fyDim = (fN + 1) * cov_dim -
          1; // Dimension of full covariance matrix (number of measured
             // variables x cov_dim). Mother momentum is not measured

  y.ResizeTo(fyDim, 1);
  V.ResizeTo(fyDim, fyDim);

  y.Zero();
  V.Zero();

  fM.clear();

  // Set y to measurements and the covariance, set mass
  for (int ix = 0; ix < fN; ix++) // for daughters
  {
    KinFitParticle cand = fCands[ix];

    y(0 + ix * cov_dim, 0) = 1. / cand.P();
    y(1 + ix * cov_dim, 0) = cand.Theta();
    y(2 + ix * cov_dim, 0) = cand.Phi();
    y(3 + ix * cov_dim, 0) = cand.getR();
    y(4 + ix * cov_dim, 0) = cand.getZ();
    fM.push_back(cand.M());

    // Full 5x5 block (only the diagonal was used before).
    const TMatrixD &covariance = cand.getCovariance();
    for (size_t i = 0; i < cov_dim; ++i)
      for (size_t j = 0; j < cov_dim; ++j)
        V(i + ix * cov_dim, j + ix * cov_dim) = covariance(i, j);
  }

  fMother = mother;

  y(fN * cov_dim, 0) = fMother.Theta();
  y(1 + fN * cov_dim, 0) = fMother.Phi();
  y(2 + fN * cov_dim, 0) = fMother.getR();
  y(3 + fN * cov_dim, 0) = fMother.getZ();
  fM.push_back(fMother.M());

  TMatrixD covariance = fMother.getCovariance();

  if (fVerbose > 0) {
    std::cout << "Decaying Particle Covariance: " << covariance(0, 0) << " "
              << covariance(1, 1) << " " << covariance(2, 2) << " "
              << covariance(3, 3) << " " << covariance(4, 4) << std::endl;
  }

  // Mother (theta, phi, R, Z) block: rows/columns 1..4 of its covariance,
  // with all correlations (only the diagonal and theta-phi were used before).
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      V(i + fN * cov_dim, j + fN * cov_dim) = covariance(i + 1, j + 1);

  if (fVerbose > 0) {
    std::cout << "Fitter" << std::endl;
    std::cout << "Diag elements: " << V(0 + fN * cov_dim, 0 + fN * cov_dim)
              << " " << V(1 + fN * cov_dim, 1 + fN * cov_dim) << std::endl;
    std::cout << "Off diag elements: " << V(0 + fN * cov_dim, 1 + fN * cov_dim)
              << " " << V(1 + fN * cov_dim, 0 + fN * cov_dim) << std::endl;
  }

  f3Constraint = true;
}

void KinFitter::addVertexConstraint() {
  if (fVerbose > 0) {
    std::cout
        << "--------------- KinFitter::addVertexConstraint() -----------------"
        << std::endl;
  }

  fVtxConstraint = true;
}

// Configure fit with one missing particle.
void KinFitter::addMissingParticleConstraint(TLorentzVector lv, double mass) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::addMissingParticleConstraint() "
                 "-----------------"
              << std::endl;
  }

  fInit = lv;
  fMassMissingParticle = mass;

  if (fN == 1) {
    // With one measured track the 4x4 matrix D V D^T has rank <= 3 (only
    // 1/p, theta, phi enter), so the 4C + 3 unknowns formulation cannot be
    // solved. Eliminating p_miss gives the equivalent single equation
    // (P_init - p)^2 = m^2, i.e. the missing-mass constraint (NDF = 1 in both).
    fMass = mass;
    fMissingMassConstraint = true;
    fMissingFromMissingMass = true;
    return;
  }

  x.ResizeTo(3, 1);
  Vx.ResizeTo(3, 3);
  x.Zero();
  Vx.Zero();

  if (!fMissingParticleConstraint)
    fM.push_back(fMassMissingParticle);

  fMissingParticleConstraint = true;
}

TMatrixD KinFitter::calcMissingMom(const TMatrixD &m_iter) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::calcMissingMom() -----------------"
              << std::endl;
  }

  TMatrixD xi(3, 1); // was TMatrix = TMatrixF: float precision

  xi(0, 0) = fInit.Px();
  xi(1, 0) = fInit.Py();
  xi(2, 0) = fInit.Pz();

  for (int q = 0; q < fN; q++) {
    xi(0, 0) -= 1. / m_iter(0 + q * cov_dim, 0) *
                std::sin(m_iter(1 + q * cov_dim, 0)) *
                std::cos(m_iter(2 + q * cov_dim, 0));
    xi(1, 0) -= 1. / m_iter(0 + q * cov_dim, 0) *
                std::sin(m_iter(1 + q * cov_dim, 0)) *
                std::sin(m_iter(2 + q * cov_dim, 0));
    xi(2, 0) -=
        1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0));
  }

  return xi;
}

std::vector<int> KinFitter::massParticles(size_t k) const {
  if (!fMassFitPair[k].empty())
    return fMassFitPair[k];
  const int n = (fNExclusive > 0 && fNExclusive <= fN) ? fNExclusive : fN;
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  return idx;
}

int KinFitter::nUnmeasured() const {
  if (f3Constraint)
    return 1; // 1/p of the mother
  if (fMissingParticleConstraint)
    return 3; // momentum of the missing particle
  return 0;
}

// Evaluate constraint equations.
TMatrixD KinFitter::f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::f_eval() -----------------"
              << std::endl;
  }

  std::vector<KinFitTrackKinematics> k(fN);
  for (int q = 0; q < fN; q++)
    k[q] = KinFitTrackKinematics(m_iter, q * cov_dim, fM[q]);

  TMatrixD d;

  // Invariant mass + vertex constraint
  if (fMassVtxConstraint) {
    d.ResizeTo(2, 1);
    KinFitFourMomentum P;
    for (int q = 0; q < fN; q++)
      P.add(k[q]);
    d(0, 0) = P.m2() - fMass * fMass;
    d(1, 0) = KinFitTrackLine::vertexConstraint(
        KinFitTrackLine(m_iter, 0, k[0]), KinFitTrackLine(m_iter, cov_dim, k[1]));
  }

  // Missing mass constraint: (P_init - sum p_i)^2 - m^2
  if (fMissingMassConstraint) {
    d.ResizeTo(1, 1);
    KinFitFourMomentum P(fInit);
    for (int q = 0; q < fN; q++)
      P.add(k[q], -1.);
    d(0, 0) = P.m2() - fMass * fMass;
  }

  // Vertex constraint
  if (fVtxConstraint) {
    d.ResizeTo(1, 1);
    d(0, 0) = KinFitTrackLine::vertexConstraint(
        KinFitTrackLine(m_iter, 0, k[0]), KinFitTrackLine(m_iter, cov_dim, k[1]));
  }

  // Four-momentum constraints: sum p_i - (mother | initial system - missing)
  if (f3Constraint || f4Constraint || fMissingParticleConstraint) {
    d.ResizeTo(4, 1);
    KinFitFourMomentum P;
    for (int q = 0; q < fN; q++)
      P.add(k[q]);
    if (f3Constraint) { // mother: unmeasured 1/p, measured theta, phi
      const double p = 1. / xi_iter(0, 0);
      const double th = m_iter(0 + fN * cov_dim, 0);
      const double ph = m_iter(1 + fN * cov_dim, 0);
      P.px -= p * std::sin(th) * std::cos(ph);
      P.py -= p * std::sin(th) * std::sin(ph);
      P.pz -= p * std::cos(th);
      P.E -= std::sqrt(p * p + fM[fN] * fM[fN]);
    } else {
      P.px -= fInit.Px();
      P.py -= fInit.Py();
      P.pz -= fInit.Pz();
      P.E -= fInit.E();
      if (fMissingParticleConstraint) {
        const double mx = xi_iter(0, 0), my = xi_iter(1, 0),
                     mz = xi_iter(2, 0);
        P.px += mx;
        P.py += my;
        P.pz += mz;
        P.E += std::sqrt(mx * mx + my * my + mz * mz + fM[fN] * fM[fN]);
      }
    }
    d(0, 0) = P.px;
    d(1, 0) = P.py;
    d(2, 0) = P.pz;
    d(3, 0) = P.E;
  }

  // Mass constraints: one row E^2 - p^2 - m^2 per addMassConstraint() call,
  // appended to the rows of the other constraint (if any).
  if (fMassConstraint) {
    const int row0 = d.GetNrows();
    d.ResizeTo(row0 + (int)fMassVec.size(), 1);
    for (size_t c = 0; c < fMassVec.size(); ++c) {
      KinFitFourMomentum P;
      for (int q : massParticles(c))
        P.add(k[q]);
      d(row0 + c, 0) = P.m2() - fMassVec[c] * fMassVec[c];
    }
  }

  return d;
}

// Jacobian (derivative of constraint equations with respect to measured
// variables)
TMatrixD KinFitter::Feta_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::Feta_eval() -----------------"
              << std::endl;
  }

  std::vector<KinFitTrackKinematics> k(fN);
  for (int q = 0; q < fN; q++)
    k[q] = KinFitTrackKinematics(m_iter, q * cov_dim, fM[q]);

  TMatrixD H;
  double dInvP, dTheta, dPhi;

  if (fMassVtxConstraint) {
    H.ResizeTo(2, fyDim);
    H.Zero();
    KinFitFourMomentum P;
    for (int q = 0; q < fN; q++)
      P.add(k[q]);
    for (int q = 0; q < fN; q++) {
      k[q].massDerivatives(P, 1., dInvP, dTheta, dPhi);
      H(0, 0 + q * cov_dim) = dInvP;
      H(0, 1 + q * cov_dim) = dTheta;
      H(0, 2 + q * cov_dim) = dPhi;
    }
    KinFitTrackLine::vertexJacobian(H, 1, KinFitTrackLine(m_iter, 0, k[0]), 0,
                                    KinFitTrackLine(m_iter, cov_dim, k[1]),
                                    cov_dim);
  }

  if (fMissingMassConstraint) {
    H.ResizeTo(1, fyDim);
    H.Zero();
    KinFitFourMomentum P(fInit);
    for (int q = 0; q < fN; q++)
      P.add(k[q], -1.);
    for (int q = 0; q < fN; q++) {
      k[q].massDerivatives(P, -1., dInvP, dTheta, dPhi);
      H(0, 0 + q * cov_dim) = dInvP;
      H(0, 1 + q * cov_dim) = dTheta;
      H(0, 2 + q * cov_dim) = dPhi;
    }
  }

  if (fVtxConstraint) {
    H.ResizeTo(1, fyDim);
    H.Zero();
    KinFitTrackLine::vertexJacobian(H, 0, KinFitTrackLine(m_iter, 0, k[0]), 0,
                                    KinFitTrackLine(m_iter, cov_dim, k[1]),
                                    cov_dim);
  }

  if (f3Constraint || f4Constraint || fMissingParticleConstraint) {
    H.ResizeTo(4, fyDim);
    H.Zero();
    for (int q = 0; q < fN; q++)
      k[q].momentumJacobian(H, q * cov_dim);
    if (f3Constraint) { // mother theta, phi enter with -1/xi
      const double p = 1. / xi_iter(0, 0);
      const double th = m_iter(0 + fN * cov_dim, 0);
      const double ph = m_iter(1 + fN * cov_dim, 0);
      const double st = std::sin(th), ct = std::cos(th);
      const double sp = std::sin(ph), cp = std::cos(ph);
      H(0, fN * cov_dim) = -p * ct * cp;
      H(1, fN * cov_dim) = -p * ct * sp;
      H(2, fN * cov_dim) = p * st;
      H(0, 1 + fN * cov_dim) = p * st * sp;
      H(1, 1 + fN * cov_dim) = -p * st * cp;
    }
  }

  // Mass constraints (see f_eval).
  if (fMassConstraint) {
    const int row0 = H.GetNrows();
    H.ResizeTo(row0 + (int)fMassVec.size(), fyDim);
    for (size_t c = 0; c < fMassVec.size(); ++c) {
      const int row = row0 + c;
      for (int col = 0; col < fyDim; ++col)
        H(row, col) = 0.;
      const std::vector<int> idx = massParticles(c);
      KinFitFourMomentum P;
      for (int q : idx)
        P.add(k[q]);
      for (int q : idx) {
        k[q].massDerivatives(P, 1., dInvP, dTheta, dPhi);
        H(row, 0 + q * cov_dim) = dInvP;
        H(row, 1 + q * cov_dim) = dTheta;
        H(row, 2 + q * cov_dim) = dPhi;
      }
    }
  }

  return H;
}

// Jacobian (derivative of constraint equations with respect to unmeasured
// variables)
TMatrixD KinFitter::Fxi_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::Fxi_eval() -----------------"
              << std::endl;
  }

  if (fVerbose > 4) {
    // m_iter.Print();

    // xi_iter.Print();
  }

  // std::cout << "CovDim: " << cov_dim << std::endl;
  // std::cout << "fN: " << fN << std::endl;

  TMatrixD H;

  if (f3Constraint) {
    H.ResizeTo(4, 1);
    H.Zero();

    // d(1/p)
    H(0, 0) = 1. / std::pow(xi_iter(0, 0), 2) *
              std::sin(m_iter(0 + fN * cov_dim, 0)) *
              std::cos(m_iter(1 + fN * cov_dim, 0));
    H(1, 0) = 1. / std::pow(xi_iter(0, 0), 2) *
              std::sin(m_iter(0 + fN * cov_dim, 0)) *
              std::sin(m_iter(1 + fN * cov_dim, 0));
    H(2, 0) =
        1. / std::pow(xi_iter(0, 0), 2) * std::cos(m_iter(0 + fN * cov_dim, 0));
    H(3, 0) = 1. / std::pow(xi_iter(0, 0), 3) * 1. /
              sqrt(pow(1. / xi_iter(0, 0), 2) + std::pow(fM[fN], 2));
  }

  if (fMissingParticleConstraint) {
    // std::cout << "xi_iter: " << xi_iter(0, 0) << " " << xi_iter(1, 0) << " "
    // << xi_iter(2, 0) << std::endl;

    H.ResizeTo(4, 3);
    H.Zero();

    H(0, 0) = 1;
    H(1, 1) = 1;
    H(2, 2) = 1;

    H(3, 0) = xi_iter(0, 0) /
              sqrt(pow(xi_iter(0, 0), 2) + std::pow(xi_iter(1, 0), 2) +
                   std::pow(xi_iter(2, 0), 2) + std::pow(fM[fN], 2));
    H(3, 1) = xi_iter(1, 0) /
              sqrt(pow(xi_iter(0, 0), 2) + std::pow(xi_iter(1, 0), 2) +
                   std::pow(xi_iter(2, 0), 2) + std::pow(fM[fN], 2));
    H(3, 2) = xi_iter(2, 0) /
              sqrt(pow(xi_iter(0, 0), 2) + std::pow(xi_iter(1, 0), 2) +
                   std::pow(xi_iter(2, 0), 2) + std::pow(fM[fN], 2));
  }
  // Mass constraints do not depend on the unmeasured parameters: zero rows.
  if (fMassConstraint && H.GetNoElements() > 0) {
    const int row0 = H.GetNrows();
    H.ResizeTo(row0 + (int)fMassVec.size(), H.GetNcols());
    for (int r = row0; r < H.GetNrows(); ++r)
      for (int c = 0; c < H.GetNcols(); ++c)
        H(r, c) = 0.;
  }

  return H;
}

Bool_t KinFitter::fit() {

  if (fVerbose > 0) {
    std::cout << " ----------- KinFitter::fit() -----------" << std::endl;
    std::cout << "Vertex constraint set: " << fVtxConstraint << std::endl;
    std::cout << "3C set: " << f3Constraint << std::endl;
    std::cout << "4C set: " << f4Constraint << std::endl;
    std::cout << "Momentum constraint set: " << fMissingParticleConstraint
              << std::endl;
    std::cout << "Mass constraint set: " << fMassConstraint << std::endl;
    std::cout << "Missing Mass constraint set: " << fMissingMassConstraint
              << std::endl;
    std::cout << "Mass + Vertex constraint set: " << fMassVtxConstraint
              << std::endl;
    std::cout << "" << std::endl;
  }

  fConverged = false;
  fIteration = 0;

  // --- Configuration checks (used to abort() or silently overwrite d) ---
  const int nBase = int(f3Constraint) + int(f4Constraint) +
                    int(fVtxConstraint) + int(fMissingParticleConstraint) +
                    int(fMissingMassConstraint) + int(fMassVtxConstraint);
  if (nBase == 0 && !fMassConstraint) {
    std::cerr << "KinFitter::fit(): no constraint added" << std::endl;
    return kFALSE;
  }
  if (nBase > 1) {
    std::cerr << "KinFitter::fit(): 3C, 4C, Vertex, MassVtx, MissingMass and "
                 "MissingParticle are exclusive (mass constraints can be "
                 "added to any of them)"
              << std::endl;
    return kFALSE;
  }
  if ((fVtxConstraint || fMassVtxConstraint) && fN < 2) {
    std::cerr << "KinFitter::fit(): vertex constraint needs two particles"
              << std::endl;
    return kFALSE;
  }
  for (size_t k = 0; k < fMassVec.size(); ++k)
    for (int q : massParticles(k))
      if (q < 0 || q >= fN) {
        std::cerr << "KinFitter::fit(): particle index " << q
                  << " in mass constraint out of range" << std::endl;
        return kFALSE;
      }

  const bool hasXi = nUnmeasured() > 0;
  const TMatrixD alpha0(y); // measured values, kept unchanged
  const TMatrixD V0(V);     // their covariance
  TMatrixD alpha(alpha0);
  TMatrixD xi0(x);
  if (f3Constraint)
    xi0(0, 0) = 1. / fMother.P();
  if (fMissingParticleConstraint)
    assign(xi0, calcMissingMom(alpha0));
  TMatrixD xi(xi0), neu_xi(xi0);
  double chi2 = 1e6;

  TMatrixD D = Feta_eval(alpha, xi);
  TMatrixD d = f_eval(alpha, xi);
  fNdf = d.GetNrows() - nUnmeasured();
  if (fNdf < 1) {
    std::cerr << "KinFitter::fit(): NDF = " << fNdf
              << " (fewer constraints than unknowns)" << std::endl;
    return kFALSE;
  }
  TMatrixD D_xi, DT_xi;
  if (hasXi)
    assign(D_xi, Fxi_eval(alpha, xi));

  TMatrixD DT(TMatrixD::kTransposed, D);
  TMatrixD V0DT, S, VD, VDD;

  for (int q = 0; q < fNumIterations; q++) {
    if (fVerbose > 0) {
      std::cout << "Number of iterations: " << q << std::endl;
    }
    fIteration = q;

    // Linearised constraint r = f(alpha) + D (alpha0 - alpha)
    const TMatrixD r = d + D * (alpha0 - alpha);
    DT.Transpose(D);
    assign(V0DT, V0 * DT);
    assign(S, D * V0DT); // D V0 D^T
    assign(VD, S);
    if (!invertSPD(VD)) {
      if (q == 0)
        return kFALSE; // protect against failed matrix inversion
      break;           // keep the last good iteration
    }

    TMatrixD lambda; // Lagrange multipliers
    if (hasXi) {
      DT_xi.ResizeTo(D_xi.GetNcols(), D_xi.GetNrows());
      DT_xi.Transpose(D_xi);
      assign(VDD, DT_xi * VD * D_xi);
      if (!invertSPD(VDD)) {
        if (q == 0)
          return kFALSE;
        break;
      }
      neu_xi = xi - VDD * DT_xi * VD * r;
      assign(lambda, VD * (r + D_xi * (neu_xi - xi)));
    } else {
      assign(lambda, VD * r);
    }

    const TMatrixD neu_alpha = alpha0 - V0DT * lambda;

    // chi2 = da^T V0^-1 da + 2 lambda^T f(alpha_new, xi_new). With
    // da = V0 D^T lambda the first term is lambda^T (D V0 D^T) lambda, so V0
    // never has to be inverted. f(alpha_new) is the next iteration's d.
    const TMatrixD d_new = f_eval(neu_alpha, neu_xi);
    const TMatrixD lambdaT(TMatrixD::kTransposed, lambda);
    const double chi2_new = (lambdaT * S * lambda)(0, 0) +
                            2. * (lambdaT * d_new)(0, 0);

    // Convergence. The parameter change is measured in units of the
    // measurement uncertainty; it used to be divided by the measured value,
    // giving NaN (never converged) for any parameter that is exactly 0.
    double dNorm = 0;
    for (int i = 0; i < d.GetNrows(); i++)
      dNorm += d(i, 0) * d(i, 0);
    dNorm = std::sqrt(dNorm);
    double alphaNorm = 0;
    for (int i = 0; i < fyDim; i++) {
      const double var = V0(i, i);
      if (var > 0) {
        const double da = alpha(i, 0) - neu_alpha(i, 0);
        alphaNorm += da * da / var;
      }
    }
    if (hasXi) {
      for (int i = 0; i < xi.GetNrows(); i++) {
        const double scale = std::abs(xi0(i, 0));
        if (scale > 0) {
          const double dx = (neu_xi(i, 0) - xi(i, 0)) / scale;
          alphaNorm += dx * dx;
        }
      }
    }
    alphaNorm = std::sqrt(alphaNorm);

    const bool converged = std::abs(chi2 - chi2_new) < fConvergenceCriterionChi2 &&
                           dNorm < fConvergenceCriterionD &&
                           alphaNorm < fConvergenceCriterionAlpha;
    chi2 = chi2_new;
    alpha = neu_alpha;
    if (hasXi)
      xi = neu_xi;
    if (converged) {
      fConverged = true;
      break;
    }

    D = Feta_eval(alpha, xi);
    if (hasXi)
      D_xi = Fxi_eval(alpha, xi);
    d = d_new;
  }

  if (fVerbose > 0) {
    std::cout << " Exit iterations" << std::endl;
  }

  // --- Covariance of the fitted parameters ---
  // Linearised at the final point. Before, D had already been advanced to the
  // next point while DT and VD were from the previous one, giving an
  // inconsistent (non-symmetric) matrix whenever the iteration limit was hit.
  assign(D, Feta_eval(alpha, xi));
  DT.ResizeTo(D.GetNcols(), D.GetNrows());
  DT.Transpose(D);
  assign(V0DT, V0 * DT);
  assign(VD, D * V0DT);
  if (!invertSPD(VD))
    return kFALSE;
  const TMatrixD DV0(TMatrixD::kTransposed, V0DT); // D V0 (V0 symmetric)
  assign(fVFit, V0 - V0DT * VD * DV0);
  if (hasXi) {
    assign(D_xi, Fxi_eval(alpha, xi));
    DT_xi.ResizeTo(D_xi.GetNcols(), D_xi.GetNrows());
    DT_xi.Transpose(D_xi);
    TMatrixD VxFit = DT_xi * VD * D_xi;
    if (!invertSPD(VxFit))
      return kFALSE;
    const TMatrixD G = V0DT * VD * D_xi; // fyDim x n_xi
    const TMatrixD GT(TMatrixD::kTransposed, G);
    fVFit += G * VxFit * GT;
    assign(Vx, VxFit);
    assign(fVxy, -1. * (VxFit * GT)); // cov(xi, alpha)
  }

  assign(fYFit, alpha);
  assign(x, xi);
  if (fMissingParticleConstraint) {
    fMissDaughter.SetXYZM(x(0, 0), x(1, 0), x(2, 0), fM[fN]);
    assign(fVMiss, Vx);
  }
  if (fMissingFromMissingMass) {
    // p_miss = P_init - p of the single track; its covariance follows from
    // the fitted (1/p, theta, phi) block of that track.
    const double p = 1. / alpha(0, 0), th = alpha(1, 0), ph = alpha(2, 0);
    const double st = std::sin(th), ct = std::cos(th);
    const double sp = std::sin(ph), cp = std::cos(ph);
    const TVector3 pTrack(p * st * cp, p * st * sp, p * ct);
    fMissDaughter.SetVectM(fInit.Vect() - pTrack, fMass);
    TMatrixD J(3, 3); // d(p_miss)/d(1/p, theta, phi) = -d(p_track)/d(...)
    J(0, 0) = p * p * st * cp;
    J(1, 0) = p * p * st * sp;
    J(2, 0) = p * p * ct;
    J(0, 1) = -p * ct * cp;
    J(1, 1) = -p * ct * sp;
    J(2, 1) = p * st;
    J(0, 2) = p * st * sp;
    J(1, 2) = -p * st * cp;
    J(2, 2) = 0;
    const TMatrixD Vtrack = fVFit.GetSub(0, 2, 0, 2);
    assign(fVMiss, J * Vtrack * TMatrixD(TMatrixD::kTransposed, J));
  }
  fChi2 = chi2;
  fProb = TMath::Prob(chi2, fNdf);

  // -----------------------------------------
  // Pulls
  // -----------------------------------------
  if (fVerbose > 1) {
    std::cout << " Calculate pulls" << std::endl;
  }
  fPull.assign(fyDim, -10000.);
  for (int b = 0; b < fyDim; b++) {
    const double num = alpha0(b, 0) - alpha(b, 0);
    const double dem = V0(b, b) - fVFit(b, b);
    if (dem > 0)
      fPull[b] = num / std::sqrt(dem);
  }

  updateDaughters();
  if (f3Constraint)
    updateMother();

  if (fNumIterations == 1)
    return kTRUE; // single linearised step: no convergence test
  return fConverged;
}

void KinFitter::updateDaughters() {
  if (fVerbose > 1) {
    std::cout << " Update daughters" << std::endl;
  }
  for (int val = 0; val < fN; ++val) {
    KinFitParticle &cand = fCands[val];
    const double p = 1. / fYFit(0 + val * cov_dim, 0);
    double theta = fYFit(1 + val * cov_dim, 0);
    double phi = fYFit(2 + val * cov_dim, 0);
    double R = fYFit(3 + val * cov_dim, 0);

    // Full 5x5 block of the fitted covariance.
    TMatrixD cov(cov_dim, cov_dim);
    for (size_t i = 0; i < cov_dim; ++i)
      for (size_t j = 0; j < cov_dim; ++j)
        cov(i, j) = fVFit(i + val * cov_dim, j + val * cov_dim);
    // The fit works with unwrapped angles (pulls stay correct); only the
    // returned track is brought back to the standard ranges.
    normaliseAngles(theta, phi, R, cov);

    cand.SetXYZM(p * std::sin(theta) * std::cos(phi),
                 p * std::sin(theta) * std::sin(phi), p * std::cos(theta),
                 fM[val]);
    cand.setMomentum(p);
    cand.setThetaRad(theta);
    cand.setPhiRad(phi);
    cand.setR(R);
    cand.setZ(fYFit(4 + val * cov_dim, 0));
    cand.setCovariance(cov);
  }
}

void KinFitter::updateMother() {
  if (fVerbose > 1) {
    std::cout << " Update mother" << std::endl;
  }

  KinFitParticle &mother = fMother;
  const double p = 1. / x(0, 0);
  double theta = fYFit(0 + fN * cov_dim, 0);
  double phi = fYFit(1 + fN * cov_dim, 0);
  double R = fYFit(2 + fN * cov_dim, 0);

  // (1/p, theta, phi, R, Z): 1/p is the unmeasured parameter, the rest the
  // measured mother block, including their correlations.
  TMatrixD cov(cov_dim, cov_dim);
  cov(0, 0) = Vx(0, 0);
  for (int i = 0; i < 4; ++i) {
    cov(0, i + 1) = cov(i + 1, 0) = fVxy(0, i + fN * cov_dim);
    for (int j = 0; j < 4; ++j)
      cov(i + 1, j + 1) = fVFit(i + fN * cov_dim, j + fN * cov_dim);
  }
  normaliseAngles(theta, phi, R, cov);

  mother.SetXYZM(p * std::sin(theta) * std::cos(phi),
                 p * std::sin(theta) * std::sin(phi), p * std::cos(theta),
                 fM[fN]);
  mother.setMomentum(p);
  mother.setThetaRad(theta);
  mother.setPhiRad(phi);
  mother.setR(R);
  mother.setZ(fYFit(3 + fN * cov_dim, 0));
  mother.setCovariance(cov);
}
