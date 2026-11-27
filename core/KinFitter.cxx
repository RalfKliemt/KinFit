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
#include "KinFitTrajectory.h"
#include "KinFitVertexFinderField.h"
#include "RtypesCore.h"

#include <TVector2.h>

#include <cmath>
#include <limits>
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
KinFitter::KinFitter(const std::vector<KinFitParticle> &cands)
    : fCands(cands), fCandsIn(cands) {
  // fN is the number of daughters, e.g. (Lambda -> p pi-) gives n=2.
  fN = fCands.size();
  fyDim = fN * cov_dim; // 5 parameters per particle

  y.ResizeTo(fyDim, 1);
  x.ResizeTo(1, 1);
  V.ResizeTo(fyDim, fyDim);
  Vx.ResizeTo(1, 1);
  y.Zero();
  x.Zero();
  V.Zero();
  Vx.Zero();

  for (int ix = 0; ix < fN; ix++)
    fillParticle(ix);
}

void KinFitter::fillParticle(int ix) {
  const KinFitParticle &cand = fCands[ix];
  const int o = ix * cov_dim;
  const bool neutral = cand.isNeutralHit();
  if (neutral) { // (E, X, Y, Z, unused)
    const TVector3 hit = cand.getHit();
    y(o + 0, 0) = cand.E();
    y(o + 1, 0) = hit.X();
    y(o + 2, 0) = hit.Y();
    y(o + 3, 0) = hit.Z();
    y(o + 4, 0) = 0.;
  } else { // (1/p, theta, phi, R, Z)
    y(o + 0, 0) = 1. / cand.P();
    y(o + 1, 0) = cand.Theta();
    y(o + 2, 0) = cand.Phi();
    y(o + 3, 0) = cand.getR();
    y(o + 4, 0) = cand.getZ();
  }
  const double m2 = cand.M2();
  fM.push_back(m2 > 0 ? std::sqrt(m2) : 0.);
  fIsNeutral.push_back(neutral);
  fCharge.push_back(neutral ? 0. : cand.getCharge());

  // Production point of a neutral hit when no vertex is fitted: the
  // particle's vertex if set, else the origin.
  TVector3 origin;
  if (cand.getVertexX() > -1e8 && cand.getVertexY() > -1e8 &&
      cand.getVertexZ() > -1e8)
    origin.SetXYZ(cand.getVertexX(), cand.getVertexY(), cand.getVertexZ());
  fOrigin.push_back(origin);

  // Reference axis for the two vertex residual directions: the coordinate
  // axis least aligned with the measured direction (fixed during the fit).
  const TVector3 dir = neutral ? TVector3(0, 0, 1) : cand.Vect().Unit();
  const double ax = std::abs(dir.X()), ay = std::abs(dir.Y()),
               az = std::abs(dir.Z());
  fVtxAxis.push_back(ax <= ay && ax <= az   ? TVector3(1, 0, 0)
                     : ay <= ax && ay <= az ? TVector3(0, 1, 0)
                                            : TVector3(0, 0, 1));

  const TMatrixD &covariance = cand.getCovariance();
  for (size_t i = 0; i < cov_dim; ++i)
    for (size_t j = 0; j < cov_dim; ++j)
      V(i + o, j + o) = covariance(i, j);
}

void KinFitter::addMassConstraint(double mass) {
  if (fVerbose > 0) {
    std::cout
        << "--------------- KinFitter::addMassConstraint() -----------------"
        << std::endl;
  }

  // One equation per call; the particles selected with setUseCandInFit()
  // since the previous call belong to this constraint.
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
  addCommonVertexConstraint();
  addMassConstraint(mass);
}

void KinFitter::addVertexConstraint() { addCommonVertexConstraint(); }

void KinFitter::addCommonVertexConstraint() {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::addCommonVertexConstraint() "
                 "-----------------"
              << std::endl;
  }
  fCommonVertex = true;
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

  // Daughters keep their 5 parameters; the mother adds (theta, phi, R, Z),
  // its 1/p is unmeasured.
  fyDim = (fN + 1) * cov_dim - 1;
  y.ResizeTo(fyDim, 1);
  V.ResizeTo(fyDim, fyDim);
  for (int i = fN * cov_dim; i < fyDim; ++i) {
    y(i, 0) = 0.;
    for (int j = 0; j < fyDim; ++j)
      V(i, j) = V(j, i) = 0.;
  }

  fMother = mother;
  y(fN * cov_dim, 0) = fMother.Theta();
  y(1 + fN * cov_dim, 0) = fMother.Phi();
  y(2 + fN * cov_dim, 0) = fMother.getR();
  y(3 + fN * cov_dim, 0) = fMother.getZ();
  fM.resize(fN);
  fM.push_back(fMother.M());

  // Mother (theta, phi, R, Z) block: rows/columns 1..4 of its covariance.
  const TMatrixD &covariance = fMother.getCovariance();
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      V(i + fN * cov_dim, j + fN * cov_dim) = covariance(i + 1, j + 1);

  f3Constraint = true;
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
    // With one measured particle the 4x4 matrix D V D^T has rank <= 3, so
    // the 4C + 3 unknowns formulation cannot be solved. Eliminating p_miss
    // gives the equivalent single equation (P_init - p)^2 = m^2, i.e. the
    // missing-mass constraint (NDF = 1 in both forms).
    fMass = mass;
    fMissingMassConstraint = true;
    fMissingFromMissingMass = true;
    return;
  }

  fM.resize(fN);
  fM.push_back(fMassMissingParticle);
  fMissingParticleConstraint = true;
}

std::vector<int> KinFitter::massParticles(size_t k) const {
  if (!fMassFitPair[k].empty())
    return fMassFitPair[k];
  const int n = (fNExclusive > 0 && fNExclusive <= fN) ? fNExclusive : fN;
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  return idx;
}

int KinFitter::nBaseUnmeasured() const {
  if (f3Constraint)
    return 1; // 1/p of the mother
  if (fMissingParticleConstraint)
    return 3; // momentum of the missing particle
  return 0;
}

int KinFitter::nUnmeasured() const {
  return nBaseUnmeasured() + (fCommonVertex ? 3 : 0);
}

int KinFitter::nCharged() const {
  int n = 0;
  for (int q = 0; q < fN; ++q)
    n += !fIsNeutral[q];
  return n;
}

bool KinFitter::transportToVertex() const {
  if (!fCommonVertex || !fField)
    return false;
  // a uniform zero field is a straight line: nothing to transport
  return !(fField->isUniform() && fField->getField(0., 0., 0.).Mag2() == 0);
}

TVector3 KinFitter::momentumAtVertex(const double eta[5], int q,
                                     const TVector3 &v, double *sInOut) const {
  const KinFitTrajectory t = KinFitTrajectory::fromTrackParameters(
      eta[0], eta[1], eta[2], eta[3], eta[4], fCharge[q], fField.get());
  KinFitTrajectory::Vec pos, dir;
  // warm start from a nearby solution (finite differences); NaN = none
  const double s = t.closestApproach({v.X(), v.Y(), v.Z()}, pos, dir,
                                     sInOut ? *sInOut
                                            : std::numeric_limits<double>::quiet_NaN());
  if (sInOut)
    *sInOut = s;
  const double p = 1. / eta[0];
  return TVector3(p * dir[0], p * dir[1], p * dir[2]);
}

std::vector<KinFitTrackKinematics>
KinFitter::kinematics(const TMatrixD &m, const TMatrixD &xi,
                      bool derivatives) const {
  const bool useVertex = fCommonVertex && xi.GetNrows() >= vertexOffset() + 3;
  TVector3 v;
  if (useVertex)
    v.SetXYZ(xi(vertexOffset(), 0), xi(vertexOffset() + 1, 0),
             xi(vertexOffset() + 2, 0));
  const bool transport = useVertex && transportToVertex();
  std::vector<KinFitTrackKinematics> k(fN);
  for (int q = 0; q < fN; q++) {
    if (fIsNeutral[q]) {
      k[q] = KinFitTrackKinematics::neutralHit(m, q * cov_dim, fM[q],
                                               useVertex ? v : fOrigin[q]);
      continue;
    }
    k[q] = KinFitTrackKinematics(m, q * cov_dim, fM[q]);
    if (!transport || fCharge[q] == 0)
      continue;
    // In a field the momentum direction changes along the track: momentum
    // conservation holds at the vertex, so use the momentum at the closest
    // approach to the vertex (|p| and E unchanged). Derivatives by central
    // differences.
    double eta[5];
    for (int j = 0; j < 5; ++j)
      eta[j] = m(j + q * cov_dim, 0);
    double s0 = std::numeric_limits<double>::quiet_NaN();
    const TVector3 p = momentumAtVertex(eta, q, v, &s0);
    k[q].px = p.X();
    k[q].py = p.Y();
    k[q].pz = p.Z();
    if (!derivatives)
      continue;
    for (int j = 0; j < 5; ++j) {
      const double h = j == 0   ? 1e-6 * std::abs(eta[0])
                       : j <= 2 ? 1e-6
                                : 1e-4;
      double ep[5], em[5];
      std::copy(eta, eta + 5, ep);
      std::copy(eta, eta + 5, em);
      ep[j] += h;
      em[j] -= h;
      double sp = s0, sm = s0;
      const TVector3 dp = (momentumAtVertex(ep, q, v, &sp) -
                           momentumAtVertex(em, q, v, &sm)) *
                          (0.5 / h);
      for (int i = 0; i < 3; ++i)
        k[q].dEta[i][j] = dp[i];
    }
    for (int j = 0; j < 3; ++j) {
      TVector3 vp(v), vm(v);
      vp[j] += 1e-4;
      vm[j] -= 1e-4;
      double sp = s0, sm = s0;
      const TVector3 dp = (momentumAtVertex(eta, q, vp, &sp) -
                           momentumAtVertex(eta, q, vm, &sm)) *
                          (0.5 / 1e-4);
      for (int i = 0; i < 3; ++i)
        k[q].dVtx[i][j] = dp[i];
    }
  }
  return k;
}

void KinFitter::vertexResiduals(const double eta[5], int q, const TVector3 &v,
                                double res[2], double *sInOut) const {
  const KinFitTrajectory t = KinFitTrajectory::fromTrackParameters(
      eta[0], eta[1], eta[2], eta[3], eta[4], fCharge[q], fField.get());
  KinFitTrajectory::Vec pos, dir;
  const KinFitTrajectory::Vec vv = {v.X(), v.Y(), v.Z()};
  // warm start from a nearby solution (finite differences); NaN = none
  const double s =
      t.closestApproach(vv, pos, dir,
                        sInOut ? *sInOut
                               : std::numeric_limits<double>::quiet_NaN());
  if (sInOut)
    *sInOut = s;
  // Distance vector at the point of closest approach (perpendicular to the
  // track), in two directions e1 = dir x a / |.|, e2 = dir x e1.
  const TVector3 &a = fVtxAxis[q];
  double e1[3] = {dir[1] * a.Z() - dir[2] * a.Y(), dir[2] * a.X() - dir[0] * a.Z(),
                  dir[0] * a.Y() - dir[1] * a.X()};
  const double n1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  for (double &c : e1)
    c /= n1;
  const double e2[3] = {dir[1] * e1[2] - dir[2] * e1[1],
                        dir[2] * e1[0] - dir[0] * e1[2],
                        dir[0] * e1[1] - dir[1] * e1[0]};
  const double d[3] = {pos[0] - vv[0], pos[1] - vv[1], pos[2] - vv[2]};
  res[0] = d[0] * e1[0] + d[1] * e1[1] + d[2] * e1[2];
  res[1] = d[0] * e2[0] + d[1] * e2[1] + d[2] * e2[2];
}

bool KinFitter::startValues(TMatrixD &xi0) const {
  xi0.ResizeTo(std::max(nUnmeasured(), 1), 1);
  xi0.Zero();
  if (fCommonVertex) {
    // Unweighted vertex of all charged tracks in the field
    const KinFitVertexFinderField finder(fCandsIn, fField.get());
    TVector3 v;
    if (finder.isValid())
      v = finder.getVertex();
    else if (nCharged() == 1) { // one track: its reference point
      for (int q = 0; q < fN; ++q)
        if (!fIsNeutral[q])
          v.SetXYZ(-y(3 + q * cov_dim, 0) * std::sin(y(2 + q * cov_dim, 0)),
                   y(3 + q * cov_dim, 0) * std::cos(y(2 + q * cov_dim, 0)),
                   y(4 + q * cov_dim, 0));
    } else if (nCharged() > 1)
      return false; // e.g. parallel tracks
    for (int j = 0; j < 3; ++j)
      xi0(vertexOffset() + j, 0) = v[j];
  }
  if (f3Constraint)
    xi0(0, 0) = 1. / fMother.P();
  if (fMissingParticleConstraint) {
    const std::vector<KinFitTrackKinematics> k = kinematics(y, xi0);
    xi0(0, 0) = fInit.Px();
    xi0(1, 0) = fInit.Py();
    xi0(2, 0) = fInit.Pz();
    for (int q = 0; q < fN; q++) {
      xi0(0, 0) -= k[q].px;
      xi0(1, 0) -= k[q].py;
      xi0(2, 0) -= k[q].pz;
    }
  }
  return true;
}

bool KinFitter::getStartValues(TMatrixD &yOut, TMatrixD &xiOut) const {
  yOut.ResizeTo(y);
  yOut = y;
  return startValues(xiOut);
}

namespace {
/** Number of rows of the kinematic constraint. */
int nBaseRows(bool missingMass, bool fourRows) {
  return missingMass ? 1 : (fourRows ? 4 : 0);
}
// Finite-difference steps for the vertex Jacobians: 1/p (relative), angles
// [rad], lengths [mm]; neutral hits never enter the vertex rows.
constexpr double kStepInvPRel = 1e-6, kStepAngle = 1e-6, kStepLength = 1e-4;
} // namespace

// Evaluate constraint equations: [kinematic rows][mass rows][vertex rows].
TMatrixD KinFitter::f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::f_eval() -----------------"
              << std::endl;
  }

  const std::vector<KinFitTrackKinematics> k =
      kinematics(m_iter, xi_iter, false);
  const bool fourRows = f3Constraint || f4Constraint || fMissingParticleConstraint;
  const int nVtxRows = fCommonVertex ? 2 * nCharged() : 0;
  const int nRows = nBaseRows(fMissingMassConstraint, fourRows) +
                    (int)fMassVec.size() + nVtxRows;
  TMatrixD d(std::max(nRows, 1), 1);
  int row = 0;

  // Missing mass constraint: (P_init - sum p_i)^2 - m^2
  if (fMissingMassConstraint) {
    KinFitFourMomentum P(fInit);
    for (int q = 0; q < fN; q++)
      P.add(k[q], -1.);
    d(row++, 0) = P.m2() - fMass * fMass;
  }

  // Four-momentum constraints: sum p_i - (mother | initial system - missing)
  if (fourRows) {
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
    d(row++, 0) = P.px;
    d(row++, 0) = P.py;
    d(row++, 0) = P.pz;
    d(row++, 0) = P.E;
  }

  // Mass constraints: one row E^2 - p^2 - m^2 per addMassConstraint() call.
  for (size_t c = 0; c < fMassVec.size(); ++c) {
    KinFitFourMomentum P;
    for (int q : massParticles(c))
      P.add(k[q]);
    d(row++, 0) = P.m2() - fMassVec[c] * fMassVec[c];
  }

  // Common vertex: two residuals per charged track
  if (fCommonVertex) {
    const int vo = vertexOffset();
    const TVector3 v(xi_iter(vo, 0), xi_iter(vo + 1, 0), xi_iter(vo + 2, 0));
    for (int q = 0; q < fN; ++q) {
      if (fIsNeutral[q])
        continue;
      double eta[5], res[2];
      for (int j = 0; j < 5; ++j)
        eta[j] = m_iter(j + q * cov_dim, 0);
      vertexResiduals(eta, q, v, res);
      d(row++, 0) = res[0];
      d(row++, 0) = res[1];
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

  const std::vector<KinFitTrackKinematics> k = kinematics(m_iter, xi_iter);
  const bool fourRows = f3Constraint || f4Constraint || fMissingParticleConstraint;
  const int nVtxRows = fCommonVertex ? 2 * nCharged() : 0;
  const int nRows = nBaseRows(fMissingMassConstraint, fourRows) +
                    (int)fMassVec.size() + nVtxRows;
  TMatrixD H(std::max(nRows, 1), fyDim);
  int row = 0;
  double de[5], dv[3];

  if (fMissingMassConstraint) {
    KinFitFourMomentum P(fInit);
    for (int q = 0; q < fN; q++)
      P.add(k[q], -1.);
    for (int q = 0; q < fN; q++) {
      k[q].massDerivatives(P, -1., de, dv);
      for (int j = 0; j < 5; ++j)
        H(row, j + q * cov_dim) = de[j];
    }
    ++row;
  }

  if (fourRows) {
    for (int q = 0; q < fN; q++)
      for (int r = 0; r < 4; ++r)
        for (int j = 0; j < 5; ++j)
          H(row + r, j + q * cov_dim) = k[q].dEta[r][j];
    if (f3Constraint) { // mother theta, phi enter with -1/xi
      const double p = 1. / xi_iter(0, 0);
      const double th = m_iter(0 + fN * cov_dim, 0);
      const double ph = m_iter(1 + fN * cov_dim, 0);
      const double st = std::sin(th), ct = std::cos(th);
      const double sp = std::sin(ph), cp = std::cos(ph);
      H(row + 0, fN * cov_dim) = -p * ct * cp;
      H(row + 1, fN * cov_dim) = -p * ct * sp;
      H(row + 2, fN * cov_dim) = p * st;
      H(row + 0, 1 + fN * cov_dim) = p * st * sp;
      H(row + 1, 1 + fN * cov_dim) = -p * st * cp;
    }
    row += 4;
  }

  for (size_t c = 0; c < fMassVec.size(); ++c) {
    const std::vector<int> idx = massParticles(c);
    KinFitFourMomentum P;
    for (int q : idx)
      P.add(k[q]);
    for (int q : idx) {
      k[q].massDerivatives(P, 1., de, dv);
      for (int j = 0; j < 5; ++j)
        H(row, j + q * cov_dim) = de[j];
    }
    ++row;
  }

  // Vertex rows: central differences of the residuals of each track (the
  // closest approach in a general field has no closed form).
  if (fCommonVertex) {
    const int vo = vertexOffset();
    const TVector3 v(xi_iter(vo, 0), xi_iter(vo + 1, 0), xi_iter(vo + 2, 0));
    for (int q = 0; q < fN; ++q) {
      if (fIsNeutral[q])
        continue;
      double eta[5];
      for (int j = 0; j < 5; ++j)
        eta[j] = m_iter(j + q * cov_dim, 0);
      const double s0 = std::numeric_limits<double>::quiet_NaN();
      for (int j = 0; j < 5; ++j) {
        const double h = j == 0   ? kStepInvPRel * std::abs(eta[0])
                         : j <= 2 ? kStepAngle
                                  : kStepLength;
        double ep[5], em[5], rp[2], rm[2];
        std::copy(eta, eta + 5, ep);
        std::copy(eta, eta + 5, em);
        ep[j] += h;
        em[j] -= h;
        double sp = s0, sm = s0;
        vertexResiduals(ep, q, v, rp, &sp);
        vertexResiduals(em, q, v, rm, &sm);
        H(row, j + q * cov_dim) = (rp[0] - rm[0]) / (2 * h);
        H(row + 1, j + q * cov_dim) = (rp[1] - rm[1]) / (2 * h);
      }
      row += 2;
    }
  }

  return H;
}

// Jacobian (derivative of constraint equations with respect to unmeasured
// variables): columns [kinematic unknowns][vertex].
TMatrixD KinFitter::Fxi_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter) {
  if (fVerbose > 0) {
    std::cout << "--------------- KinFitter::Fxi_eval() -----------------"
              << std::endl;
  }

  const int nXi = nUnmeasured();
  if (nXi == 0)
    return TMatrixD();
  const std::vector<KinFitTrackKinematics> k = kinematics(m_iter, xi_iter);
  const bool fourRows = f3Constraint || f4Constraint || fMissingParticleConstraint;
  const int nVtxRows = fCommonVertex ? 2 * nCharged() : 0;
  const int nRows = nBaseRows(fMissingMassConstraint, fourRows) +
                    (int)fMassVec.size() + nVtxRows;
  const int vo = vertexOffset();
  TMatrixD H(std::max(nRows, 1), nXi);
  int row = 0;
  double de[5], dv[3];

  if (fMissingMassConstraint) {
    if (fCommonVertex) {
      KinFitFourMomentum P(fInit);
      for (int q = 0; q < fN; q++)
        P.add(k[q], -1.);
      for (int q = 0; q < fN; q++) {
        k[q].massDerivatives(P, -1., de, dv);
        for (int j = 0; j < 3; ++j)
          H(row, vo + j) += dv[j];
      }
    }
    ++row;
  }

  if (fourRows) {
    if (f3Constraint) { // d/d(1/p) of the mother term
      const double xi = xi_iter(0, 0);
      const double th = m_iter(0 + fN * cov_dim, 0);
      const double ph = m_iter(1 + fN * cov_dim, 0);
      H(row + 0, 0) = 1. / (xi * xi) * std::sin(th) * std::cos(ph);
      H(row + 1, 0) = 1. / (xi * xi) * std::sin(th) * std::sin(ph);
      H(row + 2, 0) = 1. / (xi * xi) * std::cos(th);
      H(row + 3, 0) = 1. / (xi * xi * xi) /
                      std::sqrt(1. / (xi * xi) + fM[fN] * fM[fN]);
    }
    if (fMissingParticleConstraint) {
      const double mx = xi_iter(0, 0), my = xi_iter(1, 0), mz = xi_iter(2, 0);
      const double E = std::sqrt(mx * mx + my * my + mz * mz + fM[fN] * fM[fN]);
      H(row + 0, 0) = H(row + 1, 1) = H(row + 2, 2) = 1.;
      H(row + 3, 0) = mx / E;
      H(row + 3, 1) = my / E;
      H(row + 3, 2) = mz / E;
    }
    if (fCommonVertex) // neutral hits: direction depends on the vertex
      for (int q = 0; q < fN; q++)
        for (int r = 0; r < 4; ++r)
          for (int j = 0; j < 3; ++j)
            H(row + r, vo + j) += k[q].dVtx[r][j];
    row += 4;
  }

  for (size_t c = 0; c < fMassVec.size(); ++c) {
    if (fCommonVertex) {
      const std::vector<int> idx = massParticles(c);
      KinFitFourMomentum P;
      for (int q : idx)
        P.add(k[q]);
      for (int q : idx) {
        k[q].massDerivatives(P, 1., de, dv);
        for (int j = 0; j < 3; ++j)
          H(row, vo + j) += dv[j];
      }
    }
    ++row;
  }

  if (fCommonVertex) {
    const TVector3 v(xi_iter(vo, 0), xi_iter(vo + 1, 0), xi_iter(vo + 2, 0));
    for (int q = 0; q < fN; ++q) {
      if (fIsNeutral[q])
        continue;
      double eta[5];
      for (int j = 0; j < 5; ++j)
        eta[j] = m_iter(j + q * cov_dim, 0);
      const double s0 = std::numeric_limits<double>::quiet_NaN();
      for (int j = 0; j < 3; ++j) {
        TVector3 vp(v), vm(v);
        vp[j] += kStepLength;
        vm[j] -= kStepLength;
        double rp[2], rm[2], sp = s0, sm = s0;
        vertexResiduals(eta, q, vp, rp, &sp);
        vertexResiduals(eta, q, vm, rm, &sm);
        H(row, vo + j) = (rp[0] - rm[0]) / (2 * kStepLength);
        H(row + 1, vo + j) = (rp[1] - rm[1]) / (2 * kStepLength);
      }
      row += 2;
    }
  }

  return H;
}

Bool_t KinFitter::fit() {

  if (fVerbose > 0) {
    std::cout << " ----------- KinFitter::fit() -----------" << std::endl;
    std::cout << "Common vertex constraint set: " << fCommonVertex
              << std::endl;
    std::cout << "3C set: " << f3Constraint << std::endl;
    std::cout << "4C set: " << f4Constraint << std::endl;
    std::cout << "Momentum constraint set: " << fMissingParticleConstraint
              << std::endl;
    std::cout << "Mass constraint set: " << fMassConstraint << std::endl;
    std::cout << "Missing Mass constraint set: " << fMissingMassConstraint
              << std::endl;
    std::cout << "" << std::endl;
  }

  fConverged = false;
  fIteration = 0;

  // --- Configuration checks (used to abort() or silently overwrite d) ---
  const int nBase = int(f3Constraint) + int(f4Constraint) +
                    int(fMissingParticleConstraint) +
                    int(fMissingMassConstraint);
  if (nBase == 0 && !fMassConstraint && !fCommonVertex) {
    std::cerr << "KinFitter::fit(): no constraint added" << std::endl;
    return kFALSE;
  }
  if (nBase > 1) {
    std::cerr << "KinFitter::fit(): 3C, 4C, MissingMass and MissingParticle "
                 "are exclusive (common-vertex and mass constraints can be "
                 "added to any of them)"
              << std::endl;
    return kFALSE;
  }
  if (f3Constraint && fCommonVertex) {
    std::cerr << "KinFitter::fit(): 3C cannot be combined with the common "
                 "vertex constraint"
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
  TMatrixD xi0;
  if (!startValues(xi0)) {
    std::cerr << "KinFitter::fit(): no start vertex (parallel tracks?)"
              << std::endl;
    return kFALSE;
  }
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
        // relative change; vertex coordinates on a scale of at least 1 mm
        const bool isVertex = fCommonVertex && i >= vertexOffset();
        const double scale = isVertex ? std::max(std::abs(xi0(i, 0)), 1.)
                                      : std::abs(xi0(i, 0));
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
    assign(fVMiss, Vx.GetSub(0, 2, 0, 2));
  }
  if (fCommonVertex) {
    const int vo = vertexOffset();
    fVertex.SetXYZ(x(vo, 0), x(vo + 1, 0), x(vo + 2, 0));
    assign(fVertexCov, Vx.GetSub(vo, vo + 2, vo, vo + 2));
  }
  if (fMissingFromMissingMass) {
    // p_miss = P_init - p of the single particle; its covariance from the
    // fitted parameter block of that particle.
    const KinFitTrackKinematics k = kinematics(alpha, xi)[0];
    fMissDaughter.SetVectM(fInit.Vect() - TVector3(k.px, k.py, k.pz), fMass);
    TMatrixD J(3, cov_dim);
    for (int i = 0; i < 3; ++i)
      for (size_t j = 0; j < cov_dim; ++j)
        J(i, j) = -k.dEta[i][j];
    const TMatrixD Vpart = fVFit.GetSub(0, cov_dim - 1, 0, cov_dim - 1);
    assign(fVMiss, J * Vpart * TMatrixD(TMatrixD::kTransposed, J));
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
    const int o = val * cov_dim;

    // Full 5x5 block of the fitted covariance.
    TMatrixD cov(cov_dim, cov_dim);
    for (size_t i = 0; i < cov_dim; ++i)
      for (size_t j = 0; j < cov_dim; ++j)
        cov(i, j) = fVFit(i + o, j + o);

    if (fIsNeutral[val]) {
      // (E, X, Y, Z): direction from the fitted vertex (or fixed origin)
      const TVector3 origin = fCommonVertex ? fVertex : fOrigin[val];
      const TVector3 hit(fYFit(o + 1, 0), fYFit(o + 2, 0), fYFit(o + 3, 0));
      const double E = fYFit(o, 0);
      const double p2 = E * E - fM[val] * fM[val];
      cand.SetVectM((p2 > 0 ? std::sqrt(p2) : 0.) * (hit - origin).Unit(),
                    fM[val]);
      cand.setMomentum(cand.P());
      cand.setThetaRad(cand.Theta());
      cand.setPhiRad(cand.Phi());
      cand.setHit(hit);
      cand.setVertexX(origin.X());
      cand.setVertexY(origin.Y());
      cand.setVertexZ(origin.Z());
      cand.setCovariance(cov);
      continue;
    }

    const double p = 1. / fYFit(0 + o, 0);
    double theta = fYFit(1 + o, 0);
    double phi = fYFit(2 + o, 0);
    double R = fYFit(3 + o, 0);
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
    cand.setZ(fYFit(4 + o, 0));
    cand.setCovariance(cov);
    if (transportToVertex() && fCharge[val] != 0) {
      // four-momentum at the fitted vertex (track parameters stay at the
      // reference point), so that daughters add up to the mother
      double eta[5];
      for (int j = 0; j < 5; ++j)
        eta[j] = fYFit(j + o, 0);
      cand.SetVectM(momentumAtVertex(eta, val, fVertex), fM[val]);
    }
    if (fCommonVertex) {
      cand.setVertexX(fVertex.X());
      cand.setVertexY(fVertex.Y());
      cand.setVertexZ(fVertex.Z());
    }
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
