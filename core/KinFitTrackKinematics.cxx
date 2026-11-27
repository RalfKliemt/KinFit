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

#include "KinFitTrackKinematics.h"
#include "KinFitFourMomentum.h"

#include <cmath>

KinFitTrackKinematics::KinFitTrackKinematics(const TMatrixD &m, int offset,
                                             double mass) {
  const double th = m(offset + 1, 0), ph = m(offset + 2, 0);
  const double p = 1. / m(offset, 0);
  const double st = std::sin(th), ct = std::cos(th);
  const double sp = std::sin(ph), cp = std::cos(ph);
  E = std::sqrt(p * p + mass * mass);
  px = p * st * cp;
  py = p * st * sp;
  pz = p * ct;
  const double p2 = p * p;
  // d/d(1/p): dp/d(1/p) = -p^2 u, dE/d(1/p) = -p^3/E
  dEta[0][0] = -p2 * st * cp;
  dEta[1][0] = -p2 * st * sp;
  dEta[2][0] = -p2 * ct;
  dEta[3][0] = -p2 * p / E;
  // d/d(theta), d/d(phi); E depends on 1/p only; R, Z do not enter
  dEta[0][1] = p * ct * cp;
  dEta[1][1] = p * ct * sp;
  dEta[2][1] = -p * st;
  dEta[0][2] = -p * st * sp;
  dEta[1][2] = p * st * cp;
}

KinFitTrackKinematics KinFitTrackKinematics::neutralHit(const TMatrixD &m,
                                                        int offset,
                                                        double mass,
                                                        const TVector3 &vertex) {
  KinFitTrackKinematics k;
  const double e = m(offset, 0);
  const TVector3 d = TVector3(m(offset + 1, 0), m(offset + 2, 0),
                              m(offset + 3, 0)) -
                     vertex;
  const double L = d.Mag();
  const TVector3 u = d * (1. / L);
  const double p2 = e * e - mass * mass;
  const double p = p2 > 0 ? std::sqrt(p2) : 0.;
  k.E = e;
  k.px = p * u.X();
  k.py = p * u.Y();
  k.pz = p * u.Z();
  // d/dE: dp/dE = E/p (photon: 1), dE/dE = 1
  const double dpdE = p > 0 ? e / p : 0.;
  const double uu[3] = {u.X(), u.Y(), u.Z()};
  for (int i = 0; i < 3; ++i)
    k.dEta[i][0] = dpdE * uu[i];
  k.dEta[3][0] = 1.;
  // d(p u)/d(hit) = p (I - u u^T) / L = -d(p u)/d(vertex)
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const double v = p * ((i == j ? 1. : 0.) - uu[i] * uu[j]) / L;
      k.dEta[i][1 + j] = v;
      k.dVtx[i][j] = -v;
    }
  return k;
}

void KinFitTrackKinematics::massDerivatives(const KinFitFourMomentum &P,
                                            double sign, double dEtaOut[5],
                                            double dVtxOut[3]) const {
  for (int k = 0; k < 5; ++k)
    dEtaOut[k] = sign * 2. *
                 (P.E * dEta[3][k] - P.px * dEta[0][k] - P.py * dEta[1][k] -
                  P.pz * dEta[2][k]);
  for (int j = 0; j < 3; ++j)
    dVtxOut[j] = sign * 2. *
                 (P.E * dVtx[3][j] - P.px * dVtx[0][j] - P.py * dVtx[1][j] -
                  P.pz * dVtx[2][j]);
}
