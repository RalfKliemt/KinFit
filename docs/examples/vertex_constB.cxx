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
 * vertex_constB.cxx
 *
 * Self-contained example for vertex finding in a constant magnetic field.
 * Charged tracks are generated from a known vertex in Bz = 1 T, propagated
 * along their helices to the point of closest approach (PCA) to the z axis,
 * and stored as KinFitParticles (R, Z, theta, phi at the PCA). Then the
 * straight-line KFitVertexFinder and the helix KFitVertexFinderConstB are
 * compared.
 *
 * Units: MeV/c, mm, Tesla.
 */

// In an interpreted ROOT session the KinFit classes are autoloaded from the
// library dictionary; including the headers again would clash with it.
#ifndef __CLING__
#include "KFitVertexFinder.h"
#include "KFitVertexFinderConstB.h"
#include "KinFitParticle.h"
#endif

#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <TVector3.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

/** Propagate a helix from vertex `v` with momentum `p` and charge `q` in a
 *  field Bz [T] to its PCA to the z axis and return it as KinFitParticle. */
KinFitParticle trackAtPCA(const TVector3 &v, const TVector3 &p, int q,
                          double bz, double mass) {
  const double pt = p.Perp();
  const double rho = pt / (0.299792458 * std::abs(q * bz)); // mm
  // Turning sense in the xy plane: +1 counter-clockwise. A positive charge
  // in +Bz turns clockwise (F = q v x B).
  const double h = (q * bz > 0) ? -1.0 : 1.0;
  const double phi0 = p.Phi();

  // Circle centre and PCA of the circle to the origin.
  const TVector2 c(v.X() - h * rho * std::sin(phi0),
                   v.Y() + h * rho * std::cos(phi0));
  const TVector2 pca = c - rho * c.Unit();

  // Turning angle from the vertex to the PCA (shortest way round).
  const double aV = std::atan2(v.Y() - c.Y(), v.X() - c.X());
  const double aP = std::atan2(pca.Y() - c.Y(), pca.X() - c.X());
  const double dA = TVector2::Phi_mpi_pi(aP - aV);
  const double s = h * dA * rho; // signed transverse arc length

  const double phiPCA = phi0 + dA;
  const double zPCA = v.Z() + s * p.Z() / pt;
  const double R = -pca.X() * std::sin(phiPCA) + pca.Y() * std::cos(phiPCA);

  TLorentzVector lv;
  lv.SetXYZM(pt * std::cos(phiPCA), pt * std::sin(phiPCA), p.Z(), mass);
  KinFitParticle cand(lv, R, zPCA);
  cand.setCharge(q);
  return cand;
}

} // namespace

int vertex_constB(int nEvents = 10) {
  const double bz = 1.0; // Tesla
  TRandom3 rnd(42);

  std::printf("Vertex finding in Bz = %.1f T, 3 tracks/event, no smearing\n",
              bz);
  std::printf("%-28s %16s %16s\n", "true vertex [mm]", "|d| helix [mm]",
              "|d| lines [mm]");

  for (int ev = 0; ev < nEvents; ++ev) {
    const TVector3 vTrue(rnd.Uniform(-30, 30), rnd.Uniform(-30, 30),
                         rnd.Uniform(-50, 50));
    std::vector<KinFitParticle> cands;
    for (int t = 0; t < 3; ++t) {
      TVector3 p;
      p.SetPtThetaPhi(rnd.Uniform(150, 800), rnd.Uniform(0.3, 1.3),
                      rnd.Uniform(-TMath::Pi(), TMath::Pi()));
      cands.push_back(trackAtPCA(vTrue, p, (t % 2) ? -1 : +1, bz, 139.570));
    }

    KFitVertexFinderConstB helixFinder(cands, bz);
    KFitVertexFinder lineFinder(cands);

    std::printf("(%7.2f, %7.2f, %7.2f)    %16.6f %16.3f\n", vTrue.X(),
                vTrue.Y(), vTrue.Z(),
                (helixFinder.getVertex() - vTrue).Mag(),
                (lineFinder.getVertex() - vTrue).Mag());
  }
  return 0;
}

#ifndef __CLING__
int main(int argc, char **argv) {
  return vertex_constB(argc > 1 ? std::atoi(argv[1]) : 10);
}
#endif
