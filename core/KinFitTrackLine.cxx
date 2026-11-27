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

#include "KinFitTrackLine.h"
#include "KinFitTrackKinematics.h"

KinFitTrackLine::KinFitTrackLine(const TMatrixD &m, int offset,
                                 const KinFitTrackKinematics &k) {
  const double R = m(offset + 3, 0), Z = m(offset + 4, 0);
  b.SetXYZ(-R * k.sp, R * k.cp, Z);
  d.SetXYZ(k.st * k.cp, k.st * k.sp, k.ct);
  dd_dTheta.SetXYZ(k.ct * k.cp, k.ct * k.sp, -k.st);
  dd_dPhi.SetXYZ(-k.st * k.sp, k.st * k.cp, 0.);
  db_dPhi.SetXYZ(-R * k.cp, -R * k.sp, 0.);
  db_dR.SetXYZ(-k.sp, k.cp, 0.);
}

double KinFitTrackLine::vertexConstraint(const KinFitTrackLine &l1,
                                         const KinFitTrackLine &l2) {
  return l1.d.Cross(l2.d).Dot(l1.b - l2.b);
}

void KinFitTrackLine::vertexJacobian(TMatrixD &H, int row,
                                     const KinFitTrackLine &l1, int col1,
                                     const KinFitTrackLine &l2, int col2) {
  const TVector3 n = l1.d.Cross(l2.d);
  const TVector3 db = l1.b - l2.b;
  // 1/p does not enter the vertex constraint.
  H(row, col1 + 1) = l1.dd_dTheta.Cross(l2.d).Dot(db);
  H(row, col1 + 2) = l1.dd_dPhi.Cross(l2.d).Dot(db) + n.Dot(l1.db_dPhi);
  H(row, col1 + 3) = n.Dot(l1.db_dR);
  H(row, col1 + 4) = n.Z();
  H(row, col2 + 1) = l1.d.Cross(l2.dd_dTheta).Dot(db);
  H(row, col2 + 2) = l1.d.Cross(l2.dd_dPhi).Dot(db) - n.Dot(l2.db_dPhi);
  H(row, col2 + 3) = -n.Dot(l2.db_dR);
  H(row, col2 + 4) = -n.Z();
}
