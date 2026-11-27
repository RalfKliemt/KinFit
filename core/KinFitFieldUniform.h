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
 * KinFitFieldUniform.h
 *
 * Constant magnetic field (any direction), including B = 0. Covers a
 * solenoid-like constant Bz and a dipole-like constant By.
 */

#pragma once

#include "KinFitField.h"

class KinFitFieldUniform : public KinFitField {
private:
  TVector3 fB; // Field in Tesla

public:
  /** @brief B = 0 (straight tracks). */
  KinFitFieldUniform() = default;
  /** @brief Constant field (bx, by, bz) in Tesla. */
  KinFitFieldUniform(double bx, double by, double bz) : fB(bx, by, bz) {}

  /** @brief Constant Bz in Tesla. */
  static KinFitFieldUniform Bz(double bz) { return {0., 0., bz}; }
  /** @brief Constant By in Tesla. */
  static KinFitFieldUniform By(double by) { return {0., by, 0.}; }

  TVector3 getField(double, double, double) const override { return fB; }
  using KinFitField::getField;
  bool isUniform() const override { return true; }
  KinFitField *clone() const override { return new KinFitFieldUniform(*this); }
};
