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
 * KinFitFieldFunction.h
 *
 * Field from a user function B(x, y, z) (mm -> Tesla), e.g. a wrapper
 * around an experiment's field map.
 */

#pragma once

#include "KinFitField.h"

#include <functional>

class KinFitFieldFunction : public KinFitField {
public:
  using Function = std::function<TVector3(double x, double y, double z)>;

private:
  Function fFunc;
  bool fUniform = false;

public:
  /** @param func returns B [T] at (x, y, z) [mm]
   *  @param uniform set true only if func is constant (enables helices) */
  explicit KinFitFieldFunction(Function func, bool uniform = false)
      : fFunc(std::move(func)), fUniform(uniform) {}

  TVector3 getField(double x, double y, double z) const override {
    return fFunc(x, y, z);
  }
  using KinFitField::getField;
  bool isUniform() const override { return fUniform; }
  KinFitField *clone() const override { return new KinFitFieldFunction(*this); }
};
