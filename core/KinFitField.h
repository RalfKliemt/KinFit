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
 * KinFitField.h
 *
 * Abstract magnetic field B(x, y, z) used for track propagation.
 * Units: positions in mm, field in Tesla.
 */

#pragma once

#include <TVector3.h>

class KinFitField {
public:
  virtual ~KinFitField() = default;

  /** @brief Field in Tesla at (x, y, z) in mm. */
  virtual TVector3 getField(double x, double y, double z) const = 0;

  /** @brief Field in Tesla at a position in mm. */
  TVector3 getField(const TVector3 &pos) const {
    return getField(pos.X(), pos.Y(), pos.Z());
  }

  /** @brief True if the field is the same everywhere. Tracks in a uniform
   * field are exact helices (no numerical propagation needed). */
  virtual bool isUniform() const { return false; }

  /** @brief Polymorphic copy. */
  virtual KinFitField *clone() const = 0;
};
