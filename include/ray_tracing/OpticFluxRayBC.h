//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "OpticRayBC.h"

class OpticFluxRayBC : public OpticRayBC
{
public:
  OpticFluxRayBC(const InputParameters & params);

  static InputParameters validParams();

  void preExecuteStudy() override;

  void apply(const Elem * elem,
             const unsigned short intersected_side,
             const BoundaryID bnd_id,
             const Point & intersection_point,
             const std::shared_ptr<Ray> & ray,
             const bool applying_at_corner) override;

protected:
  /// Storage for accumulating the flux
  std::vector<Real> _values;
};
