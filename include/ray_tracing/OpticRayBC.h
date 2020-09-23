//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RayBC.h"

// Forward declarations
class OpticRayStudyBase;

class OpticRayBC : public RayBC
{
public:
  OpticRayBC(const InputParameters & params);

  static InputParameters validParams();

  virtual void apply(const Elem * elem,
                     const unsigned short intersected_side,
                     const BoundaryID bnd_id,
                     const Point & intersection_point,
                     const std::shared_ptr<Ray> & ray,
                     const bool applying_at_corner) = 0;

protected:
  /// The OpticRayStudyBase
  OpticRayStudyBase & _optic_study;
  /// The number of energy groups in the optic study
  const unsigned int _num_energy_groups;
  /// The index into the Ray's data for the optic data
  const RayDataIndex _energy_index;
};
