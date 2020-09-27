//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OpticFluxRayBC.h"

registerMooseObject("PikaApp", OpticFluxRayBC);

InputParameters
OpticFluxRayBC::validParams()
{
  return OpticRayBC::validParams();
}

OpticFluxRayBC::OpticFluxRayBC(const InputParameters & params)
  : OpticRayBC(params), _values(_num_energy_groups)
{
}

void
OpticFluxRayBC::apply(const Elem * /* elem */,
                      const unsigned short /* intersected_side */,
                      const BoundaryID /* bnd_id */,
                      const Point & /* intersection_point */,
                      const std::shared_ptr<Ray> & ray,
                      const unsigned int num_applying)
{
  for (unsigned int g = 0; g < _num_energy_groups; ++g)
    _values[g] += ray->data(_energy_index + g) / (Real)num_applying;
}

void
OpticFluxRayBC::preExecuteStudy()
{
  // Zero the values before execution
  for (auto & value : _values)
    value = 0;
}
