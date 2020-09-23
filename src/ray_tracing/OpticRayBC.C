//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OpticRayBC.h"

// Local includes
#include "OpticRayStudyBase.h"

InputParameters
OpticRayBC::validParams()
{
  return RayBC::validParams();
}

OpticRayBC::OpticRayBC(const InputParameters & params)
  : RayBC(params),
    _optic_study(*dynamic_cast<OpticRayStudyBase *>(&_study)),
    _num_energy_groups(_optic_study.numEnergyGroups()),
    _energy_index(_optic_study.energyDataIndex())
{
}
