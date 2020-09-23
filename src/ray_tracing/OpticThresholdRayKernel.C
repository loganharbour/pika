#include "OpticThresholdRayKernel.h"

// Local includes
#include "OpticRayStudyBase.h"

registerMooseObject("PikaApp", OpticThresholdRayKernel);

InputParameters
OpticThresholdRayKernel::validParams()
{
  auto params = RayKernelBase::validParams();

  params.addRequiredParam<std::vector<Real>>("threshold", "The threshold at which to kill a Ray");

  return params;
}

OpticThresholdRayKernel::OpticThresholdRayKernel(const InputParameters & params)
  : RayKernelBase(params),
    _optic_study(*dynamic_cast<OpticRayStudyBase *>(&_study)),
    _num_energy_groups(_optic_study.numEnergyGroups()),
    _energy_index(_optic_study.energyDataIndex()),
    _threshold(getParam<std::vector<Real>>("threshold"))
{
  if (_threshold.size() != _num_energy_groups)
    paramError("threshold", "Must be the length of the number of groups: ", _num_energy_groups);
}

void
OpticThresholdRayKernel::onSegment(const Elem * /* elem */,
                                   const Point & /* start */,
                                   const Point & /* end */,
                                   bool /* ends_in_elem */)
{
  for (unsigned int g = 0; g < _num_energy_groups; ++g)
    if (_ray->data(_energy_index + g) > _threshold[g])
      return;

  _ray->setShouldContinue(false);
}
