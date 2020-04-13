#include "OpticDepositionRayKernel.h"

#include "OpticRayStudyBase.h"

registerMooseObject("PikaApp", OpticDepositionRayKernel);

InputParameters
OpticDepositionRayKernel::validParams()
{
  auto params = AuxRayKernel<Real>::validParams();

  return params;
}

OpticDepositionRayKernel::OpticDepositionRayKernel(const InputParameters & params)
  : AuxRayKernel<Real>(params),
    _optic_study(*dynamic_cast<OpticRayStudyBase *>(&_study)),
    _num_energy_groups(_optic_study.numEnergyGroups()),
    _optic_data(declareRayVariable(_optic_study.opticRayVariableName(), _num_energy_groups))
{
  if (_num_energy_groups == 0)
    mooseError("RayKernel ", name(), " is not needed with zero energy groups");
}

void
OpticDepositionRayKernel::onSegment(const Elem * /* elem */,
                                    const Point & start,
                                    const Point & end,
                                    bool /* ends_in_elem */)
{
  // 'start' and 'end' are the true traced start and end points before
  // refraction (if any). If a refraction kernel just changed a Ray, we now have
  // two segments:
  // [start -> _ray->start()] and [_ray->start() -> new end]
  //     ^ tracing now ^             ^ tracing next ^
  const auto distance =
      _ray->trajectoryChanged() ? (start - _ray->start()).norm() : (end - start).norm();

  Real value = 0;
  for (unsigned int g = 0; g < _num_energy_groups; ++g)
    value += distance * _optic_data[g];

  std::cerr << "[" << name() << "_" << _ray->id() << "]: distance " << distance << " adding "
            << value << std::endl;

  addValue(value);
}
