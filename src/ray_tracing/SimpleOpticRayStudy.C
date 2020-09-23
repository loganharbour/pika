#include "SimpleOpticRayStudy.h"

registerMooseObject("PikaApp", SimpleOpticRayStudy);

InputParameters
SimpleOpticRayStudy::validParams()
{
  auto params = OpticRayStudyBase::validParams();

  params.addRequiredParam<std::vector<Point>>("start_points", "The points to start Rays from");
  params.addParam<std::vector<Point>>("directions",
                                      "The directions to spawn Rays in (they do not need to be "
                                      "normalized)");
  params.addParam<std::vector<Real>>("data", "The data to set on the Rays");

  return params;
}

SimpleOpticRayStudy::SimpleOpticRayStudy(const InputParameters & parameters)
  : OpticRayStudyBase(parameters),
    _start_points(getParam<std::vector<Point>>("start_points")),
    _directions(getParam<std::vector<Point>>("directions")),
    _data(parameters.isParamSetByUser("data") ? &getParam<std::vector<Real>>("data") : nullptr)
{
  if (_start_points.size() != _directions.size())
    paramError("start_points", "Not the same size as directions");
  if (_data && _num_energy_groups != _data->size())
    paramError("data", "size does not equal number of energy groups");
  if (!_data && _num_energy_groups > 0)
    paramError("data", "not supplied but num_energy_groups > 0");
  if (!isRectangularDomain())
    mooseError("Optic ray study ", name(), " does not support non-rectangular domains");
}

void
SimpleOpticRayStudy::defineRays()
{
  for (unsigned int i = 0; i < _start_points.size(); ++i)
  {
    Ray & ray = defineRay(_start_points[i], _directions[i], (RayID)i);

    // Set the optic data if there any
    if (_data)
      for (unsigned int g = 0; g < _num_energy_groups; ++g)
        ray.data(_energy_data_index + g) = (*_data)[g];
  }
}
