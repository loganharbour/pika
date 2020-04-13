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
  // Start with an ID for Rays that is unique to this processor
  _next_id = DofObject::invalid_id / (dof_id_type)_comm.size() * (dof_id_type)_pid;

  // Starting index for the data if we need to fill it
  const auto optic_data_index = _data ? rayVariableIndex("optic_data", 0) : 0;

  for (unsigned int i = 0; i < _start_points.size(); ++i)
  {
    std::shared_ptr<Ray> ray = _ray_pool.acquire();
    _rays.push_back(ray);

    const auto start = _start_points[i];
    const auto direction = _directions[i].unit();
    const auto end = rayDomainIntersection(_start_points[i], direction);

    ray->setStart(start);
    ray->setEnd(end);
    ray->setDirection(direction);
    ray->setID(_next_id++);

    // Resize data regardless in case other things declared ray variables
    ray->data().resize(rayVariableSize(), 0);
    // Set the optic data if there any
    if (_data)
      for (unsigned int i = optic_data_index; i < optic_data_index + _data->size(); ++i)
        ray->data(i) = (*_data)[i];
  }
}
