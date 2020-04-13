#include "OpticRayStudyBase.h"

// MOOSE incldues
#include "RayStudyTraceRay.h"

InputParameters
OpticRayStudyBase::validParams()
{
  auto params = RepeatableRayStudyBase::validParams();

  params.addParam<unsigned int>("energy_groups", 0, "Number of energy groups");

  return params;
}

OpticRayStudyBase::OpticRayStudyBase(const InputParameters & parameters)
  : RepeatableRayStudyBase(parameters),
    _num_energy_groups(getParam<unsigned int>("energy_groups")),
    _optic_ray_variable_name("optic_data")
{
  // Disable ray registration as we don't want to name/register every Ray
  _use_ray_registration = false;
  // We will resize the data manually when we generate the Rays
  _set_data_at_generation = false;

  // Subdomain setup in the trace does not depend on the Ray here
  for (auto & trace_ray : _threaded_trace_ray)
    trace_ray->setRayDependentSubdomainSetup(false);

  if (_num_energy_groups)
    declareRayVariable(_optic_ray_variable_name, 0, _num_energy_groups);
}

void
OpticRayStudyBase::buildSegmentQuadrature(const Point & start,
                                          const Point & end,
                                          std::vector<Point> & points,
                                          std::vector<Real> & weights) const
{
  // For an optic study, we only care about reinit at both ends of the segment
  // instead of the standard 1D quadrature along the segment that comes from
  // RayTracingStudy
  points = {start, end};
  weights = {1, 1};
}
