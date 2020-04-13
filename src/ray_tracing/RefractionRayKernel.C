#include "RefractionRayKernel.h"

// Local includes
#include "BasicLineSearch.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseUtils.h"

registerMooseObject("PikaApp", RefractionRayKernel);

InputParameters
RefractionRayKernel::validParams()
{
  auto params = RayKernelBase::validParams();

  params.addRequiredCoupledVar("phase", "The field variable that contains the phase");
  params.addCoupledVar("refractive_index", "The field variable that contains the refractive index");

  return params;
}

RefractionRayKernel::RefractionRayKernel(const InputParameters & params)
  : RayKernelBase(params),
    _phase(coupledValue("phase")),
    _grad_phase(coupledGradient("phase")),
    _otf_phase(_mesh.dimension(),
               getVar("phase", 0)->feType(),
               coupledDofValues("phase"),
               _fe_problem.assembly(_tid).elem()),
    _refractive_index(coupledValue("refractive_index")),
    _q_point(_fe_problem.assembly(_tid).qPoints()),
    _just_refracted(false),
    _threshold(0.5)
{
}

void
RefractionRayKernel::onSegment(const Elem * /* elem */,
                               const Point & start,
                               const Point & end,
                               bool /* ends_in_elem */)
{
  mooseAssert(start.absolute_fuzzy_equals(_q_point[0]), "Quadrature point 0 must be segment start");
  mooseAssert(end.absolute_fuzzy_equals(_q_point[1]), "Quadrature point 1 must be segment end");

  std::cerr << std::endl
            << "[" << name() << "_" << _ray->id() << "]: " << start << " -> " << end << std::endl;
  std::cerr << "  phase start = " << _phase[0] << std::endl;
  std::cerr << "  phase end = " << _phase[1] << std::endl;

  // If we just refracted, it means the start point is in the middle of an
  // element so there's nothing to do here except remove the flag
  if (_just_refracted)
  {
    _just_refracted = false;
    return;
  }

  // Phase at the start and end points on the segment
  const auto start_phase = _phase[0];
  const auto end_phase = _phase[1];

  // If we refract at the end point, let the next segment handle it. Without
  // doing this, the refraction will happen both in this elem and the next elem
  if (MooseUtils::absoluteFuzzyEqual(end_phase, _threshold))
  {
    std::cerr << "  refracting on next segment" << std::endl;
    return;
  }

  // Phase changes happens in this element
  if ((start_phase >= _threshold && end_phase <= _threshold) ||
      (end_phase >= _threshold && start_phase <= _threshold))
  {
    // Normal to the phase change surface
    const auto phase_normal = _grad_phase[1].unit();
    // Indices of refraction at the start and end points on the segment
    const auto start_r = _refractive_index[0];
    const auto end_r = _refractive_index[1];

    // Point in this element at which we refract
    const auto refracted_point = phaseChange(start, end, _ray->direction(), start_phase, end_phase);

    // Direction we refract to
    const auto refracted_direction = snell(_ray->direction(), phase_normal, start_r, end_r);

    // Refract!
    changeRay(refracted_point, refracted_direction);

    std::cerr << "  refracted at point " << refracted_point << std::endl;
    std::cerr << "  refracted to direction " << refracted_direction << std::endl;

    _just_refracted = true;
  }
}

void
RefractionRayKernel::preTrace(const std::shared_ptr<Ray> & ray)
{
  RayKernelBase::preTrace(ray);

  // On a new Ray, reset the refraction flag
  _just_refracted = false;
}

Point
RefractionRayKernel::snell(const Point & direction,
                           const Point & normal,
                           const Real r1,
                           const Real r2) const
{
  mooseAssert(r1 > 0, "Invalid index of refraction");
  mooseAssert(r2 > 0, "Invalid index of refraction");

  // Dot product between normal and direction
  const Real c = std::abs(normal * direction);

  // Direction and normal are parallel
  if (c > 1.0 - TOLERANCE)
    return direction;

  const Real r = r1 / r2;
  return (r * direction + (r * c - std::sqrt(1 - r * r * (1 - c * c))) * normal).unit();
}

Point
RefractionRayKernel::phaseChange(const Point & start,
                                 const Point & end,
                                 const Point & direction,
                                 const Real & start_phase,
                                 const Real & end_phase)
{
  // No need to line search if it happens at the beginning
  if (MooseUtils::absoluteFuzzyEqual(start_phase, _threshold))
    return start;

  // Segment distance used for parameterization
  const Real distance = (end - start).norm();

  // Parameterized function for the line search
  const auto f = [this, &start, &distance, &direction](const Real t) {
    return _otf_phase(start + t * distance * direction) - _threshold;
  };

  // Location on segment (start is t = 0, end is t = 1) where phase = _threshold
  const auto t = BasicLineSearch::brents(0, 1, start_phase - _threshold, end_phase - _threshold, f);

  // Change back into physical space
  return start + t * distance * direction;
}
