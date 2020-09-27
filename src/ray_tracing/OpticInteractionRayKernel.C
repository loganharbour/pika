#include "OpticInteractionRayKernel.h"

// Local includes
#include "BasicLineSearch.h"
#include "OpticRayStudyBase.h"
#include "OpticsUtils.h"

// MOOSE includes
#include "Assembly.h"
#include "MathUtils.h"
#include "MooseUtils.h"
#include "ReflectRayBC.h"

registerMooseObject("PikaApp", OpticInteractionRayKernel);

InputParameters
OpticInteractionRayKernel::validParams()
{
  auto params = AuxRayKernel::validParams();

  params.addRequiredCoupledVar("phase", "The field variable that contains the phase");

  params.addRequiredRangeCheckedParam<Real>(
      "refractive_index_0",
      "refractive_index_0 > 0",
      "The refractive index for phase 0 (the phase that is < phase_change_threshold)");
  params.addRequiredRangeCheckedParam<Real>(
      "refractive_index_1",
      "refractive_index_1 > 0",
      "The refractive index for phase 1 (the phase that is > phase_change_threshold)");

  params.addParam<Real>("phase_change_threshold", 0.5, "The threshold to trigger a phase change");

  params.addRequiredParam<std::vector<Real>>(
      "attenuation_coefficient_0",
      "The linear attenuation coefficient for phase 0 (the phase "
      "that is < phase_change_threshold)");
  params.addRequiredParam<std::vector<Real>>(
      "attenuation_coefficient_1",
      "The linear attenuation coefficient for phase 1 (the phase "
      "that is < phase_change_threshold)");

  params.addParam<std::vector<Real>>("threshold", "The threshold at which to kill a Ray");

  return params;
}

OpticInteractionRayKernel::OpticInteractionRayKernel(const InputParameters & params)
  : AuxRayKernel(params),
    _optic_study(*dynamic_cast<OpticRayStudyBase *>(&_study)),
    _num_energy_groups(_optic_study.numEnergyGroups()),
    _phase(coupledValue("phase")),
    _otf_phase(_mesh.dimension(),
               getVar("phase", 0)->feType(),
               coupledDofValues("phase"),
               _fe_problem.assembly(_tid).elem(),
               {OnTheFlyVariableValue::VALUE}),
    _otf_grad_phase(_mesh.dimension(),
                    getVar("phase", 0)->feType(),
                    coupledDofValues("phase"),
                    _fe_problem.assembly(_tid).elem(),
                    {OnTheFlyVariableValue::GRAD_VALUE}),
    _refractive_index_0(getParam<Real>("refractive_index_0")),
    _refractive_index_1(getParam<Real>("refractive_index_1")),
    _phase_change_threshold(getParam<Real>("phase_change_threshold")),
    _attenuation_coefficient_0(getParam<std::vector<Real>>("attenuation_coefficient_0")),
    _attenuation_coefficient_1(getParam<std::vector<Real>>("attenuation_coefficient_1")),
    _threshold(params.isParamSetByUser("threshold")
                   ? getParam<std::vector<Real>>("threshold")
                   : std::vector<Real>(std::numeric_limits<Real>::max(), _num_energy_groups)),
    _just_interacted_aux_index(_study.registerRayAuxData("just_interacted")),
    _energy_index(_optic_study.energyDataIndex())
{
  std::string groups_error =
      "Must be the length of the number of groups, " + std::to_string(_num_energy_groups);
  if (_threshold.size() != _num_energy_groups)
    paramError("threshold", groups_error);
  if (_attenuation_coefficient_0.size() != _num_energy_groups)
    paramError("attenuation_coefficient_0", groups_error);
  if (_attenuation_coefficient_1.size() != _num_energy_groups)
    paramError("attenuation_coefficient_1", groups_error);

  if (MooseUtils::absoluteFuzzyEqual(_refractive_index_0, _refractive_index_1, 1.e-6))
    paramError("refractive_index_0", "should not be equal to refractive_index_1!");
}

void
OpticInteractionRayKernel::onSegment(const Elem * elem,
                                     const Point & start,
                                     const Point & end,
                                     const Real length,
                                     bool /* ends_in_elem */)
{
  if (underThreshold())
  {
    _ray->setShouldContinue(false);
    return;
  }

  // Phase at the start and end points on the segment
  const auto start_phase = _phase[0];
  const auto end_phase = _phase[1];

  // Whether or not we start at phase 0 (used for figuring out properties)
  const bool start_phase_0 = start_phase < _phase_change_threshold;

  // The attenuation coefficient
  const auto & attenuation_coefficient =
      start_phase_0 ? _attenuation_coefficient_0 : _attenuation_coefficient_1;

  // If we just interacted, it means the start point is in the middle of an
  // element so there's nothing to do here except remove the flag
  auto & just_interacted = _ray->auxData(_just_interacted_aux_index);

  // Phase changes happens in this element
  if (!just_interacted &&
      ((start_phase >= _phase_change_threshold && end_phase <= _phase_change_threshold) ||
       (end_phase >= _phase_change_threshold && start_phase <= _phase_change_threshold)))
  {
    // If we interact at the end point, let the next segment handle it. Without
    // doing this, the interaction will happen both in this elem and the next elem
    if (MooseUtils::absoluteFuzzyEqual(end_phase, _phase_change_threshold))
      return;

    // Current direction of the ray
    const auto & direction = _ray->direction();

    // Normal to the phase change surface at the midpoint of the segment
    const auto phase_normal = _otf_grad_phase.gradValue(0.5 * (start + end)).unit();

    // Indices of refraction at the start and end points on the segment
    const auto n1 = start_phase_0 ? _refractive_index_0 : _refractive_index_1;
    const auto n2 = start_phase_0 ? _refractive_index_1 : _refractive_index_0;

    // Point in this element at which we interact (the phase changes)
    const auto interaction_point = phaseChangePoint(start, end, direction, start_phase, end_phase);

    // The refraction direction and type
    Point refracted_direction;
    const auto snell_type =
        OpticsUtils::snell(direction, phase_normal, n1, n2, refracted_direction);
    // The reflected direction
    const auto reflected_direction = ReflectRayBC::reflectedDirection(direction, phase_normal);

    attenuateRay(start, interaction_point, attenuation_coefficient);

    // Without internal reflection, refract the ray and potentially create a reflected ray
    if (snell_type != OpticsUtils::SNELL_INTERNAL_REFLECTION)
    {
      // Refract the current ray
      changeRay(interaction_point, refracted_direction);

      // Reflect the ray, only if we're above the threshold
      if (snell_type != OpticsUtils::SNELL_PARALLEL)
      {
        const auto R = OpticsUtils::reflectionCoefficient(
            direction, reflected_direction, phase_normal, n1, n2);

        // Created the reflected ray if over the threshold
        if (!underThreshold(R))
        {
          auto & reflected_ray = _optic_study.addReflectedRay(
              interaction_point, elem, reflected_direction, *_ray, R, _tid);
          reflected_ray.setAuxData(_just_interacted_aux_index, 1);
        }

        // Energy loss in the refracted ray
        for (unsigned int g = 0; g < _num_energy_groups; ++g)
          _ray->data(_energy_index + g) *= (1 - R);
      }
    }
    // Angle above the critical angle: we only reflect
    else
      changeRay(interaction_point, reflected_direction);

    // Mark that this Ray just interacted
    just_interacted = 1;
  }
  // Phase change didn't happen in this element
  else
  {
    just_interacted = 0;
    attenuateRay(start, end, attenuation_coefficient);
  }

  if (underThreshold())
    _ray->setShouldContinue(false);
}

Point
OpticInteractionRayKernel::phaseChangePoint(const Point & start,
                                            const Point & end,
                                            const Point & direction,
                                            const Real start_phase,
                                            const Real end_phase)
{
  // No need to line search if it happens at the beginning
  if (MooseUtils::absoluteFuzzyEqual(start_phase, _phase_change_threshold))
    return start;

  // Segment distance used for parameterization
  const Real distance = (end - start).norm();

  // Parameterized function for the line search
  const auto f = [this, &start, &distance, &direction](const Real t) {
    return _otf_phase.value(start + t * distance * direction) - _phase_change_threshold;
  };

  // Location on segment (start is t = 0, end is t = 1) where phase = _threshold
  const auto t = BasicLineSearch::brents(
      0, 1, start_phase - _phase_change_threshold, end_phase - _phase_change_threshold, f);

  // Change back into physical space
  return start + t * distance * direction;
}

void
OpticInteractionRayKernel::attenuateRay(const Point & start,
                                        const Point & end,
                                        const std::vector<Real> & attenuation_coefficient)
{
  const auto distance = (end - start).norm();
  const auto attenuation = std::exp(-distance * attenuation_coefficient[0]);
  const Real deposition = _ray->data(_energy_index) * (1 - attenuation);
  _ray->data(_energy_index) *= attenuation;
  addValue(deposition);
}

bool
OpticInteractionRayKernel::underThreshold(const Real factor) const
{
  for (unsigned int g = 0; g < _num_energy_groups; ++g)
    if (factor * _ray->data(_energy_index + g) > _threshold[g])
      return false;

  return true;
}
