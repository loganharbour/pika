#pragma once

#include "AuxRayKernel.h"

#include "OnTheFlyVariableValue.h"

class OpticRayStudyBase;

/**
 * RayKernel that physically refracts a Ray at the point a phase hits a
 * threadhold within an element.
 */
class OpticInteractionRayKernel : public AuxRayKernel
{
public:
  OpticInteractionRayKernel(const InputParameters & params);

  static InputParameters validParams();

  virtual void onSegment(const Elem * elem,
                         const Point & start,
                         const Point & end,
                         const Real length,
                         bool ends_in_elem) override;

protected:
  /**
   * The point at which the phase changes on a line segment within element,
   * which is the point on the segment where the phase value is 0.5
   * @param start The start point of the segment
   * @param end The end point of the segment
   * @param direction The direction of the segment
   * @param start_phase The phase value at the start point of the segment
   * @param end_phase The phase value at the end point of the segment
   * @return The point on the segment where the phase is 0.5
   */
  Point phaseChangePoint(const Point & start,
                         const Point & end,
                         const Point & direction,
                         const Real start_phase,
                         const Real end_phase);

  void attenuateRay(const Point & start,
                    const Point & end,
                    const std::vector<Real> & attenuation_coefficient);

  bool underThreshold(const Real factor = 1) const;

  /// The OpticRayStudyBase
  OpticRayStudyBase & _optic_study;
  /// The number of energy groups in the study
  const unsigned int _num_energy_groups;

  /// The field variable that contains the phase
  const VariableValue & _phase;
  /// On the fly value for the phase for quick reinits
  OnTheFlyVariableValue _otf_phase;
  /// On the fly grad value for the phase for quick reinits
  OnTheFlyVariableValue _otf_grad_phase;

  /// The refractive index for phase 0
  const Real _refractive_index_0;
  /// The refractive index for phase 1
  const Real _refractive_index_1;

  /// The phase threshold for refraction
  const Real _phase_change_threshold;

  /// The group-wise linear attenuation coefficients for phase 0
  const std::vector<Real> _attenuation_coefficient_0;
  /// The group-wise linear attenuation coefficients for phase 1
  const std::vector<Real> _attenuation_coefficient_1;

  const std::vector<Real> _threshold;

  /// The index into the Ray's aux data for the just interacted flag
  const RayDataIndex _just_interacted_aux_index;
  /// The index into the Ray's data for the start of the group-wis eenergy data
  const RayDataIndex _energy_index;
};
