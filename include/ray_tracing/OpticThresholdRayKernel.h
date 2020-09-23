#pragma once

#include "RayKernelBase.h"

// Forward declarations
class OpticRayStudyBase;

/**
 * RayKernel that physically refracts a Ray at the point a phase hits a
 * threadhold within an element.
 */
class OpticThresholdRayKernel : public RayKernelBase
{
public:
  OpticThresholdRayKernel(const InputParameters & params);

  static InputParameters validParams();

  virtual void
  onSegment(const Elem * elem, const Point & start, const Point & end, bool ends_in_elem) override;

protected:
  /// The OpticRayStudyBase
  OpticRayStudyBase & _optic_study;
  /// The number of energy groups in the optic study
  const unsigned int _num_energy_groups;
  /// The index into the Ray's data for the optic data
  const RayDataIndex _energy_index;

  /// The threshold in each group that must be met to kill a Ray
  const std::vector<Real> _threshold;
};
