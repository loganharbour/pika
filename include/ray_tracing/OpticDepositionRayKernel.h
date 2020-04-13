#pragma once

#include "AuxRayKernel.h"

class OpticRayStudyBase;

/**
 * RayKernel that deposits energy from Rays into an aux variable
 */
class OpticDepositionRayKernel : public AuxRayKernel<Real>
{
public:
  OpticDepositionRayKernel(const InputParameters & params);

  static InputParameters validParams();

  virtual void
  onSegment(const Elem * elem, const Point & start, const Point & end, bool ends_in_elem) override;

protected:
  /// The optic study
  const OpticRayStudyBase & _optic_study;

  // The number of energy groups in the optic study
  const unsigned int _num_energy_groups;

  /// The first value in the optic data for the current Ray
  Real * const & _optic_data;
};
