#pragma once

#include "OpticRayStudyBase.h"

class SimpleOpticRayStudy : public OpticRayStudyBase
{
public:
  SimpleOpticRayStudy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void defineRays() override;

  /// The points to start the Rays from
  const std::vector<Point> & _start_points;
  /// The Ray directions
  const std::vector<Point> & _directions;
  /// The data to fill on each Ray (optional)
  const std::vector<Real> * _data;

  /// The next available ID to assign to a Ray
  dof_id_type _next_id;
};
