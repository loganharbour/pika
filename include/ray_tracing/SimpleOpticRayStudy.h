#pragma once

#include "OpticRayStudyBase.h"

class SimpleOpticRayStudy : public OpticRayStudyBase
{
public:
  SimpleOpticRayStudy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  void defineRays() override;

  /// The points to start the Rays from
  const std::vector<Point> & _start_points;
  /// The Ray directions
  const std::vector<Point> & _directions;
  /// The data to fill on each Ray (optional)
  const std::vector<Real> * _data;
};
