#pragma once

#include "OpticRayStudyBase.h"

class UniformOpticRayStudy : public OpticRayStudyBase
{
public:
  UniformOpticRayStudy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  void defineRays() override;
  void buildBoundingBoxes();

  const std::vector<BoundaryID> _bnd_ids;
  const Point _direction;
  const std::vector<Real> _intensity;
  std::vector<std::vector<Real>> _splits;

  std::vector<BoundingBox> _local_bboxes;
  std::vector<BoundingBox> _bboxes;
  std::vector<std::vector<Real>> _widths;
  std::vector<std::vector<unsigned int>> _dims;
};
