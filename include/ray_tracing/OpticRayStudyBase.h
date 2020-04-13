#pragma once

#include "RepeatableRayStudyBase.h"

/**
 * Base optic ray study that overrides the segment quadrature build so that
 * reinits on Ray segments are done at the beginning and the end of said
 * segments.
 */
class OpticRayStudyBase : public RepeatableRayStudyBase
{
public:
  OpticRayStudyBase(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void buildSegmentQuadrature(const Point & start,
                                      const Point & end,
                                      std::vector<Point> & points,
                                      std::vector<Real> & weights) const override;

  /**
   * Number of energy groups
   */
  unsigned int numEnergyGroups() const { return _num_energy_groups; }

  /**
   * The name of the ray variable for the optic group-wise data
   */
  std::string opticRayVariableName() const { return _optic_ray_variable_name; }

protected:
  virtual void defineRays() override = 0;

  /// Number of energy groups
  const unsigned int _num_energy_groups;

  /// Name of the ray variable that holds the optic, group-wise data (if any)
  const std::string _optic_ray_variable_name;
};
