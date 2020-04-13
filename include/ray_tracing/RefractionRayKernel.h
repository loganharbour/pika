#pragma once

#include "RayKernelBase.h"

#include "OnTheFlyVariableValue.h"

/**
 * RayKernel that physically refracts a Ray at the point a phase hits a
 * threadhold within an element.
 */
class RefractionRayKernel : public RayKernelBase
{
public:
  RefractionRayKernel(const InputParameters & params);

  static InputParameters validParams();

  virtual void
  onSegment(const Elem * elem, const Point & start, const Point & end, bool ends_in_elem) override;
  virtual void preTrace(const std::shared_ptr<Ray> & ray) override;

protected:
  /**
   * The direction at which the wave refracts when passing through a boundary
   * between two different isotropic media, as described by Snell's law
   * @param direction The incoming direction
   * @param normal The normal of the boundary between the two media
   * @param r1 The index of refraction for the incoming medium
   * @param r2 The index of refraction for the outgoing medium
   * @return The refracted direction
   */
  Point snell(const Point & direction, const Point & normal, const Real r1, const Real r2) const;

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
  Point phaseChange(const Point & start,
                    const Point & end,
                    const Point & direction,
                    const Real & start_phase,
                    const Real & end_phase);

  /// The field variable that contains the phase
  const VariableValue & _phase;
  /// The gradient of the field variable that contains the phase
  const VariableGradient & _grad_phase;
  /// On the fly value for the phase for quick reinits
  OnTheFlyVariableValue _otf_phase;

  /// The field variable that contains the refractive index
  const VariableValue & _refractive_index;

  /// The physical location of the segment's quadrature points
  const MooseArray<Point> & _q_point;

  /// Whether or not the current Ray refracted on the last segment
  bool _just_refracted;
  /// The phase threshold for refraction
  const Real _threshold;
};
