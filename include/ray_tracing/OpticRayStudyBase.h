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

  void generateRays() override;

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
  RayDataIndex energyDataIndex() const { return _energy_data_index; }

  /**
   * Adds a ray to be traced while tracing has already begun, specifically for use in adding
   * reflected rays while also tracing
   * @param start The starting point
   * @param elem The starting elem
   * @param direction The direction of travel
   * @param tid The thread
   * /return The created Ray (for modification: changing data, etc)
   */
  Ray & addReflectedRay(const Point & start,
                        const Elem * elem,
                        const Point & direction,
                        const Ray & parent_ray,
                        const Real energy_factor,
                        const THREAD_ID tid);

protected:
  virtual void defineRays() override = 0;

  void generateThreadedNextIDs();

  /**
   * Defines a Ray for generation
   *
   * Users should call this in their OpticRayStudyBase derived classes to define individual rays
   *
   * @param start The start point
   * @param direction The direction
   * @param elem The element the Ray starts in (optional, claiming will find this if not given)
   * @param incoming_side The side on elem the Ray starts on (optional, claiming will find this if
   * not given)
   * @param size_data Whether or not to the size the Ray data now (optional, default: true)
   * @param size_aux_data Whether or not to size the Ray aux data now (optional, default: true)
   * /return The created Ray (for modification: changing data, etc)
   */
  Ray & defineRay(const Point & start,
                  const Point & direction,
                  const RayID id,
                  const Elem * elem = nullptr,
                  const unsigned short incoming_side = RayTracingCommon::invalid_side,
                  const bool size_data = true,
                  const bool size_aux_data = true);

  /// Number of energy groups
  const unsigned int _num_energy_groups;

  /// Name of the ray variable that holds the optic, group-wise data (if any)
  const RayDataIndex _energy_data_index;

  std::vector<RayID> _threaded_next_ids;
  std::vector<RayID> _threaded_rays_added;
  RayID _max_thread_ids;
};
