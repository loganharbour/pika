#include "MooseTypes.h"

namespace OpticsUtils
{

/**
 * The interaction types for Snell's law:
 * - Internal reflection: When refraction does not occur; happens when light travels from a medium
 * with ahigher refractive index to one with a lower refractive index and the incident angle is >
 * the critical angle
 * - Parallel: When the incident direction and the normal are parallel, so the refracted/reflected
 * directions are the same
 * - Refract: Standard refraction
 */
enum SnellInteractionType
{
  SNELL_INTERNAL_REFLECTION = 0,
  SNELL_PARALLEL = 1,
  SNELL_REFRACT = 2
};

/**
 * The direction at which the wave refracts when passing through a boundary between two different
 * isotropic media, as described by Snell's law
 * @param direction The incident direction with the surface
 * @param normal The normal of the boundary between the two media
 * @param n1 The index of refraction for the incoming medium
 * @param n2 The index of refraction for the outgoing medium
 * @param refracted_direction Filled with the refracted direction, if any
 * @return The interaction type
 */
SnellInteractionType snell(const Point & direction,
                           const Point & normal,
                           const Real n1,
                           const Real n2,
                           Point & refracted_direction);

/**
 * Approximation for the reflection coefficient, that is, how much of a wave is reflected
 * by an impedance discontinuity in the transmission medium.
 * @param incident_direction The incident direction with the surface
 * @param reflected_direction The reflected direction
 @param normal The normal of the boundary between the two media
 * @param n1 The index of refraction for the incoming medium
 * @param n2 The index of refraction for the outgoing medium
 * @return The approximate reflection coefficient
 */
Real reflectionCoefficient(const Point & incident_direction,
                           const Point & reflected_direction,
                           const Point & normal,
                           const Real n1,
                           const Real n2);

/**
 * Computes the reflected direction given a direction and an outward normal for the surface it
 * reflects off of.
 * @param direction The incident direction
 * @param normal The outward normal
 * @return The reflected direction
 */

Point reflectedDirection(const Point & direction, const Point & normal);

}
