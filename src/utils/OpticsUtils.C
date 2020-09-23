#include "OpticsUtils.h"

#include "MathUtils.h"
#include "MooseUtils.h"

namespace OpticsUtils
{

SnellInteractionType
snell(const Point & direction,
      const Point & normal,
      const Real n1,
      const Real n2,
      Point & refracted_direction)
{
  mooseAssert(n1 > 0 && n2 > 0, "Invalid index of refraction");
  mooseAssert(MooseUtils::absoluteFuzzyEqual(direction.norm(), 1), "Unnormalized direction");
  mooseAssert(MooseUtils::absoluteFuzzyEqual(normal.norm(), 1), "Unnormalized normal");

  // Dot product between the normal and the incident direction
  const auto c = std::abs(direction * normal);

  // Normal and direction are parallel, continues in the same direction
  if (MooseUtils::absoluteFuzzyEqual(c, 1.0))
  {
    refracted_direction = direction;
    return SNELL_PARALLEL;
  }

  const auto r = n1 / n2;

  // Refracts: n2 > n1, or n2 < n1 and below the critical angle
  if (r <= 1 || (1 - r * r * (1 - c * c)) >= 0)
  {
    refracted_direction = r * direction;
    refracted_direction += (r * c - std::sqrt(1 - (r * r) * (1 - c * c))) * normal;
    refracted_direction = refracted_direction.unit();
    return SNELL_REFRACT;
  }

  // Otherwise, reflect
  return SNELL_INTERNAL_REFLECTION;
}

Real
reflectionCoefficient(const Point & incident_direction,
                      const Point & reflected_direction,
                      const Point & normal,
                      const Real n1,
                      const Real n2)
{
  mooseAssert(n1 > 0 && n2 > 0, "Invalid index of refraction");
  mooseAssert(MooseUtils::absoluteFuzzyEqual(incident_direction.norm(), 1),
              "Unnormalized incident direction");
  mooseAssert(MooseUtils::absoluteFuzzyEqual(reflected_direction.norm(), 1),
              "Unnormalized reflected direction");
  mooseAssert(MooseUtils::absoluteFuzzyEqual(normal.norm(), 1), "Unnormalized normal");

  const auto ci = std::abs(incident_direction * normal);
  const auto ct = std::abs(reflected_direction * normal);
  const auto Rs = MathUtils::pow((n1 * ci - n2 * ct) / (n1 * ci + n2 * ct), 2);
  const auto Rp = MathUtils::pow((n1 * ct - n2 * ci) / (n1 * ct + n2 * ci), 2);
  return 0.5 * (Rs + Rp);
}

Point
reflectedDirection(const Point & direction, const Point & normal)
{
  mooseAssert(MooseUtils::absoluteFuzzyEqual(direction.norm(), 1), "Unnormalized direction");
  mooseAssert(MooseUtils::absoluteFuzzyEqual(normal.norm(), 1), "Unnormalized normal");

  Point reflected_direction = direction;
  reflected_direction -= 2.0 * (reflected_direction * normal) * normal;
  return reflected_direction;
}

};
