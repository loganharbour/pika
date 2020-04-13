#include <functional>

#include "MooseTypes.h"

namespace BasicLineSearch
{

/**
 * Brent's root-finding algorithm.
 *
 * Combines the bisection method, the secant method, and inverse quadratic
 * interpolation.
 *
 * @param a0 The lower bounds
 * @param b0 The upper bounds
 * @param f_a0 The function evaluated at a0
 * @param f_b0 The function evaluated at b0
 * @param f Function for evaluating the function
 * @return The root
 */
Real brents(const Real a0,
            const Real b0,
            const Real f_a0,
            const Real f_b0,
            const std::function<Real(const Real)> & f);

void bracketedError(const Real a0, const Real b0, const Real f_a0, const Real f_b0);

} // namespace BasicLineSearch
